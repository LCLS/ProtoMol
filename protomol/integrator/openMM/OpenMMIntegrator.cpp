#include <protomol/integrator/openMM/OpenMMIntegrator.h>
#include <protomol/base/Report.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/type/Vector3DBlock.h>
#include <protomol/force/ForceGroup.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/topology/TopologyUtilities.h>
#include <protomol/ProtoMolApp.h>
#include <protomol/base/PMConstants.h>
#include <protomol/base/Zap.h>
#include <protomol/topology/LennardJonesParameters.h>

#include <openmm/LocalEnergyMinimizer.h>

#include <vector>
#include <algorithm>

#include <iostream>

using namespace std;
using namespace ProtoMol::Report;
using namespace ProtoMol;

const string OpenMMIntegrator::keyword( "OpenMM" );

OpenMMIntegrator::OpenMMIntegrator() : STSIntegrator(), isLTMD( false ) {
	system = 0;
	integrator = 0;
	context = 0;
}

OpenMMIntegrator::OpenMMIntegrator( const std::vector<Value>& params, ForceGroup *overloadedForces )
	: STSIntegrator( params[0], overloadedForces ), isLTMD( false ) {

	mTemperature = params[1];
	mGamma = params[2];
	mSeed = params[3];

	isUsingHarmonicBondForce = params[4];
	isUsingHarmonicAngleForce = params[5];
	isUsingRBDihedralForce = params[6];
	isUsingPeriodicTorsionForce = params[7];
	isUsingNonBondedForce = params[8];
	isUsingGBForce = params[9];
	isUsingSCPISMForce = params[10];
	isUsingImproperTorsionForce = params[11];
	isUsingUreyBradleyForce = params[12];
	isUsingCHARMMLennardJonesForce = params[13];

	mCommonMotionRate = params[14];
	mGBSAEpsilon = params[15];
	mGBSASolvent = params[16];

	mPlatform = params[17];
	mMinSteps = params[18];
	mTolerance = params[19];

	mNonbondedCutoff = params[20];
	mGBCutoff = params[21];

	mPropagatorDevice = params[22];
	mBlockDevice = params[23];

	system = 0;
	integrator = 0;
	context = 0;
}

OpenMMIntegrator::~OpenMMIntegrator() {
	zap( context );
	zap( integrator );
	zap( system );
}

struct NBForce {
	int atom1, atom2;
	Real charge, sigma, epsilon;

	NBForce( int a, int b, Real c, Real s, Real e )
		: atom1( a ), atom2( b ), charge( c ), sigma( s ), epsilon( e ) {
	}

	bool operator< ( const NBForce &other ) const {
		if( atom1 < other.atom1 ) return true;
		if( atom1 == other.atom1 && atom2 < other.atom2 ) return true;

		return false;
	}
};

void OpenMMIntegrator::initialize( ProtoMolApp *app ) {
	STSIntegrator::initialize( app );
	initializeForces();

	std::string sDirectory;
	std::stringstream stream( OpenMM::Platform::getDefaultPluginsDirectory() );

	while( getline( stream, sDirectory, ':' ) ){
		std::cout << "Loading OpenMM plugins from directory: " << sDirectory << std::endl;
		OpenMM::Platform::loadPluginsFromDirectory( sDirectory );
	}

	//find system size
	unsigned int sz = app->positions.size();

	//find constraint size
	const std::vector<Bond::Constraint> *myListOfConstraints =
		&( app->topology->bondRattleShakeConstraints );

	unsigned int numConstraints = ( *myListOfConstraints ).size();

	//Initialize system
	system = new OpenMM::System();
	for( unsigned int i = 0; i < sz; ++i ) {
		system->addParticle( app->topology->atoms[i].scaledMass );
	}

	//remove common motion?
	if( mCommonMotionRate > 0 ) {
		mForceList.push_back( "CenterOfMass" );
		system->addForce( new OpenMM::CMMotionRemover( mCommonMotionRate ) );
	}

	//openMM forces
	if( isUsingHarmonicBondForce ) {
		mForceList.push_back( "Bond" );
		unsigned int numBonds = app->topology->bonds.size();

		unsigned int numConstBonds = 0;

		if( numConstraints ) {
			for( unsigned int i = 0; i < numBonds; ++i ) {
				if( ( app->topology->atoms[ app->topology->bonds[i].atom1 ].name[0] == 'H' ) ||
						( app->topology->atoms[ app->topology->bonds[i].atom2 ].name[0] == 'H' ) ) {
					numConstBonds++;
				}
			}
		}

		OpenMM::HarmonicBondForce *bonds = new OpenMM::HarmonicBondForce();
		system->addForce( bonds );

		for( unsigned int i = 0; i < numBonds; ++i ) {
			unsigned int a1 = app->topology->bonds[i].atom1, a2 = app->topology->bonds[i].atom2;

			Real r_0 = app->topology->bonds[i].restLength  * Constant::ANGSTROM_NM;
			Real k = app->topology->bonds[i].springConstant
					 * Constant::KCAL_KJ * Constant::INV_ANGSTROM_NM * Constant::INV_ANGSTROM_NM * 2.0; //times 2 as Amber is 1/2 k(b-b_0)^2
			if( numConstraints ) {

			} else {
				bonds->addBond( a1, a2, r_0, k );
			}
		}
	}

	if( isUsingHarmonicAngleForce ) {
		mForceList.push_back( "Angle" );
		unsigned int numAngles = app->topology->angles.size();

		OpenMM::HarmonicAngleForce *angles = new OpenMM::HarmonicAngleForce();
		system->addForce( angles );
		for( unsigned int i = 0; i < numAngles; i++ ) {
			unsigned int a1 = app->topology->angles[i].atom1;
			unsigned int a2 = app->topology->angles[i].atom2;
			unsigned int a3 = app->topology->angles[i].atom3;
			Real theta0 = acos( cos( app->topology->angles[i].restAngle ) );
			Real k_t = app->topology->angles[i].forceConstant * Constant::KCAL_KJ * 2.0; //times 2 as Amber is 1/2 k(a-a_0)^2

			angles->addAngle( a1, a2, a3, theta0, k_t );
		}
	}

	if( isUsingPeriodicTorsionForce ) {
		mForceList.push_back( "Dihedral" );
		unsigned int numPTor = app->topology->dihedrals.size();

		OpenMM::PeriodicTorsionForce *PTorsion = new OpenMM::PeriodicTorsionForce();
		system->addForce( PTorsion );
		for( unsigned int i = 0; i < numPTor; i++ ) {
			unsigned int a1 = app->topology->dihedrals[i].atom1;
			unsigned int a2 = app->topology->dihedrals[i].atom2;
			unsigned int a3 = app->topology->dihedrals[i].atom3;
			unsigned int a4 = app->topology->dihedrals[i].atom4;

			for( int j = 0; j < app->topology->dihedrals[i].multiplicity; j++ ) {
				unsigned int mult = app->topology->dihedrals[i].periodicity[j];
				Real phiA = app->topology->dihedrals[i].phaseShift[j];
				Real cpA = app->topology->dihedrals[i].forceConstant[j] * Constant::KCAL_KJ;

				PTorsion->addTorsion( a1, a2, a3, a4, mult, phiA, cpA );
			}
		}
	}

	if( isUsingRBDihedralForce ) {
		mForceList.push_back( "Improper" );
		unsigned int numRBDih = app->topology->rb_dihedrals.size();

		OpenMM::RBTorsionForce *RBDihedral = new OpenMM::RBTorsionForce();

		system->addForce( RBDihedral );
		for( unsigned int i = 0; i < numRBDih; i++ ) {
			unsigned int a1 = app->topology->rb_dihedrals[i].atom1;
			unsigned int a2 = app->topology->rb_dihedrals[i].atom2;
			unsigned int a3 = app->topology->rb_dihedrals[i].atom3;
			unsigned int a4 = app->topology->rb_dihedrals[i].atom4;
			Real C0 = app->topology->rb_dihedrals[i].C0 * Constant::KCAL_KJ;
			Real C1 = app->topology->rb_dihedrals[i].C1 * Constant::KCAL_KJ;
			Real C2 = app->topology->rb_dihedrals[i].C2 * Constant::KCAL_KJ;
			Real C3 = app->topology->rb_dihedrals[i].C3 * Constant::KCAL_KJ;
			Real C4 = app->topology->rb_dihedrals[i].C4 * Constant::KCAL_KJ;
			Real C5 = app->topology->rb_dihedrals[i].C5 * Constant::KCAL_KJ;
			if( C0 != 0 || C1 != 0 || C2 != 0 || C3 != 0 || C4 != 0 || C5 != 0 ) {
				RBDihedral->addTorsion( a1, a2, a3, a4, C0, C1, C2, C3, C4, C5 );
			}
		}
	}

	if( isUsingNonBondedForce) {
		mForceList.push_back( "Nonbonded" );

		cout << "Nonbonded" << endl;

		//get 1-4 interaction size
		unsigned int exclSz = app->topology->exclusions.getTable().size();

		OpenMM::NonbondedForce *nonbonded = new OpenMM::NonbondedForce();
		system->addForce( nonbonded );

		//normal interactions
		for( unsigned int i = 0; i < sz; i++ ) {
			int type1 = app->topology->atoms[i].type;
			Real sigma = app->topology->atomTypes[type1].sigma;
			Real epsilon = app->topology->atomTypes[type1].epsilon;
			Real charge = app->topology->atoms[i].scaledCharge / Constant::SQRTCOULOMBCONSTANT;

			nonbonded->addParticle( charge, sigma * Constant::ANGSTROM_NM , epsilon * Constant::KCAL_KJ );
		}

		//1-4 interactions, note fudgeQQ set to 0.6059 and fudgeLJ set to 0.33333.
		for( unsigned int i = 0; i < exclSz; i++ ) {

			unsigned int atom1 = ( app->topology->exclusions.getTable() )[i].a1;
			unsigned int atom2 = ( app->topology->exclusions.getTable() )[i].a2;

			//full exclusion?
			if( ( app->topology->exclusions.getTable() )[i].excl == EXCLUSION_FULL ) {
				nonbonded->addException( atom1, atom2, 0.0, 1.0, 0.0 );
			}

			//modified exclusion?
			if( ( app->topology->exclusions.getTable() )[i].excl == EXCLUSION_MODIFIED ) {

				unsigned int type1 = app->topology->atoms[atom1].type;
				unsigned int type2 = app->topology->atoms[atom2].type;
				Real sigma = 0.5 * ( app->topology->atomTypes[type1].sigma +
									 app->topology->atomTypes[type2].sigma );
				Real epsilon = app->topology->LJScalingFactor *
							   sqrt( app->topology->atomTypes[type1].epsilon *
									 app->topology->atomTypes[type2].epsilon );
				Real chargeij =  app->topology->coulombScalingFactor * //FudgeQQ
								 ( app->topology->atoms[atom1].scaledCharge / Constant::SQRTCOULOMBCONSTANT ) *
								 ( app->topology->atoms[atom2].scaledCharge / Constant::SQRTCOULOMBCONSTANT );

				nonbonded->addException( atom1, atom2, chargeij, sigma * Constant::ANGSTROM_NM, epsilon * Constant::KCAL_KJ );

			}
		}

		if( mNonbondedCutoff != 0 ){
			nonbonded->setNonbondedMethod( OpenMM::NonbondedForce::CutoffNonPeriodic );
			nonbonded->setCutoffDistance( mNonbondedCutoff * Constant::ANGSTROM_NM );
		}
	}

	// Add GBSA if needed.
	if( app->topology->implicitSolvent  == GBSA && isUsingGBForce ) {
		mForceList.push_back( "GBSA" );
		OpenMM::GBSAOBCForce *gbsa = new OpenMM::GBSAOBCForce();
		system->addForce( gbsa );

		gbsa->setSoluteDielectric( mGBSAEpsilon );
		gbsa->setSolventDielectric( mGBSASolvent );

		vector<Real> scaleFactors;
		getObcScaleFactors( scaleFactors );

		for( unsigned int i = 0; i < sz; ++i ) {
			Real charge = app->topology->atoms[i].scaledCharge / Constant::SQRTCOULOMBCONSTANT;
			Real radius = app->topology->atoms[i].myGBSA_T->vanDerWaalRadius * Constant::ANGSTROM_NM; //0.1 factor in openMM, file in A

			gbsa->addParticle( charge, radius, scaleFactors[i] );
		}

		if( mGBCutoff != 0 ){
			gbsa->setNonbondedMethod( OpenMM::GBSAOBCForce::CutoffNonPeriodic );
			gbsa->setCutoffDistance( mGBCutoff * Constant::ANGSTROM_NM );
		}
	}

	// Add SCPISM if needed
	if (app->topology->doSCPISM && isUsingSCPISMForce) {
	  cout << "Using SCPISM" << endl;
	  OpenMM::CustomGBForce *custom = new OpenMM::CustomGBForce();
	  system->addForce(custom);

	  // global parameter names
	  custom->addGlobalParameter("epsilon", 80.0); // hard-coded in ProtoMol
	  custom->addGlobalParameter("pi", 3.14159265); // OpenMM doesn't include pi?
	  custom->addGlobalParameter("nAtoms", sz); // number of atoms
	  custom->addGlobalParameter("E", 0.80); // hard-coded in ProtoMol
	  custom->addGlobalParameter("KCAL_KJ", Constant::KCAL_KJ);

	  // per particle parameter names
	  custom->addPerParticleParameter("alpha");
	  custom->addPerParticleParameter("eta");
	  custom->addPerParticleParameter("C");
	  custom->addPerParticleParameter("q");
	  custom->addPerParticleParameter("zeta");
	  custom->addPerParticleParameter("HBondType");
	  custom->addPerParticleParameter("residue");
	  custom->addPerParticleParameter("g");
	  custom->addPerParticleParameter("isHN");
	  custom->addPerParticleParameter("isO");
	  custom->addPerParticleParameter("atomIndex");

	  //	  custom->setCutoffDistance(0.5); // 0.5 NM, 5 A
	  //custom->setNonbondedMethod(OpenMM::CustomGBForce::CutoffNonPeriodic);

	  // Per-atom or per-atom-pair computed values -- used to precompute values for energy term (e.g., Born radii)
	  // Born Radii
	  custom->addComputedValue("etasum", "eta1 * f * exp(-C1 * rScaled) + polar_term;"
				   // only add polar term if atoms are H-bonded and not in the same residue
				   // atom1 must have H-bond of type PH (1) and atom 2 must have
				   // H-bond type of PA (2)
				   "polar_term = step(HBondType1 - 1) * step(1 - HBondType1) * "
				   "step(HBondType2 - 2) * step(2 - HBondType2) *"
				   "(1 - step(residue1 - residue2) * step(residue2 - residue1) ) * gsum;"
				   "gsum = g_i * g2 * f * exp(-E * rScaled);"
				   // g_i = -0.378 if atom1 == HN and atom2 == O
				   // otherwise stick with given value
				   "g_i = g1 + isHN1 * isO2 * (-0.378 - g1);"
				   "f = step(5 - rScaled) * (1 - rScaled^4/625)^4;"
				   // convert distance from NM to A
				   // all parameters assume A
				   "rScaled = 10 * r;",
				   OpenMM::CustomGBForce::ParticlePair);
	  custom->addComputedValue("R", "zeta + etasum;",
	  			   OpenMM::CustomGBForce::SingleParticle);
	  // Energy terms
	  // Born Self
	  custom->addEnergyTerm("KCAL_KJ * 0.5 * q^2 * (1 / D_s - 1) / R;"
				"D_s = (1 + epsilon) / (1 + K * exp( - alpha * R)) - 1;"
				"K = 0.5 * (epsilon - 1);",
				OpenMM::CustomGBForce::SingleParticle);
	  // SCPISM Coulomb force
	  custom->addEnergyTerm("KCAL_KJ * q1 * q2 * scalingFactors(atomIndex1 * nAtoms + atomIndex2) / (rScaled * D_s);"
				"D_s = (1 + epsilon) / (1 + K * exp( - sqrt(alpha1) * sqrt(alpha2) * rScaled)) - 1;"
				"K = 0.5 * (epsilon - 1);"
				// all parameters assume A,
			        // convert distance from A to NM
				"rScaled = 10 * r;",
				OpenMM::CustomGBForce::ParticlePair);

	  for(unsigned int i = 0; i < sz; i++)
	    {
	      vector<double> parameters(11);
	      int type = app->topology->atoms[i].type;
	      parameters[0] = app->topology->atomTypes[type].mySCPISM_T->alpha;
	      parameters[1] = app->topology->atoms[i].mySCPISM_A->eta;
	      parameters[2] = app->topology->atomTypes[type].mySCPISM_T->C_i;
	      parameters[3] = app->topology->atoms[i].scaledCharge;
	      parameters[4] = app->topology->atoms[i].mySCPISM_A->zeta;
	      parameters[5] = app->topology->atomTypes[type].mySCPISM_T->isHbonded;
	      parameters[6] = app->topology->atoms[i].residue_seq;
	      parameters[7] = app->topology->atomTypes[type].mySCPISM_T->g_i;
	      parameters[8] = app->topology->atoms[i].name == "HN" ? 1.0 : 0.0;
	      parameters[9] = app->topology->atoms[i].name == "O" ? 1.0 : 0.0;
	      parameters[10] = i;

	      custom->addParticle(parameters);
	    }

	  // default scaling factor is 1.0 (no scaling)
	  // CustomGBForce doesn't provide a way to control
	  // modified interactions through the API
	  // Thus, we're (ab)using the tabulated function to
	  // provide similar functionality.
	  vector<double> scalingFactors(sz * sz, 1.0);

	  unsigned int exclSz = app->topology->exclusions.getTable().size();
	  for( unsigned int i = 0; i < exclSz; i++ ) {

	    const unsigned int atom1 = ( app->topology->exclusions.getTable() )[i].a1;
	    const unsigned int atom2 = ( app->topology->exclusions.getTable() )[i].a2;

	    const ExclusionClass excl = app->topology->exclusions.check(atom1, atom2);

	    if( excl == EXCLUSION_FULL ) {
	      custom->addExclusion( atom1, atom2);
	    }
	    else if ( excl == EXCLUSION_MODIFIED ) {
	      scalingFactors[atom1 * sz + atom2] = app->topology->coulombScalingFactor;
	      scalingFactors[atom2 * sz + atom1] = app->topology->coulombScalingFactor;
	    }
	  }

	  custom->addFunction("scalingFactors", scalingFactors, 0, sz * sz - 1);

	}

	// LJ as a custom force to go with SCPISM
	if( isUsingCHARMMLennardJonesForce) {
	  std::cout << "Performing Custom LJ" << std::endl;
	  mForceList.push_back( "Nonbonded" );

	  //get 1-4 interaction size
	  unsigned int exclSz = app->topology->exclusions.getTable().size();

	  string energy = "A / rScaled^12 - B / rScaled^6;"
	    // scale NM to A to match parameter units
	    "rScaled = 10.0 * r;"
	    "A = aTable(atomIndex1 * nAtoms + atomIndex2);"
	    "B = bTable(atomIndex1 * nAtoms + atomIndex2);";
	  OpenMM::CustomNonbondedForce *nonbonded = new OpenMM::CustomNonbondedForce(energy);
	  system->addForce( nonbonded );

	  nonbonded->addGlobalParameter("nAtoms", sz);

	  nonbonded->addPerParticleParameter("atomIndex");

	  // Store A and B parameters as tables
	  vector<double> aTable(sz * sz, 0.0);
	  vector<double> bTable(sz * sz, 0.0);

	  for(unsigned int i = 0; i < sz; i++)
	    {
	      vector<double> parameters(1);

	      const unsigned int atomType1 = app->topology->atoms[i].type;

	      parameters[0] = i;

	      nonbonded->addParticle(parameters);

	      for(unsigned int j = i+1; j < sz; j++)
		{
		  const unsigned int atomType2 = app->topology->atoms[j].type;
		  const LennardJonesParameters ljParams =
		    app->topology->lennardJonesParameters(atomType1, atomType2);

		  const unsigned int atomIndex1 = i;
		  const unsigned int atomIndex2 = j;

		  aTable[atomIndex1 * sz + atomIndex2] = ljParams.A * Constant::KCAL_KJ;
		  aTable[atomIndex2 * sz + atomIndex1] = ljParams.A * Constant::KCAL_KJ;
		  bTable[atomIndex1 * sz + atomIndex2] = ljParams.B * Constant::KCAL_KJ;
		  bTable[atomIndex2 * sz + atomIndex1] = ljParams.B * Constant::KCAL_KJ;
		}
	    }

	  for( unsigned int i = 0; i < exclSz; i++ )
	    {


	      const unsigned int atomIndex1 = ( app->topology->exclusions.getTable() )[i].a1;
	      const unsigned int atomIndex2 = ( app->topology->exclusions.getTable() )[i].a2;

	      const ExclusionClass excl = app->topology->exclusions.check(atomIndex1, atomIndex2);

	      if(excl == EXCLUSION_FULL)
		{
		  nonbonded->addExclusion(atomIndex1, atomIndex2);
		  aTable[atomIndex1 * sz + atomIndex2] = 0.0;
		  aTable[atomIndex2 * sz + atomIndex1] = 0.0;
		  bTable[atomIndex1 * sz + atomIndex2] = 0.0;
		  bTable[atomIndex2 * sz + atomIndex1] = 0.0;
		}
	      else if(excl == EXCLUSION_MODIFIED)
		{
		  const unsigned int atomType1 = app->topology->atoms[atomIndex1].type;
		  const unsigned int atomType2 = app->topology->atoms[atomIndex2].type;

		  const LennardJonesParameters ljParams =
		    app->topology->lennardJonesParameters(atomType1, atomType2);

		  aTable[atomIndex1 * sz + atomIndex2] = ljParams.A14 * Constant::KCAL_KJ;
		  aTable[atomIndex2 * sz + atomIndex1] = ljParams.A14 * Constant::KCAL_KJ;
		  bTable[atomIndex1 * sz + atomIndex2] = ljParams.B14 * Constant::KCAL_KJ;
		  bTable[atomIndex2 * sz + atomIndex1] = ljParams.B14 * Constant::KCAL_KJ;
		}
	    }
	  nonbonded->addFunction("aTable", aTable, 0, sz * sz - 1);
	  nonbonded->addFunction("bTable", bTable, 0, sz * sz - 1);
	}

	if(isUsingUreyBradleyForce)
	  {
	    string energy = "forceConstant * (r - restLength) * (r - restLength);";
	    OpenMM::CustomBondForce *bondForce = new OpenMM::CustomBondForce(energy);
	    system->addForce(bondForce);

	    bondForce->addPerBondParameter("forceConstant");
	    bondForce->addPerBondParameter("restLength");

	    for(unsigned int i = 0; i < app->topology->angles.size(); i++)
	      {
		const Angle angle = app->topology->angles[i];

		vector<double> parameters(2);
		parameters[0] = angle.ureyBradleyConstant * Constant::KCAL_KJ * Constant::NM_ANGSTROM * Constant::NM_ANGSTROM;
		parameters[1] = angle.ureyBradleyRestLength * Constant::ANGSTROM_NM;

		bondForce->addBond(angle.atom1, angle.atom3, parameters);
	      }
	  }

	if(isUsingImproperTorsionForce && app->topology->impropers.size() > 0)
	  {
	    cout << "Impropers" << endl;

	    string energy = "forceConstant * (theta - phaseShift)^2";
	    OpenMM::CustomTorsionForce *torsion = new OpenMM::CustomTorsionForce(energy);
	    system->addForce(torsion);

	    torsion->addPerTorsionParameter("forceConstant");
	    torsion->addPerTorsionParameter("phaseShift");

	    for(unsigned int i = 0; i < app->topology->impropers.size(); i++)
	      {
		const Torsion currTorsion = app->topology->impropers[i];

		vector<double> parameters(2);
		parameters[0] = currTorsion.forceConstant[0] * Constant::KCAL_KJ;
		parameters[1] = currTorsion.phaseShift[1];

		torsion->addTorsion(currTorsion.atom1, currTorsion.atom2, currTorsion.atom3, currTorsion.atom4, parameters);
	      }
	  }


	// Set constraints.
	for( unsigned int i = 0; i < numConstraints; ++i ) {

		int atom1 = ( *myListOfConstraints )[i].atom1;
		int atom2 = ( *myListOfConstraints )[i].atom2;
		Real restLength = ( *myListOfConstraints )[i].restLength * Constant::ANGSTROM_NM;

		system->addConstraint( atom1, atom2, restLength );
	}

#ifdef HAVE_OPENMM_FBM
		for( unsigned int i = 0; i < mForceList.size(); i++ ){
		  cout << "Adding force " << mForceList[i] << " to list." << endl;
			mFBMParameters.forces.push_back( FBMForce( mForceList[i], i ) );
		}

#endif

	if( isLTMD ){
#ifdef HAVE_OPENMM_LTMD
		for( unsigned int i = 0; i < mForceList.size(); i++ ){
			mLTMDParameters.forces.push_back( OpenMM::LTMD::Force( mForceList[i], i ) );
		}

		OpenMM::LTMD::Integrator *ltmd = new OpenMM::LTMD::Integrator( mTemperature, mGamma, getTimestep() * Constant::FS_PS, mLTMDParameters );
		ltmd->setRandomNumberSeed( mSeed );

		integrator = ltmd;
#endif
	}else{
		OpenMM::LangevinIntegrator *langevin = new OpenMM::LangevinIntegrator( mTemperature, mGamma, getTimestep() * Constant::FS_PS );
		langevin->setRandomNumberSeed( mSeed );

		integrator = langevin;
	}

	std::string sPlatform = "Reference";

	switch( mPlatform ){
		case 0:
			sPlatform = "Reference";
			break;
		case 1:
			sPlatform = "OpenCL";
			break;
		case 2:
			sPlatform = "CUDA";
			break;
	}

	std::cout << "OpenMM Propagation Platform: " << sPlatform << std::endl;
	OpenMM::Platform& platform = OpenMM::Platform::getPlatformByName( sPlatform );

	if( mPlatform != 0 && mPropagatorDevice != -1 ){
        std::ostringstream stream;
        stream << mPropagatorDevice;

        std::cout << "OpenMM Propagation Device: " << mPropagatorDevice << std::endl;

        if( mPlatform == 1 ){
            platform.setPropertyDefaultValue("OpenCLDeviceIndex", stream.str() );
        }else{
            platform.setPropertyDefaultValue("CudaDeviceIndex", stream.str() );
        }
    }

	if( mPlatform == 2 ){
		//platform.setPropertyDefaultValue("CudaPrecision", "double");
	}

	std::cout << "Creating context" << std::endl;
	context = new OpenMM::Context( *system, *integrator, platform );

	std::cout << "Initializing positions and velocities" << std::endl;
	std::vector<OpenMM::Vec3> positions, velocities;
	positions.reserve( sz );
	velocities.reserve( sz );

	OpenMM::Vec3 openMMvecp, openMMvecv;
	for( unsigned int i = 0; i < sz; ++i ) {
		for( int j = 0; j < 3; j++ ) {
			openMMvecp[j] = app->positions[i].c[j] * Constant::ANGSTROM_NM;
			openMMvecv[j] = app->velocities[i].c[j] * Constant::ANGSTROM_NM
							* Constant::INV_TIMEFACTOR * Constant::PS_FS;
		}
		positions.push_back( openMMvecp );
		velocities.push_back( openMMvecv );
	}

	context->setPositions( positions );
	context->setVelocities( velocities );

	//print platform
	report << plain << "OpenMM platform is: '" << context->getPlatform().getName()
		   << "' Integrator " << (int)isLTMD << "." << endr;
	const OpenMM::State state = context->getState( OpenMM::State::Positions |
								OpenMM::State::Velocities |
								OpenMM::State::Forces |
								OpenMM::State::Energy );

	if( mMinSteps > 0 ) {
		OpenMM::LocalEnergyMinimizer lem;
		lem.minimize( *context, mTolerance, mMinSteps );
	}
}

long OpenMMIntegrator::run( const long numTimesteps ) {
	preStepModify();

    integrator->step( numTimesteps );

	// Retrive data
	const OpenMM::State state = context->getState( OpenMM::State::Positions |
								OpenMM::State::Velocities |
								OpenMM::State::Forces |
								OpenMM::State::Energy );

	const std::vector<OpenMM::Vec3> positions( state.getPositions() );
	const std::vector<OpenMM::Vec3> velocities( state.getVelocities() );
	const std::vector<OpenMM::Vec3> forces( state.getForces() );

	const unsigned int sz = app->positions.size();
	for( unsigned int i = 0; i < sz; ++i ) {
		for( int j = 0; j < 3; j++ ) {
			app->positions[i].c[j] = positions[i][j] * Constant::NM_ANGSTROM; //nm to A
			app->velocities[i].c[j] = velocities[i][j] * Constant::NM_ANGSTROM *
									  Constant::TIMEFACTOR * Constant::FS_PS; //nm/ps to A/fs?
			( *myForces )[i].c[j] = forces[i][j] * Constant::INV_NM_ANGSTROM * Constant::KJ_KCAL; //KJ/nm to Kcal/A
		}
	}

	app->energies.clear();

	//clear old energies
	app->energies[ScalarStructure::COULOMB] =
		app->energies[ScalarStructure::LENNARDJONES] =
			app->energies[ScalarStructure::BOND] =
				app->energies[ScalarStructure::ANGLE] =
					app->energies[ScalarStructure::DIHEDRAL] =
						app->energies[ScalarStructure::IMPROPER] = 0.0;

	//save total potential energy
	app->energies[ScalarStructure::OTHER] = state.getPotentialEnergy() * Constant::KJ_KCAL;

	//fix time as no forces calculated
	app->topology->time += numTimesteps * getTimestep();

	postStepModify();

    return numTimesteps;
}

void OpenMMIntegrator::getParameters( vector<Parameter> &parameters ) const {
	STSIntegrator::getParameters( parameters );
	parameters.push_back( Parameter( "temperature", Value( mTemperature, ConstraintValueType::NotNegative() ) ) );
	parameters.push_back( Parameter( "gamma", Value( mGamma, ConstraintValueType::NotNegative() ) ) );
	parameters.push_back( Parameter( "seed", Value( mSeed, ConstraintValueType::NotNegative() ), 1234 ) );

	//Force Switches
	parameters.push_back( Parameter( "HarmonicBondForce", Value( isUsingHarmonicBondForce, ConstraintValueType::NoConstraints() ), false ) );
	parameters.push_back( Parameter( "HarmonicAngleForce", Value( isUsingHarmonicAngleForce, ConstraintValueType::NoConstraints() ), false ) );
	parameters.push_back( Parameter( "RBDihedralForce", Value( isUsingRBDihedralForce, ConstraintValueType::NoConstraints() ), false ) );
	parameters.push_back( Parameter( "PeriodicTorsion", Value( isUsingPeriodicTorsionForce, ConstraintValueType::NoConstraints() ), false ) );
	parameters.push_back( Parameter( "NonbondedForce", Value( isUsingNonBondedForce, ConstraintValueType::NoConstraints() ), false ) );
	parameters.push_back( Parameter( "GBForce", Value( isUsingGBForce, ConstraintValueType::NoConstraints() ), false ) );
	parameters.push_back( Parameter( "SCPISMForce", Value ( isUsingSCPISMForce, ConstraintValueType::NoConstraints() ), false) );
	parameters.push_back( Parameter( "ImproperTorsionForce", Value ( isUsingImproperTorsionForce, ConstraintValueType::NoConstraints() ), false) );
	parameters.push_back( Parameter( "UreyBradleyForce", Value ( isUsingUreyBradleyForce, ConstraintValueType::NoConstraints() ), false) );
	parameters.push_back( Parameter( "CHARMMLennardJonesForce", Value ( isUsingCHARMMLennardJonesForce, ConstraintValueType::NoConstraints() ), false) );

	//Implicit solvent parameters
	parameters.push_back( Parameter( "commonmotion", Value( mCommonMotionRate, ConstraintValueType::NotNegative() ), 0.0 ) );
	parameters.push_back( Parameter( "GBSAEpsilon", Value( mGBSAEpsilon, ConstraintValueType::NotNegative() ), 1.0 ) );
	parameters.push_back( Parameter( "GBSASolvent", Value( mGBSASolvent, ConstraintValueType::NotNegative() ), 78.3 ) );

	// OpenMM Parameters
	parameters.push_back( Parameter( "platform", Value( mPlatform, ConstraintValueType::NoConstraints() ), 2 ) );
	parameters.push_back( Parameter( "minSteps", Value( mMinSteps, ConstraintValueType::NotNegative() ), 0 ) );
	parameters.push_back( Parameter( "tolerance", Value( mTolerance, ConstraintValueType::NotNegative() ), 1.0 ) );

	// Cutoff Parameters
	parameters.push_back( Parameter( "CutoffNonBonded", Value( mNonbondedCutoff, ConstraintValueType::NotNegative() ), 0.0 ) );
	parameters.push_back( Parameter( "CutoffGB", Value( mGBCutoff, ConstraintValueType::NotNegative() ), 0.0 ) );

	parameters.push_back( Parameter( "DeviceID", Value( mPropagatorDevice, ConstraintValueType::NoConstraints() ), -1 ) );
	parameters.push_back( Parameter( "BlockDeviceID", Value( mBlockDevice, ConstraintValueType::NoConstraints() ), -1 ) );
}

STSIntegrator *OpenMMIntegrator::doMake( const vector<Value> &values, ForceGroup *fg ) const {
	return ( STSIntegrator * ) new OpenMMIntegrator( values, fg );
}

// 1 for STS Integrator + 22 for OpenMM Integrator
unsigned int OpenMMIntegrator::getParameterSize() const {
	return 1 + 23;
}

// Figure out OBC scale factors based on the atomic masses.
void OpenMMIntegrator::getObcScaleFactors( std::vector<Real>& scaleFactors ) {
	unsigned int numAtoms = app->positions.size();

	scaleFactors.resize( numAtoms );
	for( unsigned int atomI = 0; atomI < numAtoms; atomI++ ) {
		Real scaleFactor = 0.8;
		Real mass        = app->topology->atoms[atomI].scaledMass;

		if( mass < 1.2 && mass >= 1.0 ) {        // hydrogen
			scaleFactor  = 0.85;
		} else if( mass > 11.8 && mass < 12.2 ) { // carbon
			scaleFactor  = 0.72;
		} else if( mass > 14.0 && mass < 15.0 ) { // nitrogen
			scaleFactor  = 0.79;
		} else if( mass > 15.5 && mass < 16.5 ) { // oxygen
			scaleFactor  = 0.85;
		} else if( mass > 31.5 && mass < 32.5 ) { // sulphur
			scaleFactor  = 0.96;
		} else if( mass > 29.5 && mass < 30.5 ) { // phosphorus
			scaleFactor  = 0.86;
		} else {
			report << plain << " Warning: mass for atom with mass = " << mass << " not recognized." << endr;
		}

		scaleFactors[atomI] = scaleFactor;
	}
}


