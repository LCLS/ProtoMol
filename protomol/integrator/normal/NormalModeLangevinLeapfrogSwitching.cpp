#include <math.h>
#include <protomol/ProtoMolApp.h>
#include <protomol/base/PMConstants.h>
#include <protomol/base/Report.h>
#include <protomol/integrator/normal/ModifierForceProjection.h>
#include <protomol/integrator/normal/NormalModeLangevinLeapfrogSwitching.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/topology/TopologyUtilities.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/type/Vector3DBlock.h>
#include <stddef.h>

#include "protomol/base/MathUtilities.h"
#include "protomol/config/ConstraintValueType.h"
#include "protomol/config/Parameter.h"
#include "protomol/config/Value.h"
#include "protomol/integrator/MTSIntegrator.h"
#include "protomol/integrator/StandardIntegrator.h"
#include "protomol/integrator/normal/NormalModeUtilities.h"
#include "protomol/topology/Atom.h"
#include "protomol/type/EigenvectorInfo.h"
#include "protomol/type/SimpleTypes.h"

using namespace std;

using namespace ProtoMol::Report;

using std::string;
using std::vector;

namespace ProtoMol {
  //__________________________________________________ NormalModeLangevinLeapfrogSwitching

  const string NormalModeLangevinLeapfrogSwitching::keyword( "NormalModeLangevinLeapfrogSwitching" );

  NormalModeLangevinLeapfrogSwitching::NormalModeLangevinLeapfrogSwitching() : MTSIntegrator(), NormalModeUtilities()
  {
  }

  NormalModeLangevinLeapfrogSwitching::NormalModeLangevinLeapfrogSwitching(int cycles, int firstmode, int nummode, Real gammaInside, Real gammaOutside, int seed, Real temperature, bool gencn,
									   Real switchon, Real cutoff,
                     ForceGroup *overloadedForces, StandardIntegrator *nextIntegrator) 
    : MTSIntegrator(cycles, overloadedForces, nextIntegrator), 
      NormalModeUtilities( firstmode, nummode, gammaInside, seed, temperature), genCompNoise(gencn),
      myGammaInside(gammaInside / (1000 * Constant::INV_TIMEFACTOR)),
      myGammaOutside(gammaOutside / (1000 * Constant::INV_TIMEFACTOR)),
      mySwitchOn(switchon),
      mySwitchOn2(switchon * switchon),
      myCutoff(cutoff),
      myCutoff2(cutoff * cutoff),
      mySwitch1(1.0 / power<3>(cutoff * cutoff - switchon * switchon)),
      mySwitch2(cutoff * cutoff - 3.0 * switchon * switchon),
      mySwitch3(4.0 / power<3>(cutoff * cutoff - switchon * switchon))
  {
  }


  NormalModeLangevinLeapfrogSwitching::~NormalModeLangevinLeapfrogSwitching() 
  {  
  }

  void NormalModeLangevinLeapfrogSwitching::initialize(ProtoMolApp *app){
    MTSIntegrator::initialize(app);
    //check valid eigenvectors
    //NM initialization if OK
    int nm_flags = NO_NM_FLAGS;
    if(genCompNoise) nm_flags |= GEN_COMP_NOISE;
    NormalModeUtilities::initialize((int)app->positions.size(), app, myForces, nm_flags); //last int for no complimentary forces or gen noise: GEN_COMP_NOISE
    //
    //do first force calculation, and remove non sub-space part
    app->energies.clear();	//Need this or initial error, due to inner integrator energy?
    initializeForces();
    //
    //take initial C velocites from system and remove non-subspace part
    if(*Q != NULL) subspaceVelocity(&app->velocities, &app->velocities);
    //

  }

  //*************************************************************************************
  //****Normal run routine***************************************************************
  //*************************************************************************************

  long NormalModeLangevinLeapfrogSwitching::run(const long numTimesteps) {
    Real actTime;

    if( numTimesteps < 1 ) return 0;

    //check valid eigenvectors
    if(*Q == NULL)
        report << error << "No Eigenvectors for NormalMode integrator."<<endr;
    //
    //time calculated in forces! so fix here
    actTime = app->topology->time + numTimesteps * getTimestep();
    //
    //main loop
    for( int i = 0; i < numTimesteps; i++ ) {
      //****main loop*************************************
      preStepModify();
      genProjGauss(&gaussRandCoord1, app->topology);
      doHalfKick();
      //
      //nmlDrift(&app->positions, &app->velocities, h, app->topology);
      doDrift();
      //constraints?
      app->energies.clear();
      //run minimizer if any remaining modes
      if(testRemainingModes()) myNextIntegrator->run(myCycleLength); //cyclelength 
      if(app->eigenInfo.reDiagonalize){	//rediagonalize?
            app->topology->time = actTime + (i - numTimesteps) * getTimestep();
            if(myPreviousIntegrator == NULL) 
                report << error << "[NormalModeLangevinLeapfrogSwitching::Run] Re-diagonalization forced with NormalModeLangevinLeapfrogSwitching as outermost Integrator. Aborting."<<endr;
            return i;
      }
      //calculate sub space forces
      app->energies.clear();
      calculateForces();
      //
      genProjGauss(&gaussRandCoord1, app->topology);
      doHalfKick();
      //
      postStepModify();
    }
    //fix time
    app->topology->time = actTime;
    
    return numTimesteps;
  }  

  void NormalModeLangevinLeapfrogSwitching::doDrift() {
      const Real h = getTimestep() * Constant::INV_TIMEFACTOR;
      app->positions.intoWeightedAdd(h, app->velocities);
      buildMolecularCenterOfMass(&app->positions, app->topology);
  }

  void NormalModeLangevinLeapfrogSwitching::doHalfKick() {
    for(unsigned int i = 0; i < _N; i++) 
      {
	const Vector3D diff = (app->positions)[i] - myCenterOfMass;
	const Real distance = diff.norm();
	const Real distSquared = distance * distance;
	
	// switch
	Real switchValue;
	if(distSquared > myCutoff2)
	  {
	    switchValue = 0.0;
	  }
	else if(distSquared >= mySwitchOn2)
	  {
	    const Real c2 = myCutoff2 - distSquared;
	    const Real c4 = c2 * (mySwitch2 + 2.0 * distSquared);
	    switchValue = mySwitch1 * c2 * c4;
	  }
	else
	  {
	    switchValue = 1.0;
	  }

	const Real complimentSwitchValue = 1.0 - switchValue;

	const Real myGamma = switchValue * myGammaInside + (1.0 - switchValue) * myGammaOutside;

	//	cout << "DistGamma " << distance << " " << myGamma << endl;

	const Real dt = getTimestep() * Constant::INV_TIMEFACTOR; // in fs
	const Real fdt = ( 1.0 - exp( -0.5 * myGamma * dt ) ) / myGamma;
	const Real vdt = exp(-0.5*myGamma*dt);
	const Real ndt = sqrt( ( 1.0 - exp( -myGamma * dt ) ) / (2.0 * myGamma) ); //was sqrt( fdt );
	const Real sqrtFCoverM = sqrt( 2.0 * Constant::BOLTZMANN * myTemp * myGamma );

        app->velocities[i] = app->velocities[i]*vdt
                                +(*myForces)[i] * fdt / app->topology->atoms[i].scaledMass
                                    +gaussRandCoord1[i]*sqrtFCoverM*ndt;
	

	
      }

    subspaceVelocity(&app->velocities, &app->velocities);
    buildMolecularMomentum(&app->velocities, app->topology);
  }

  //*************************************************************************************
  //****Output int paramiters************************************************************
  //*************************************************************************************

  void NormalModeLangevinLeapfrogSwitching::getParameters(vector<Parameter>& parameters) const {
    MTSIntegrator::getParameters(parameters);
    parameters.push_back(Parameter("firstmode",Value(firstMode,ConstraintValueType::NotNegative()),1,Text("First mode to use in set")));
    parameters.push_back(Parameter("numbermodes",Value(numMode,ConstraintValueType::NotNegative()),1,Text("Number of modes propagated")));
    parameters.push_back(Parameter("gammaInside",Value(myGammaInside*(1000 * Constant::INV_TIMEFACTOR),ConstraintValueType::NotNegative()),80.0,Text("Langevin Gamma")));
    parameters.push_back(Parameter("gammaOutside",Value(myGammaInside*(1000 * Constant::INV_TIMEFACTOR),ConstraintValueType::NotNegative()),80.0,Text("Langevin Gamma")));
    parameters.push_back(Parameter("seed",Value(mySeed,ConstraintValueType::NotNegative()),1234,Text("Langevin random seed")));
    parameters.push_back(Parameter("temperature",Value(myTemp,ConstraintValueType::NotNegative()),300.0,Text("Langevin temperature")));
    parameters.push_back(Parameter("gencompnoise",Value(genCompNoise,ConstraintValueType::NoConstraints()),false,Text("Generate complimentary noise")));
    parameters.push_back(Parameter("switchon", Value(mySwitchOn, ConstraintValueType::NotNegative()), 10, Text("Switch on")));
    parameters.push_back(Parameter("cutoff", Value(myCutoff, ConstraintValueType::NotNegative()), 10, Text("Cutoff")));
 }

  MTSIntegrator* NormalModeLangevinLeapfrogSwitching::doMake(const vector<Value>& values,ForceGroup* fg, StandardIntegrator *nextIntegrator)const{
    return new NormalModeLangevinLeapfrogSwitching(values[0],values[1],values[2],values[3],values[4],values[5],values[6],values[7],values[8],values[9],fg,nextIntegrator);
  }

  void NormalModeLangevinLeapfrogSwitching::addModifierAfterInitialize(){
    adoptPostForceModifier(new ModifierForceProjection(this));
    MTSIntegrator::addModifierAfterInitialize();
  }

}
