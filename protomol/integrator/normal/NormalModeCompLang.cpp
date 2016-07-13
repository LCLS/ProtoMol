#include <math.h>
#include <protomol/ProtoMolApp.h>
#include <protomol/base/PMConstants.h>
#include <protomol/base/Report.h>
#include <protomol/integrator/normal/ModifierForceProjection.h>
#include <protomol/integrator/normal/NormalModeCompLang.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/topology/TopologyUtilities.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/type/Vector3DBlock.h>
#include <stddef.h>
#include <iostream>

#include "protomol/config/ConstraintValueType.h"
#include "protomol/config/Parameter.h"
#include "protomol/config/Value.h"
#include "protomol/integrator/STSIntegrator.h"
#include "protomol/integrator/normal/NormalModeUtilities.h"
#include "protomol/topology/Atom.h"
#include "protomol/type/SimpleTypes.h"
#include "protomol/type/Vector3D.h"


using namespace std;

using namespace ProtoMol::Report;

using std::string;
using std::vector;

namespace ProtoMol {
  //__________________________________________________ NormalModeCompLang

  const string NormalModeCompLang::keyword( "NormalModeCompLang" );

  NormalModeCompLang::NormalModeCompLang() : STSIntegrator(), NormalModeUtilities()
  {
  }

  NormalModeCompLang::NormalModeCompLang(Real timestep, int firstmode, int nummode, Real gamma, int seed, Real temperature, bool gencn,
                     ForceGroup *overloadedForces) 
    : STSIntegrator(timestep, overloadedForces), 
        NormalModeUtilities( firstmode, nummode, gamma, seed, temperature), genCompNoise(gencn)    
  {
  }


  NormalModeCompLang::~NormalModeCompLang() 
  {  
  }

  void NormalModeCompLang::initialize(ProtoMolApp *app){
    STSIntegrator::initialize(app);
    initializeForces();
    //check valid eigenvectors
    //NM initialization if OK
    int nm_flags = NO_NM_FLAGS;
    if(genCompNoise) nm_flags |= GEN_COMP_NOISE;
    NormalModeUtilities::initialize((int)app->positions.size(), app, myForces, nm_flags); //last int for no complimentary forces or gen noise: GEN_COMP_NOISE
    //
    //
    //take initial C velocites from system and remove non-subspace part
    if(*Q != NULL) subspaceVelocity(&app->velocities, &app->velocities);
    //

     tempV3DBlk.resize(_N);
    temp2V3DBlk.resize(_N);

  }

  //*************************************************************************************
  //****Normal run routine***************************************************************
  //*************************************************************************************

  long NormalModeCompLang::run(const long numTimesteps) {
    //Real h = getTimestep() * Constant::INV_TIMEFACTOR;
    Real actTime;

    if( numTimesteps < 1 )
      return 0;

    //
    //time calculated in forces! so fix here
    if(*Q == NULL)
        report << error << "No Eigenvectors for NormalMode integrator."<<endr;

    actTime = app->topology->time + numTimesteps * getTimestep();

    tempV3DBlk.zero(_N);
    aveForceCount = 0;

    Real ave_energy = 0;

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
      calculateForces();
      report << debug(1) <<"[NormalModeCompLang :: PE "<<app->energies.potentialEnergy()<<endr;
      ave_energy += app->energies.potentialEnergy(); 
      //Real phi = computePhiDihedral (app->topology, &(app->positions), 10);
      //Real psi = computePhiDihedral (app->topology, &(app->positions), 17);
      //cout<<"[NormalModeCompLang ::] phi "<<phi<<" psi "<<psi<<endl;
      //
      genProjGauss(&gaussRandCoord1, app->topology);
      doHalfKick();
      //
      postStepModify();
      //cout <<"Complang step"<<endl;
      
    }

     //fix average, and output
     if(aveForceCount){
       //cout<<"CompLang : aveForceCount "<<aveForceCount<<endl;
       for( unsigned int i=0;i < app->positions.size(); i++) tempV3DBlk[i] /= (Real)aveForceCount;
       //for( unsigned int i=0;i < app->positions.size(); i++) cout<<"IN Brownian : mean force for atom "<<i<<" "<<tempV3DBlk[i]<<endl;
       cout <<"Average Potential Energy "<<ave_energy/(Real)aveForceCount<<endl;
     }

    //fix time
    app->topology->time = actTime;
    //
    return numTimesteps;
  }  

  void NormalModeCompLang::doDrift() {
      const Real h = getTimestep() * Constant::INV_TIMEFACTOR;
      app->positions.intoWeightedAdd(h, app->velocities);
      buildMolecularCenterOfMass(&app->positions, app->topology);
  }

  void NormalModeCompLang::doHalfKick() {
    const Real dt = getTimestep() * Constant::INV_TIMEFACTOR; // in fs
    const Real fdt = ( 1.0 - exp( -0.5 * myGamma * dt ) ) / myGamma;
    const Real vdt = exp(-0.5*myGamma*dt);
    const Real ndt = sqrt( ( 1.0 - exp( -myGamma * dt ) ) / (2.0 * myGamma) ); //was sqrt( fdt );
    const Real sqrtFCoverM = sqrt( 2.0 * Constant::BOLTZMANN * myTemp * myGamma );

    for( int i = 0; i < _N; i++ ) {
        // semi-update velocities
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

  void NormalModeCompLang::getParameters(vector<Parameter>& parameters) const {
    STSIntegrator::getParameters(parameters);
    parameters.push_back(Parameter("firstmode",Value(firstMode,ConstraintValueType::NotNegative()),1,Text("First mode to use in set")));
    parameters.push_back(Parameter("numbermodes",Value(numMode,ConstraintValueType::NotNegative()),1,Text("Number of modes propagated")));
    parameters.push_back(Parameter("gamma",Value(myGamma*(1000 * Constant::INV_TIMEFACTOR),ConstraintValueType::NotNegative()),80.0,Text("Langevin Gamma")));
    parameters.push_back(Parameter("seed",Value(mySeed,ConstraintValueType::NotNegative()),1234,Text("Langevin random seed")));
    parameters.push_back(Parameter("temperature",Value(myTemp,ConstraintValueType::NotNegative()),300.0,Text("Langevin temperature")));
    parameters.push_back(Parameter("gencompnoise",Value(genCompNoise,ConstraintValueType::NoConstraints()),false,Text("Generate complimentary noise")));
 }

  STSIntegrator* NormalModeCompLang::doMake(const vector<Value>& values,ForceGroup* fg)const{
    return new NormalModeCompLang(values[0],values[1],values[2],values[3],values[4],values[5],values[6],fg);
  }

  void NormalModeCompLang::addModifierAfterInitialize(){
    adoptPostForceModifier(new ModifierForceProjection(this));
    STSIntegrator::addModifierAfterInitialize();
  }

        //override force projection for post force modifier, and find average in complement space
        void NormalModeCompLang::forceProjection(){
                unsigned int count = myForces->size();
                if((*Q) != NULL){
                        // myForces has total forces
                        //for( unsigned int i=0;i < count; i++) temp2V3DBlk[i] = (*myForces)[i];
                        temp2V3DBlk.resize(count);
                        temp2V3DBlk.intoAssign(*myForces);
                        // project myForces onto fast subspace
                        subspaceForce(myForces, myForces);
                        // difference between old myForces stored in temp2V3DBlk and myForces
                        // gives us instantaneous slow force
                        for( unsigned int i=0;i < count; i++) temp2V3DBlk[i] -= (*myForces)[i];
                        //for( unsigned int j = 0; j < app->positions.size(); ++j) cout<<"In Brownian in force projection [22]: atom "<<j<<" "<<temp2V3DBlk[j]<<endl;
                        // add slow force to running sum in tempV3DBlk[i]
                        for( unsigned int i=0;i < count; i++) tempV3DBlk[i] += temp2V3DBlk[i];
                        //tempV3DBlk.intoAdd(temp2V3DBlk);
                        aveForceCount++;
                        //for( unsigned int j = 0; j < app->positions.size(); ++j) cout<<"In Brownian in force projection (accumulated mean force) : atom "<<j<<" "<<tempV3DBlk[j]<<endl;

                }

        }


}
