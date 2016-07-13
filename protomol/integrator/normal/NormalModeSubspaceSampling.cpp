#include <protomol/ProtoMolApp.h>
#include <protomol/base/PMConstants.h>
#include <protomol/base/Report.h>
#include <protomol/integrator/normal/ModifierForceProjection.h>
#include <protomol/integrator/normal/NormalModeSubspaceSampling.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/type/Vector3DBlock.h>
#include <stddef.h>
#include <iostream>

#include "protomol/config/ConstraintValueType.h"
#include "protomol/config/Parameter.h"
#include "protomol/config/Value.h"
#include "protomol/integrator/Integrator.h"
#include "protomol/integrator/MTSIntegrator.h"
#include "protomol/integrator/StandardIntegrator.h"
#include "protomol/integrator/normal/NormalModeUtilities.h"
#include "protomol/topology/Atom.h"
#include "protomol/type/SimpleTypes.h"
#include "protomol/type/Vector3D.h"

using namespace std;


using namespace ProtoMol::Report;

using std::string;
using std::vector;

namespace ProtoMol {
  //__________________________________________________ NormalModeSubspaceSampling

  const string NormalModeSubspaceSampling::keyword( "NormalModeSubspaceSampling" );

  NormalModeSubspaceSampling::NormalModeSubspaceSampling() : MTSIntegrator(), NormalModeUtilities()
  {

      myBottomNormalMode = NULL;
      ex0=NULL;//####diagnostics
  }

  NormalModeSubspaceSampling::NormalModeSubspaceSampling(int cycles, int firstmode, int nummode, Real temperature,
                    bool instf, Real sl, int do_up,
                        ForceGroup *overloadedForces, StandardIntegrator *nextIntegrator) 
    : MTSIntegrator(cycles, overloadedForces, nextIntegrator), 
        NormalModeUtilities( firstmode, nummode, 80.0, 1234 , temperature),     
            instForce(instf), dh(sl), doUpdate(do_up)

  {
      myBottomNormalMode = NULL;
      ex0=NULL;//####diagnostics
  }


  NormalModeSubspaceSampling::~NormalModeSubspaceSampling() 
  {        
      if(ex0!=NULL) delete ex0;//####diagnostics
  }

  void NormalModeSubspaceSampling::initialize(ProtoMolApp* app){
    MTSIntegrator::initialize(app);
    //
    //point to bottom integrator, for average force
    myBottomNormalMode  = dynamic_cast<NormalModeUtilities*>(bottom());
    //check valid eigenvectors
    //NM initialization if OK
    NormalModeUtilities::initialize((int)app->positions.size(), app, myForces, NO_NM_FLAGS); //last for non-complimentary forces
    //
    //do first force calculation, and remove non sub-space part
    app->energies.clear();	//Need this or initial error, due to inner integrator energy?
    initializeForces();
    //
    //take initial C velocites from system and remove non-subspace part
    if(*Q != NULL) subspaceVelocity(&app->velocities, &app->velocities);

    //std::cout<<"NormalModeSubspaceSampling::initialize x_0 size "<<x_0.size()<<", "<<tmpC<<" , app size "<<app->positions.size()<<std::endl;
    //if (x_0.size()>0) {
    	//modeSpaceProj(tmpC, &app->positions, &x_0); 
    //}

  }

  //*************************************************************************************
  //****Normal run routine***************************************************************
  //*************************************************************************************

  long NormalModeSubspaceSampling::run(const long numTimesteps) {
    //Real h = getTimestep() * Constant::INV_TIMEFACTOR;
    Real h;
    Real actTime;

    if (dh >0) h = dh;
    else  h = getTimestep() * Constant::INV_TIMEFACTOR;

    if( numTimesteps < 1 )
      return 0;

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
      //
      //constraints?
      app->energies.clear();
      save_pos = *&app->positions;
      //run minimizer if any remaining modes
      if(testRemainingModes()) myNextIntegrator->run(myCycleLength); //cyclelength 
      else cout<<"Not doing brownian "<<endl;
      *&app->positions =  save_pos;

      if(*Q == NULL){	//rediagonalize?
            app->topology->time = actTime - (i - numTimesteps) * getTimestep();
            if(myPreviousIntegrator == NULL) 
                report << error << "[NormalModeSubspaceSampling::Run] Re-diagonalization forced with NormalModeSubspaceSampling as outermost Integrator. Aborting."<<endr;
            return 0;
      }
      //#################Put averaged force code here##############################
      //calculate sub space forces, just do this for the energy
      app->energies.clear();
      calculateForces();
      //transfer the real average force, Note: force fields must be identical
      if(!instForce){
        //myForces->intoAssign(myBottomNormalMode->tempV3DBlk);
        cout<<"SubspaceSampling : Grabbing mean force from CompLang..."<<endl;
        for( unsigned int i=0;i < app->positions.size(); i++) {
            (*myForces)[i] = myBottomNormalMode->tempV3DBlk[i];
            //cout<<"In Mori, mean force, Atom "<<i<<" "<<(*myForces)[i]<<", from NM Brownian "<<myBottomNormalMode->tempV3DBlk[i]<<endl;
        }
        
      }
      //cout<<"NormalModeSubspaceSampling:: h "<<h<<endl;
      //update position
      if (doUpdate) {
       cout <<"Updating"<<endl;
      for (unsigned int i=0;i<app->positions.size();i++) {
          app->positions[i] += ((*myForces)[i]*h*h)/app->topology->atoms[i].scaledMass;
          //app->positions[i] += ((*myForces)[i]*dh*dh)/app->topology->atoms[i].scaledMass;
      }
      }
      //###########################################################################
      //
      //
      postStepModify();
    }	
    //fix time
    app->topology->time = actTime;

    //project into the subspace
    //reparameterization code will be passed this array
    //modeSpaceProj(tmpC, &app->positions, &x_0);

    //std::cout<<"Projection of point in the subspace..."<<std::endl;
    //for(int i=0;i<_rfM;i++) std::cout<<tmpC[i]<<" ";
    //std::cout<<std::endl;
    //
    return numTimesteps;
  }


  //*************************************************************************************
  //****Output int paramiters************************************************************
  //*************************************************************************************

  void NormalModeSubspaceSampling::getParameters(vector<Parameter>& parameters) const {
    MTSIntegrator::getParameters(parameters);
    parameters.push_back(Parameter("firstmode",Value(firstMode,ConstraintValueType::NotNegative()),1,Text("First mode to use in set")));
    parameters.push_back(Parameter("numbermodes",Value(numMode,ConstraintValueType::NotNegative()),1,Text("Number of modes propagated")));
    parameters.push_back(Parameter("temperature",Value(myTemp,ConstraintValueType::NotNegative()),300.0,Text("Langevin temperature")));
    parameters.push_back(Parameter("instForce",Value(instForce,ConstraintValueType::NoConstraints()),false,Text("Use instantaneous force")));
    parameters.push_back(Parameter("stepLength",Value(dh,ConstraintValueType::NotNegative()),0,Text("Steplength for integration on collective variable space")));
    parameters.push_back(Parameter("doUpdate",Value(doUpdate,ConstraintValueType::NotNegative()),1,Text("Steplength for integration on collective variable space")));
    
    //####
  }

  MTSIntegrator* NormalModeSubspaceSampling::doMake(const vector<Value>& values,ForceGroup* fg, StandardIntegrator *nextIntegrator)const{
    return new NormalModeSubspaceSampling(values[0],values[1],values[2],values[3],values[4],values[5],values[6],fg,nextIntegrator);
  }

  void NormalModeSubspaceSampling::addModifierAfterInitialize(){
    adoptPostForceModifier(new ModifierForceProjection(this));
    MTSIntegrator::addModifierAfterInitialize();
  }

}
