#include <math.h>
#include <protomol/ProtoMolApp.h>
#include <protomol/base/PMConstants.h>
#include <protomol/base/Report.h>
#include <protomol/integrator/base/LangevinVVVRIntegrator.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/topology/TopologyUtilities.h>
#include <protomol/type/Vector3DBlock.h>

#include "protomol/base/MathUtilities.h"
#include "protomol/config/ConstraintValueType.h"
#include "protomol/config/Parameter.h"
#include "protomol/config/Value.h"
#include "protomol/integrator/STSIntegrator.h"
#include "protomol/topology/Atom.h"
#include "protomol/type/SimpleTypes.h"
#include "protomol/type/Vector3D.h"

using namespace std; 
using namespace ProtoMol::Report;
using namespace ProtoMol;
//____ LangevinVVVRIntegrator

const string LangevinVVVRIntegrator::keyword("LangevinVVVR");

LangevinVVVRIntegrator::LangevinVVVRIntegrator() :
  STSIntegrator(), myLangevinTemperature(-1.0), myGamma(-1.0),
  // gamma is in Kcal/ps, myGamma is in Kcal/(fs*INV_TIMEFACTOR)
  mySeed(-1) {}

LangevinVVVRIntegrator::
LangevinVVVRIntegrator(Real timestep, Real LangevinTemperature, Real gamma,
		       int seed, bool correction, ForceGroup *overloadedForces) :
  STSIntegrator(timestep, overloadedForces),
  myLangevinTemperature(LangevinTemperature),
  myGamma(gamma / (1000 * Constant::INV_TIMEFACTOR)),
  // gamma is in Kcal/ps, myGamma is in Kcal/(fs*INV_TIMEFACTOR)
  mySeed(seed),
  timescaleCorrection(correction) {}

void LangevinVVVRIntegrator::initialize(ProtoMolApp *app) {
  STSIntegrator::initialize(app);
  initializeForces();
}

void LangevinVVVRIntegrator::doDrift() {
  Real h;
  if (timescaleCorrection) {
    h = cdt2 * 2.0;
  }
  else {
    h = getTimestep() * Constant::INV_TIMEFACTOR;
  }
  app->positions.intoWeightedAdd(h, app->velocities);
  buildMolecularCenterOfMass(&app->positions, app->topology);
  buildMolecularMomentum(&app->velocities, app->topology);
}

void LangevinVVVRIntegrator::doFirstHalfKick() {
    const unsigned int count = app->positions.size();
    const Real dt = getTimestep() * Constant::INV_TIMEFACTOR; // in fs
    const Real dt2 = dt * 0.5;
    const Real gdt = myGamma*dt;
    const Real gdt2 = gdt * 0.5;
    const Real igdt2 = 1.0 / gdt2;
    const Real bdt = sqrt(exp(-gdt));
    const Real rbdt = sqrt(1-exp(-gdt));
    Real c;
    if (timescaleCorrection==true) {
      c = sqrt(igdt2 * tanh(gdt2));
      cdt2 = c * dt2;
    }
    else {
        cdt2 = dt2;
    }
    const Real variance = Constant::BOLTZMANN * myLangevinTemperature;
    report << debug(990) << "timescaleCorrection = " << timescaleCorrection <<endr;
    report << debug(990) << "dt = " << dt <<endr;
    report << debug(990) << "gdt = " << gdt <<endr;
    report << debug(990) << "bdt = " << bdt <<endr;
    report << debug(990) << "rbdt = " << rbdt <<endr;
    report << debug(990) << "c = " << c <<endr;
    report << debug(990) << "CDT2 = " << cdt2 <<endr;
    report << debug(990) << "bdt = " << bdt << "rbdt = " << rbdt << endr;

    for (unsigned int i = 0; i < count; i++ ) {
        //  Generate gaussian random numbers for each spatial direction
        //force order of generation
        Real rand1 = randomGaussianNumber(mySeed);
        Real rand2 = randomGaussianNumber(mySeed);
        Real rand3 = randomGaussianNumber(mySeed);
        
        //into vector
        Vector3D gaussRandCoord1(rand3, rand2, rand1);
        
        Real mass = app->topology->atoms[i].scaledMass;
        Real sigma = sqrt(variance / mass);
        // semi-update velocities
        app->velocities[i] = app->velocities[i]*bdt
	                     +(*myForces)[i] * cdt2 / mass
	                     +gaussRandCoord1*sigma*rbdt;
    }
    buildMolecularMomentum(&app->velocities, app->topology);
}

void LangevinVVVRIntegrator::doSecondHalfKick() {
    const unsigned int count = app->positions.size();
    const Real dt = getTimestep() * Constant::INV_TIMEFACTOR; // in fs
    const Real dt2 = dt * 0.5;
    const Real gdt = myGamma*dt;
    const Real gdt2 = gdt * 0.5;
    const Real igdt2 = 1.0 / gdt2;
    const Real bdt = sqrt(exp(-gdt));
    const Real rbdt = sqrt(1-exp(-gdt));
    if (timescaleCorrection) {
      cdt2 = sqrt(igdt2 * tanh(gdt2)) * dt2;
    }
    else {
      cdt2 = dt2;
    }
    const Real variance = Constant::BOLTZMANN * myLangevinTemperature;


    for (unsigned int i = 0; i < count; i++ ) {
        //  Generate gaussian random numbers for each spatial direction
        //force order of generation
        Real rand1 = randomGaussianNumber(mySeed);
        Real rand2 = randomGaussianNumber(mySeed);
        Real rand3 = randomGaussianNumber(mySeed);
        
        //into vector
        Vector3D gaussRandCoord1(rand3, rand2, rand1);
        
        Real mass = app->topology->atoms[i].scaledMass;
        Real sigma = sqrt(variance / mass);
        // semi-update velocities
        app->velocities[i] = (app->velocities[i]
			      +(*myForces)[i] * cdt2 / mass) * bdt
	                     +gaussRandCoord1*sigma*rbdt;
    }
    buildMolecularMomentum(&app->velocities, app->topology);
}

long LangevinVVVRIntegrator::run(const long numTimesteps) {
  for (int i = 0; i < numTimesteps; i++) {
    preStepModify();
    doFirstHalfKick();
    doDrift();
    calculateForces();
    doSecondHalfKick();
    postStepModify();
  }
  return numTimesteps;
}


void LangevinVVVRIntegrator::getParameters(vector<Parameter> &parameters)
const {
  STSIntegrator::getParameters(parameters);
  parameters.push_back
    (Parameter("temperature", Value(myLangevinTemperature,
                                    ConstraintValueType::NotNegative())));
  parameters.push_back
    (Parameter("gamma", Value(myGamma * (1000 * Constant::INV_TIMEFACTOR),
                              ConstraintValueType::NotNegative())));
  parameters.push_back
    (Parameter("seed", Value(mySeed, ConstraintValueType::NotNegative()),
               1234));
  parameters.push_back
    (Parameter("correction", Value(timescaleCorrection,ConstraintValueType::NoConstraints()),false,Text("Whether to apply timescale correction")));
}



STSIntegrator *LangevinVVVRIntegrator::doMake(const vector<Value> &values,
                                                 ForceGroup *fg) const {
  return new LangevinVVVRIntegrator(values[0], values[1], values[2],
				    values[3], values[4], fg);
}
