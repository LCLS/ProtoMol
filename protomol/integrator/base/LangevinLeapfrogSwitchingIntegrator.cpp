#include <math.h>
#include <protomol/ProtoMolApp.h>
#include <protomol/base/PMConstants.h>
#include <protomol/integrator/base/LangevinLeapfrogSwitchingIntegrator.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/topology/TopologyUtilities.h>
#include <protomol/type/Vector3DBlock.h>

#include "protomol/base/MathUtilities.h"
#include "protomol/config/ConstraintValueType.h"
#include "protomol/config/Parameter.h"
#include "protomol/config/Value.h"
#include "protomol/integrator/STSIntegrator.h"
#include "protomol/topology/Atom.h"

using namespace std; 
using namespace ProtoMol::Report;
using namespace ProtoMol;

//____ LangevinLeapfrogSwitchingIntegrator

const string LangevinLeapfrogSwitchingIntegrator::keyword("LangevinLeapfrogSwitching");

LangevinLeapfrogSwitchingIntegrator::LangevinLeapfrogSwitchingIntegrator() :
  STSIntegrator(), myLangevinTemperature(-1.0), myGammaInside(-1.0),
  // gamma is in Kcal/ps, myGamma is in Kcal/(fs*INV_TIMEFACTOR)
  mySeed(-1) {}

LangevinLeapfrogSwitchingIntegrator::
LangevinLeapfrogSwitchingIntegrator(Real timestep, Real LangevinTemperature, Real gammaInside,
				    Real gammaOutside, int seed, Real switchon,
				    Real cutoff, ForceGroup *overloadedForces) :
  STSIntegrator(timestep, overloadedForces),
  myLangevinTemperature(LangevinTemperature),
  // gamma is in Kcal/ps, myGamma is in Kcal/(fs*INV_TIMEFACTOR)
  myGammaInside(gammaInside / (1000 * Constant::INV_TIMEFACTOR)),
  myGammaOutside(gammaOutside / (1000 * Constant::INV_TIMEFACTOR)),
  mySeed(seed),
  mySwitchOn(switchon),
  mySwitchOn2(switchon * switchon),
  myCutoff(cutoff),
  myCutoff2(cutoff * cutoff),
  mySwitch1(1.0 / power<3>(cutoff * cutoff - switchon * switchon)),
  mySwitch2(cutoff * cutoff - 3.0 * switchon * switchon),
  mySwitch3(4.0 / power<3>(cutoff * cutoff - switchon * switchon)){}

void LangevinLeapfrogSwitchingIntegrator::initialize(ProtoMolApp *app) {
  STSIntegrator::initialize(app);
  initializeForces();
  
  myCenterOfMass = centerOfMass(&app->positions, app->topology);

}

void LangevinLeapfrogSwitchingIntegrator::doDrift() {
  const Real h = getTimestep() * Constant::INV_TIMEFACTOR;
  app->positions.intoWeightedAdd(h, app->velocities);
  buildMolecularCenterOfMass(&app->positions, app->topology);
  buildMolecularMomentum(&app->velocities, app->topology);
}

void LangevinLeapfrogSwitchingIntegrator::doHalfKick() {
    const unsigned int count = app->positions.size();
    
    for (unsigned int i = 0; i < count; i++ ) {
      const Vector3D diff = (app->positions)[i] - myCenterOfMass;
      const Real distance = diff.norm();
      const Real distSquared = distance * distance;
      
      // switch
      Real switchValue;
      if(distSquared > myCutoff2) {
	switchValue = 0.0;
      }
      else if(distSquared >= mySwitchOn2) {
	const Real c2 = myCutoff2 - distSquared;
	const Real c4 = c2 * (mySwitch2 + 2.0 * distSquared);
	switchValue = mySwitch1 * c2 * c4;
      } else {
	switchValue = 1.0;
      }

      Real complimentSwitchValue = 1.0 - switchValue;

      const Real myGamma = switchValue * myGammaInside + complimentSwitchValue * myGammaOutside;

      const Real dt = getTimestep() * Constant::INV_TIMEFACTOR; // in fs
      const Real fdt = ( 1.0 - exp( -0.5 * myGamma * dt ) ) / myGamma;
      const Real vdt = exp(-0.5*myGamma*dt);
      const Real ndt = sqrt( ( 1.0 - exp( -myGamma * dt ) ) / (2.0 * myGamma) );
      const Real forceConstant = 2 * Constant::BOLTZMANN * myLangevinTemperature *
	myGamma;

      //  Generate gaussian random numbers for each spatial direction
      //force order of generation
      Real rand1 = randomGaussianNumber(mySeed);
      Real rand2 = randomGaussianNumber(mySeed);
      Real rand3 = randomGaussianNumber(mySeed);
        
      //into vector
      Vector3D gaussRandCoord1(rand3, rand2, rand1);
        
      Real mass = app->topology->atoms[i].scaledMass;
      Real sqrtFCoverM = sqrt(forceConstant / mass);
      // semi-update velocities
      app->velocities[i] = app->velocities[i]*vdt
	+(*myForces)[i] * fdt / mass
	+gaussRandCoord1*sqrtFCoverM*ndt;
    }
    buildMolecularMomentum(&app->velocities, app->topology);
}

void LangevinLeapfrogSwitchingIntegrator::getParameters(vector<Parameter> &parameters)
const {
  STSIntegrator::getParameters(parameters);
  parameters.push_back
    (Parameter("temperature", Value(myLangevinTemperature,
                                    ConstraintValueType::NotNegative())));
  parameters.push_back
    (Parameter("gammaInside", Value(myGammaInside * (1000 * Constant::INV_TIMEFACTOR),
                              ConstraintValueType::NotNegative())));

  parameters.push_back
    (Parameter("gammaOutside", Value(myGammaOutside * (1000 * Constant::INV_TIMEFACTOR),
                              ConstraintValueType::NotNegative())));

  parameters.push_back
    (Parameter("seed", Value(mySeed, ConstraintValueType::NotNegative()),
               1234));

  parameters.push_back
    (Parameter("switchon", Value(mySwitchOn),
	       ConstraintValueType::NotNegative()));

  parameters.push_back
  (Parameter("cutoff", Value(myCutoff),
	     ConstraintValueType::NotNegative()));

}

STSIntegrator *LangevinLeapfrogSwitchingIntegrator::doMake(const vector<Value> &values,
                                                 ForceGroup *fg) const {
  return new LangevinLeapfrogSwitchingIntegrator(values[0], values[1], values[2],
						 values[3], values[4], values[5],
						 values[6], fg);
}
