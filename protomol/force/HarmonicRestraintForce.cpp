#include <protomol/force/HarmonicRestraintForce.h>

#include "protomol/base/MathUtilities.h"
#include "protomol/config/ConstraintValueType.h"
#include "protomol/config/Parameter.h"
#include "protomol/config/Value.h"
#include "protomol/type/SimpleTypes.h"

using namespace std;

//____ HarmonicRestraintForce
using namespace ProtoMol;
using namespace ProtoMol::Report;

HarmonicRestraintForce::HarmonicRestraintForce() : sphereK(0), sphereradius(0), firstEvaluate(true) {}
HarmonicRestraintForce::HarmonicRestraintForce(Real k, Real r, Real switchon, Real cutoff) : sphereK(k), sphereradius(r), firstEvaluate(true),
          mySwitchon(switchon), mySwitchon2(switchon * switchon), myCutoff(cutoff),
          myCutoff2(cutoff * cutoff),
          mySwitch1(1.0 / power<3>(cutoff * cutoff - switchon * switchon)),
          mySwitch2(cutoff * cutoff - 3.0 * switchon * switchon),
          mySwitch3(4.0 / power<3>(cutoff * cutoff - switchon * switchon)){}

const string HarmonicRestraintForce::keyword("HarmonicRestraint");

unsigned int HarmonicRestraintForce::getParameterSize() {
  return 4;
}

void HarmonicRestraintForce::getParameters(vector<Parameter> &parameters) const {
  parameters.push_back(
                       Parameter("-sphereconstant",Value(sphereK, ConstraintValueType::NoConstraints()),3.0, Text("Sphere constant parameter")));
  parameters.push_back(
                       Parameter("-sphereradius",Value(sphereradius, ConstraintValueType::NoConstraints()), 12.0, Text("Sphere radius parameter")));
  parameters.push_back(
                       Parameter("-switchstart", Value(mySwitchon, ConstraintValueType::NotNegative()), 12.0, Text("Switch start = radius")));
  parameters.push_back(
                       Parameter("-switchend", Value(myCutoff, ConstraintValueType::Positive()), 15.0, Text("Switch end=1 > radius")));
}

HarmonicRestraintForce HarmonicRestraintForce::make(const vector<Value> &values){
  return (HarmonicRestraintForce(values[0], values[1], values[2], values[3]));
}
