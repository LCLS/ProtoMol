#include <protomol/force/ProteinRestraintForce.h>

#include "protomol/config/ConstraintValueType.h"
#include "protomol/config/Parameter.h"
#include "protomol/config/Value.h"
#include "protomol/type/SimpleTypes.h"

using namespace std;

//____ ProteinRestraintForce
using namespace ProtoMol;
using namespace ProtoMol::Report;

ProteinRestraintForce::ProteinRestraintForce() : sphereK(0), firstEvaluate(true) {}
ProteinRestraintForce::ProteinRestraintForce(Real k, int atom) : sphereK(k), firstEvaluate(true), myAtom(atom - 1){}

const string ProteinRestraintForce::keyword("ProteinRestraint");

unsigned int ProteinRestraintForce::getParameterSize() {
  return 2;
}

void ProteinRestraintForce::getParameters(vector<Parameter> &parameters) const {
  parameters.push_back(
                       Parameter("-tetherconstant",Value(sphereK, ConstraintValueType::NoConstraints()),3.0, Text("Tether constant parameter")));
  parameters.push_back(
		       Parameter("-atom", Value(myAtom, ConstraintValueType::Positive()), 1, Text("Atom to restrain")));
}

ProteinRestraintForce ProteinRestraintForce::make(const vector<Value> &values){
  return (ProteinRestraintForce(values[0], values[1]));
}
