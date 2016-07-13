#include <protomol/force/coulomb/CoulombSCPISMForce.h>

#include "protomol/config/ConstraintValueType.h"
#include "protomol/config/Parameter.h"
#include "protomol/config/Value.h"
#include "protomol/type/SimpleTypes.h"

using namespace std;
using namespace ProtoMol;

//____ CoulombForce with DiElectric Term for Implicit Solvation
const string CoulombSCPISMForce::keyword("CoulombSCPISM");

// Default constructor
CoulombSCPISMForce::CoulombSCPISMForce() : D(80) {}

// Constructor with parameters
CoulombSCPISMForce::CoulombSCPISMForce(Real Dval) : D(Dval) {}

void CoulombSCPISMForce::getParameters(vector<Parameter> &parameters) const {
  parameters.push_back
    (Parameter("-D", Value(D, ConstraintValueType::NotNegative()),
               80, Text("Bulk solvent dielectric")));
}

CoulombSCPISMForce CoulombSCPISMForce::make(const vector<Value> &values) {
  return CoulombSCPISMForce(values[0]);
}
