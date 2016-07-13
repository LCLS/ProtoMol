#include <protomol/force/GB/GBACEForce.h>
#include <protomol/parallel/Parallel.h>
#include <cmath>

#include "protomol/base/MathUtilities.h"
#include "protomol/config/ConstraintValueType.h"
#include "protomol/config/Parameter.h"
#include "protomol/config/Value.h"
#include "protomol/topology/Atom.h"
#include "protomol/topology/ExclusionTable.h"
#include "protomol/topology/GenericTopology.h"
#include "protomol/type/ScalarStructure.h"
#include "protomol/type/SimpleTypes.h"

namespace ProtoMol {
class Vector3DBlock;
}  // namespace ProtoMol

using namespace ProtoMol;
using namespace ProtoMol::Report;

#define PI 3.14169
#define SIXPOW 6

const string GBACEForce::keyword("GBACEForce");

GBACEForce::GBACEForce() : sigma(0), rho_s(0) {}
GBACEForce::GBACEForce(Real s, Real r_s) : sigma(s), rho_s(r_s) {}

void GBACEForce::operator()(Real &energy, Real &force, Real distSquared,
                Real rDistSquared, const Vector3D &,
                const GenericTopology *topo, int atom1, int atom2,
                ExclusionClass excl) const {
  
  
  //r_{ij}
  Real dist = sqrt(distSquared);
  
  Real radius_i = topo->atoms[atom1].myGBSA_T->vanDerWaalRadius; 
  Real radius_j = topo->atoms[atom2].myGBSA_T->vanDerWaalRadius;
  
  //offset radii ({\tilde{\rho}_{j}})
  Real offsetRadius_i = radius_i - topo->atoms[atom1].myGBSA_T->offsetRadius;
  Real offsetRadius_j = radius_j - topo->atoms[atom2].myGBSA_T->offsetRadius;
  
  
  //Scaling factors
  Real S_i = topo->atoms[atom1].myGBSA_T->scalingFactor;
  Real S_j = topo->atoms[atom2].myGBSA_T->scalingFactor;
  
  Real psi_i, psi_j;
  
  Real tanhparam_i, tanhparam_j;
  
  psi_i = topo->atoms[atom1].myGBSA_T->PsiValue;
  psi_j = topo->atoms[atom2].myGBSA_T->PsiValue;
  tanhparam_i = topo->alphaObc*psi_i - topo->betaObc*psi_i*psi_i + topo->gammaObc*psi_i*psi_i*psi_i;
  tanhparam_j = topo->alphaObc*psi_j - topo->betaObc*psi_j*psi_j + topo->gammaObc*psi_j*psi_j*psi_j;
  
  Real bornRad_i = topo->atoms[atom1].myGBSA_T->bornRad;
  Real bornRad_j = topo->atoms[atom2].myGBSA_T->bornRad;
  
  Real ratio_i, ratio_j;
  
  //Check Equation (14) for ACE nonpolar solvation potential
  if (!topo->atoms[atom1].myGBSA_T->ACEPotentialCount) {
    Real p1 = radius_i/topo->atoms[atom1].myGBSA_T->bornRad;
    ratio_i = power(p1,SIXPOW);
    topo->atoms[atom1].myGBSA_T->ACEPotential = 4*PI*sigma*(radius_i + rho_s)*(radius_i + rho_s)*ratio_i;
    topo->atoms[atom1].myGBSA_T->ACEPotentialCount = 1;
  }
  
  if (!topo->atoms[atom2].myGBSA_T->ACEPotentialCount) {
    Real p2 = radius_j/topo->atoms[atom2].myGBSA_T->bornRad;
    ratio_j = power(p2,SIXPOW);
    topo->atoms[atom2].myGBSA_T->ACEPotential = 4*PI*sigma*(radius_j + rho_s)*(radius_j + rho_s)*ratio_j;
    topo->atoms[atom2].myGBSA_T->ACEPotentialCount = 1;
  }
  
  Real c1 = 24*PI*sigma;
  
  Real c_i = (radius_i + rho_s)*(radius_i + rho_s)*pow(radius_i,SIXPOW);
  Real c_j = (radius_j + rho_s)*(radius_j + rho_s)*pow(radius_j,SIXPOW);
  
  //It would be nice if we can precalculate this terms.
  Real Lij, Uij;
  
  if (offsetRadius_i >=  dist + S_j*offsetRadius_j) {
    Lij = 1;
    Uij = 1;
  }else {
    Lij =(offsetRadius_i > fabs(dist - S_j*offsetRadius_j)) ? offsetRadius_i : fabs(dist - S_j*offsetRadius_j);
    Uij = dist + S_j*offsetRadius_j;
  }
  
  Real Lji, Uji;
  
  if (offsetRadius_j >=  dist + S_i*offsetRadius_i) {
    Lji = 1;
    Uji = 1;
  }else {
    Lji = (offsetRadius_j > abs(dist - S_i*offsetRadius_i)) ? offsetRadius_j : abs(dist - S_i*offsetRadius_i);         
    Uji = dist + S_i*offsetRadius_i;
  }
  
  //Derivatives for calculation of the derivative of the born radii
  Real dLijdrij, dUijdrij, dCijdrij;
  
  //Check Equations (11-13)
  if (offsetRadius_i <= (dist - S_j*offsetRadius_j)) dLijdrij = 1;
  else dLijdrij = 0;
  
  if (offsetRadius_i < (dist + S_j*offsetRadius_j)) dUijdrij = 1;
  else dUijdrij = 0;
  
  if (offsetRadius_i <= (S_j*offsetRadius_j - dist)) dCijdrij = 2*(1/Lij)*(1/Lij)*dLijdrij;
  else dCijdrij = 0;
  
  //Derivatives for calculation of the derivative of the born radii
  Real dLjidrij, dUjidrij, dCjidrij;
  
  if (offsetRadius_j <= (dist - S_i*offsetRadius_i)) dLjidrij = 1;
  else dLjidrij = 0;
  
  if (offsetRadius_j < (dist + S_i*offsetRadius_i)) dUjidrij = 1;
  else dUjidrij = 0;
  
  if (offsetRadius_j <= (S_i*offsetRadius_i - dist)) dCjidrij = 2*(1/Lji)*(1/Lji)*dLjidrij;
  else dCjidrij = 0;
  
  //tanhk angle
  Real tanh_i = tanh(tanhparam_i);
  
  //tanhi angle
  Real tanh_j = tanh(tanhparam_j);
  
  //tanhk derivative
  Real tanhparam_derv_i = (topo->alphaObc - 2*topo->betaObc*psi_i + 3*topo->gammaObc*psi_i*psi_i);
  
  //tanhi derivative
  Real tanhparam_derv_j = (topo->alphaObc - 2*topo->betaObc*psi_j + 3*topo->gammaObc*psi_j*psi_j);
  
  Real S_i_term = (S_i*offsetRadius_i)/dist;
  Real S_j_term = (S_j*offsetRadius_j)/dist;
  
  //Check Equation (10) for the derivative of the burial term
  Real dBTidrij = -0.5*dLijdrij*(1/(Lij*Lij)) + 0.5*dUijdrij*(1/(Uij*Uij)) + 0.125*((1/(Uij*Uij)) - (1/(Lij*Lij))) + 0.125*dist*((2/(Lij*Lij*Lij))*dLijdrij - (2/(Uij*Uij*Uij))*dUijdrij) - 0.25*(1/(dist*dist))*log(Lij/Uij) + (Uij/(4*dist*Lij))*((1/Uij)*dLijdrij - (Lij/(Uij*Uij))*dUijdrij) - 0.125*power(S_j_term,2)*((1/(Lij*Lij)) - (1/(Uij*Uij))) + 0.25*((S_j*S_j*offsetRadius_j*offsetRadius_j)/(dist*Uij*Uij*Uij))*dUijdrij - 0.25*((S_j*S_j*offsetRadius_j*offsetRadius_j)/(dist*Lij*Lij*Lij))*dLijdrij + dCijdrij;
  
  //Check Equation (9)
  Real dRidrij = power(bornRad_i,2)*offsetRadius_i*(1-tanh_i*tanh_i)*tanhparam_derv_i*(1/radius_i)*dBTidrij;
  
  topo->atoms[atom1].myGBSA_T->btDerv1[atom2] = dBTidrij;
  topo->atoms[atom1].myGBSA_T->bornRadiusDerivatives[atom2] = dRidrij;
  
  
  Real dBTjdrji = -0.5*dLjidrij*(1/(Lji*Lji)) + 0.5*dUjidrij*(1/(Uji*Uji)) + 0.125*((1/(Uji*Uji)) - (1/(Lji*Lji))) + 0.125*dist*((2/(Lji*Lji*Lji))*dLjidrij - (2/(Uji*Uji*Uji))*dUjidrij) - 0.25*(1/(dist*dist))*log(Lji/Uji) + (Uji/(4*dist*Lji))*((1/Uji)*dLjidrij - (Lji/(Uji*Uji))*dUjidrij) - 0.125*power(S_i_term,2)*((1/(Lji*Lji)) - (1/(Uji*Uji))) + 0.25*((S_i*S_i*offsetRadius_i*offsetRadius_i)/(dist*Uji*Uji*Uji))*dUjidrij - 0.25*((S_i*S_i*offsetRadius_i*offsetRadius_i)/(dist*Lji*Lji*Lji))*dLjidrij + dCjidrij;
  
  Real dRjdrji = power(bornRad_j,2)*offsetRadius_j*(1-tanh_j*tanh_j)*tanhparam_derv_j*(1/radius_j)*dBTjdrji;
  
  topo->atoms[atom2].myGBSA_T->btDerv1[atom1] = dBTjdrji;
  topo->atoms[atom2].myGBSA_T->bornRadiusDerivatives[atom1] = dRjdrji;
  //Check Equation (15)
  force += c1*(c_i*(1/power(bornRad_i,7))*dRidrij*(1/dist) + c_j*(1/power(bornRad_j,7))*dRjdrji*(1/dist));
}

void GBACEForce::accumulateEnergy(ScalarStructure *energies, Real energy) {
  (*energies)[ScalarStructure::COULOMB] += energy;
}

Real GBACEForce::getEnergy(const ScalarStructure *energies) {
  return (*energies)[ScalarStructure::COULOMB];
}

void GBACEForce::postProcess(const GenericTopology *topo, ScalarStructure *energies, Vector3DBlock *forces) {
  const unsigned int atoms = topo->atoms.size();
  
  Real totalSelfEnergy = 0.0;
  
  for( unsigned int i=0; i<atoms; i++){
    totalSelfEnergy += topo->atoms[i].myGBSA_T->ACEPotential;
  }
  
  accumulateEnergy( energies, totalSelfEnergy );
}

void GBACEForce::parallelPostProcess(const GenericTopology *topo, ScalarStructure *energies) {      
  const unsigned int atomnumber = topo->atoms.size();
  
  // Copy self energy count
  Real *selfcount = new Real[ atomnumber ];
  
  //put self energy count into array
  for( unsigned int i = 0; i < atomnumber; i++ ){
    selfcount[i] = (Real)topo->atoms[i].myGBSA_T->ACEPotentialCount;
  }
  
  //sum accross nodes
  Parallel::reduce(selfcount, selfcount + atomnumber); //reduceSlaves only?
  
  //find self energies
  // Copy self energies
  Real *selfenergy = new Real[ atomnumber ];
  
  //put self energy into array
  for( unsigned int i = 0; i < atomnumber; i++ ){
    selfenergy[i] = topo->atoms[i].myGBSA_T->ACEPotential;
  }
  
  //sum accross nodes
  Parallel::reduce(selfenergy, selfenergy + atomnumber); //reduceSlaves only?
  
  //put corrected self energy back
  for( unsigned int i = 0; i < atomnumber; i++ ){
    
    if( selfcount[i] != 0.0 )
      topo->atoms[i].myGBSA_T->ACEPotential = selfenergy[i] / selfcount[i] / (Real)Parallel::getNum();
  }
  
  delete [] selfenergy;
  
  delete [] selfcount;
}

bool GBACEForce::doParallelPostProcess() {
  return true;
}

// Parsing
std::string GBACEForce::getId() {
  return keyword;
}

unsigned int GBACEForce::getParameterSize() {
  return 2;
}

void GBACEForce::getParameters(vector<Parameter> &parameters) const {
     parameters.push_back(
        Parameter("-solvationparam",Value(sigma, ConstraintValueType::NoConstraints()),2.26 / 418.4, Text("solvation parameter")));
        //kJ nm^{-2} -> KCal \AA^{-2} 1/418.4
     parameters.push_back(
        Parameter("-watersphereradius",Value(rho_s, ConstraintValueType::NoConstraints()), 1.4, Text("solvation parameter")));
}

GBACEForce GBACEForce::make(const vector<Value> &values){
   return (GBACEForce(values[0], values[1]));
}
