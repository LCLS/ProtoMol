#include <protomol/base/Report.h>
#include <protomol/force/GB/GBForce.h>
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

const string GBForce::keyword("GBForce");

GBForce::GBForce() : soluteDielec(1.0), solventDielec(80.0) {}
GBForce::GBForce(Real solute_d, Real solvent_d) : soluteDielec(solute_d), solventDielec(solvent_d) {}

void GBForce::operator()(Real &energy, Real &force, Real distSquared,
                Real rDistSquared, const Vector3D &,
                const GenericTopology *topo, int atom1, int atom2,
                ExclusionClass excl) const {
  
  if(!topo->doGBSAOpenMM)
    report << error << "GBForce requires GB  parameters." << endr;
  
  //r_{ij}
  Real dist = sqrt(distSquared);
  
  Real radius_i = topo->atoms[atom1].myGBSA_T->vanDerWaalRadius;
  Real radius_j = topo->atoms[atom2].myGBSA_T->vanDerWaalRadius;
  
  //offset radii ({\tilde{\rho}_{j}})
  Real offsetRadius_i = radius_i - topo->atoms[atom1].myGBSA_T->offsetRadius;
  Real offsetRadius_j = radius_j - topo->atoms[atom2].myGBSA_T->offsetRadius;
  
  Real bornRad_i, bornRad_j;
  Real psi_i, psi_j;
  Real tanhparam_i, tanhparam_j;
  
  
  psi_i = topo->atoms[atom1].myGBSA_T->PsiValue;
  tanhparam_i = topo->alphaObc*psi_i - topo->betaObc*psi_i*psi_i + topo->gammaObc*psi_i*psi_i*psi_i;
  psi_j = topo->atoms[atom2].myGBSA_T->PsiValue;
  tanhparam_j = topo->alphaObc*psi_j - topo->betaObc*psi_j*psi_j + topo->gammaObc*psi_j*psi_j*psi_j;
  
  bornRad_i = topo->atoms[atom1].myGBSA_T->bornRad;
  bornRad_j = topo->atoms[atom2].myGBSA_T->bornRad;
  
  //Equation (17)
  Real expterm = topo->atoms[atom1].myGBSA_T->expTerm[atom2]; //std::exp( -(dist*dist)/(4*bornRad_i*bornRad_j) );
  Real fGB = topo->atoms[atom1].myGBSA_T->filTerm[atom2];//sqrt(dist*dist + bornRad_i*bornRad_j*expterm);
  
  Real scaledCharge_i = topo->atoms[atom1].scaledCharge;
  Real scaledCharge_j = topo->atoms[atom2].scaledCharge;
  
  //Equation (16)
  energy = -(scaledCharge_i*scaledCharge_j)*(1/fGB)*((1/soluteDielec) - (1/solventDielec));
  
  //self terms (Equation (18))
  if (!topo->atoms[atom1].myGBSA_T->selfEnergyCount) {
    topo->atoms[atom1].myGBSA_T->selfEnergy = -0.5 * (scaledCharge_i*scaledCharge_i)*(1/bornRad_i)*((1/soluteDielec) - (1/solventDielec));
    topo->atoms[atom1].myGBSA_T->selfEnergyCount = 1;
  }
  
  if (!topo->atoms[atom2].myGBSA_T->selfEnergyCount) {
    topo->atoms[atom2].myGBSA_T->selfEnergy = -0.5 * (scaledCharge_j*scaledCharge_j)*(1/bornRad_j)*((1/soluteDielec) - (1/solventDielec));
    topo->atoms[atom2].myGBSA_T->selfEnergyCount = 1;
  }
  
  //Scaling factors
  Real S_i = topo->atoms[atom1].myGBSA_T->scalingFactor;
  Real S_j = topo->atoms[atom2].myGBSA_T->scalingFactor;
  
  //It would be nice if we can precalculate this terms.
  Real Lij, Uij;
  
  if (offsetRadius_i >=  dist + S_j*offsetRadius_j) {
    Lij = 1;
    Uij = 1;
  }else {
    Lij =(offsetRadius_i > fabs(dist - S_j*offsetRadius_j)) ? offsetRadius_i : fabs(dist - S_j*offsetRadius_j);
    Uij = dist + S_j*offsetRadius_j;
  }
  
  Real invLij = 1/Lij;
  //Real invUij = 1/Uij;
  
  //Derivatives for calculation of the derivative of the born radii
  Real dLijdrij, dUijdrij, dCijdrij;
  if (offsetRadius_i <= (dist - S_j*offsetRadius_j)) dLijdrij = 1;
  else dLijdrij = 0;
  
  if (offsetRadius_i < (dist + S_j*offsetRadius_j)) dUijdrij = 1;
  else dUijdrij = 0;
  
  if (offsetRadius_i <= (S_j*offsetRadius_j - dist)) dCijdrij = 2*invLij*invLij*dLijdrij;
  else dCijdrij = 0;
  
  //Lik
  Real Lji;// = topo->atoms[atom2].myGBSA_T->Lvalues[atom1];
  //Uki
  Real Uji;// = topo->atoms[atom2].myGBSA_T->Uvalues[atom1];
  
  if (offsetRadius_j >=  dist + S_i*offsetRadius_i) {
    Lji = 1;
    Uji = 1;
  }else {
    Lji = (offsetRadius_j > abs(dist - S_i*offsetRadius_i)) ? offsetRadius_j : fabs(dist - S_i*offsetRadius_i);         
    Uji = dist + S_i*offsetRadius_i;
  }
  
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
  
  
  Real dBTidrij = -0.5*dLijdrij*(1/(Lij*Lij)) + 0.5*dUijdrij*(1/(Uij*Uij)) + 0.125*((1/(Uij*Uij)) - (1/(Lij*Lij))) + 0.125*dist*((2/(Lij*Lij*Lij))*dLijdrij - (2/(Uij*Uij*Uij))*dUijdrij) - 0.25*(1/(dist*dist))*log(Lij/Uij) + (Uij/(4*dist*Lij))*((1/Uij)*dLijdrij - (Lij/(Uij*Uij))*dUijdrij) - 0.125*power(S_j_term,2)*((1/(Lij*Lij)) - (1/(Uij*Uij))) + 0.25*((S_j*S_j*offsetRadius_j*offsetRadius_j)/(dist*Uij*Uij*Uij))*dUijdrij - 0.25*((S_j*S_j*offsetRadius_j*offsetRadius_j)/(dist*Lij*Lij*Lij))*dLijdrij + dCijdrij;
  Real dRidrij = power(bornRad_i,2)*offsetRadius_i*(1-tanh_i*tanh_i)*tanhparam_derv_i*(1/radius_i)*dBTidrij;
  
  Real dBTjdrji = -0.5*dLjidrij*(1/(Lji*Lji)) + 0.5*dUjidrij*(1/(Uji*Uji)) + 0.125*((1/(Uji*Uji)) - (1/(Lji*Lji))) + 0.125*dist*((2/(Lji*Lji*Lji))*dLjidrij - (2/(Uji*Uji*Uji))*dUjidrij) - 0.25*(1/(dist*dist))*log(Lji/Uji) + (Uji/(4*dist*Lji))*((1/Uji)*dLjidrij - (Lji/(Uji*Uji))*dUjidrij) - 0.125*power(S_i_term,2)*((1/(Lji*Lji)) - (1/(Uji*Uji))) + 0.25*((S_i*S_i*offsetRadius_i*offsetRadius_i)/(dist*Uji*Uji*Uji))*dUjidrij - 0.25*((S_i*S_i*offsetRadius_i*offsetRadius_i)/(dist*Lji*Lji*Lji))*dLjidrij + dCjidrij;
  
  Real dRjdrji = power(bornRad_j,2)*offsetRadius_j*(1-tanh_j*tanh_j)*tanhparam_derv_j*(1/radius_j)*dBTjdrji;
  
  //force due to pairwise i-j term
  //Check Equation (19-20). Next line only finds the pairwise terms.
  force -= scaledCharge_i*scaledCharge_j*(1/(fGB*fGB))*0.5*(1/fGB)*((2*dist - 0.5*expterm*dist) + expterm*dRidrij*(bornRad_j + (dist*dist)/(4*bornRad_i)) + expterm*dRjdrji*(bornRad_i + (dist*dist)/(4*bornRad_j)))*(1/dist);
  
  // calculation of self (i-i and j-j) terms
  force -= 0.5 * scaledCharge_i * scaledCharge_i * dRidrij / (bornRad_i * bornRad_i) / dist;
  force -= 0.5 * scaledCharge_j * scaledCharge_j * dRjdrji / (bornRad_j * bornRad_j) / dist;
  
  //new N^2 handeling of interaction with other atoms
  if (!topo->atoms[atom1].myGBSA_T->havePartialGBForceTerms) {
    //if here we should not be //lel
    if( Parallel::isParallel() ) report << error << "GBForce: Parallel execution requires force GBPartialSum!" << endr;;
    
    topo->atoms[atom1].myGBSA_T->partialGBForceTerms = Force_i_term(topo, atom1);
    topo->atoms[atom1].myGBSA_T->havePartialGBForceTerms = true;
  }
  
  force -= (topo->atoms[atom1].myGBSA_T->partialGBForceTerms - topo->atoms[atom1].myGBSA_T->partialTerm[atom2]) *(dRidrij/dist);
  
  if (!topo->atoms[atom2].myGBSA_T->havePartialGBForceTerms) {
    //if here we should not be //lel
    if( Parallel::isParallel() ) report << error << "GBForce: Parallel execution requires force GBPartialSum!" << endr;;
    
    topo->atoms[atom2].myGBSA_T->partialGBForceTerms = Force_i_term(topo, atom2);
    topo->atoms[atom2].myGBSA_T->havePartialGBForceTerms = true;
  }
  
  force -= (topo->atoms[atom2].myGBSA_T->partialGBForceTerms - topo->atoms[atom2].myGBSA_T->partialTerm[atom1])*(dRjdrji/dist);
  //end
  
  force *= ((1/soluteDielec) - (1/solventDielec));
  
}

//estimate the force term for the sum over k,l where k=i,j and l {\neq} j if 
//k=i and l {\neq}i if k=j
Real GBForce::Force_i_term(const GenericTopology *topo, int atom1) const{
  Real bornRad_i = topo->atoms[atom1].myGBSA_T->bornRad;
  Real scaledCharge_i = topo->atoms[atom1].scaledCharge;
  
  Real force = 0;
  for (unsigned int l = 0; l < topo->atoms.size(); l++) {
    if( l != (unsigned int) atom1 ){
      Real ril = topo->atoms[atom1].myGBSA_T->distij[l];
    
      Real bornRad_l = topo->atoms[l].myGBSA_T->bornRad;
      Real scaledCharge_l = topo->atoms[l].scaledCharge;
    
      Real expterm = std::exp( -(ril*ril)/(4.0*bornRad_i*bornRad_l) );
      Real filGB = sqrt(ril*ril + bornRad_i*bornRad_l*expterm);
      
      force += scaledCharge_i*scaledCharge_l*(1/(filGB*filGB))*0.5*(1/filGB)*expterm*(bornRad_l + (ril*ril)/(4.0*bornRad_i));
    }
  }
  
  return force; 
  
}

//estimate the force term for the sum over k,l where k=i,j and l {\neq} j if 
//k=i and l {\neq}i if k=j
// modified to calculate the force only between two atoms.  will be used in
// optimization
Real GBForce::Force_i_j_term(const GenericTopology *topo, const int atom1, const int atom2, const Real ril) const{
  
  const Real bornRad_i = topo->atoms[atom1].myGBSA_T->bornRad;
  const Real scaledCharge_i = topo->atoms[atom1].scaledCharge;
  
  const Real bornRad_l = topo->atoms[atom2].myGBSA_T->bornRad;
  const Real scaledCharge_l = topo->atoms[atom2].scaledCharge;
  
  const Real expterm = std::exp( -(ril*ril)/(4.0*bornRad_i*bornRad_l) );
  const Real filGB = sqrt(ril*ril + bornRad_i*bornRad_l*expterm);
  
  return scaledCharge_i*scaledCharge_l*(1/(filGB*filGB))*0.5*(1/filGB)*expterm*(bornRad_l + (ril*ril)/(4.0*bornRad_i));
  
}

void GBForce::accumulateEnergy(ScalarStructure *energies, Real energy) {
  (*energies)[ScalarStructure::COULOMB] += energy;
}

Real GBForce::getEnergy(const ScalarStructure *energies) {
  return (*energies)[ScalarStructure::COULOMB];
}

void GBForce::postProcess(const GenericTopology *topo, ScalarStructure *energies, Vector3DBlock *forces) {
  
  const unsigned int atoms = topo->atoms.size();
  
  Real totalSelfEnergy = 0.0;
  
  for( unsigned int i=0; i<atoms; i++){
    totalSelfEnergy += topo->atoms[i].myGBSA_T->selfEnergy;
  }
  
  accumulateEnergy( energies, totalSelfEnergy );
  
}

void GBForce::parallelPostProcess(const GenericTopology *topo, ScalarStructure *energies) {
  
  //find self energy count
  
  const unsigned int atomnumber = topo->atoms.size();
  
  // Copy self energy count
  Real *selfcount = new Real[ atomnumber ];
  
  //put self energy count into array
  for( unsigned int i = 0; i < atomnumber; i++ ){
    selfcount[i] = (Real)topo->atoms[i].myGBSA_T->selfEnergyCount;
  }
  
  //sum accross nodes
  Parallel::reduce(selfcount, selfcount + atomnumber); //reduceSlaves only?
  
  //find self energies
  // Copy self energies
  Real *selfenergy = new Real[ atomnumber ];
  
  //put self energy into array
  for( unsigned int i = 0; i < atomnumber; i++ ){
    selfenergy[i] = topo->atoms[i].myGBSA_T->selfEnergy;
  }
  
  //sum accross nodes
  Parallel::reduce(selfenergy, selfenergy + atomnumber); //reduceSlaves only?
  
  //put corrected self energy back
  for( unsigned int i = 0; i < atomnumber; i++ ){
    
    if( selfcount[i] != 0.0 )
      topo->atoms[i].myGBSA_T->selfEnergy = selfenergy[i] / selfcount[i] / (Real)Parallel::getNum();
  }
  
  delete [] selfenergy;
  
  delete [] selfcount;
}

bool GBForce::doParallelPostProcess() {
  return true;
}

// Parsing
std::string GBForce::getId() {
  return keyword;
}

unsigned int GBForce::getParameterSize() {
  return 2;
}

void GBForce::getParameters(vector<Parameter> &parameters) const {
   parameters.push_back
     (Parameter("-soluteDielec", Value(soluteDielec, ConstraintValueType::NoConstraints()), 1.0, Text("Solute Dielectric")));
   parameters.push_back
     (Parameter("-solventDielec", Value(solventDielec, ConstraintValueType::NoConstraints()), 80.0, Text("Solvent Dielectric")));
}

GBForce GBForce::make(const vector<Value> &values) {
   return (GBForce(values[0],values[1]));
}
