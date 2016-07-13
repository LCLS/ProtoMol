#include <math.h>
#include <protomol/base/Report.h>
#include <protomol/force/GB/GBBornRadii.h>
#include <protomol/parallel/Parallel.h>

#include "protomol/topology/Atom.h"
#include "protomol/topology/ExclusionTable.h"
#include "protomol/topology/GenericTopology.h"
#include "protomol/topology/Molecule.h"

namespace ProtoMol {
class ScalarStructure;
class Value;
class Vector3DBlock;
struct Parameter;
}  // namespace ProtoMol

using namespace ProtoMol;
using namespace ProtoMol::Report;

const string GBBornRadii::keyword("GBBornRadii");

GBBornRadii::GBBornRadii(){

}

void GBBornRadii::operator()(Real &energy, Real &force,
                Real distSquared, Real rDistSquared, const Vector3D &diff,
                const GenericTopology *topo,
                int atom1, int atom2, ExclusionClass excl) const {
  
  Real one = (Real) 1.0 ;
  Real two = (Real) 2.0 ;
  Real four = (Real) 4.0;
  
  if (!topo->doGBSAOpenMM) {
    report << error <<"GBBornBurialTerm require GBSA set"<<endr;
  }
  
  energy = 0;
  
  force = 0;
  
  // If either molecule belongs to a water, do nothing.
  // Won't happen in most simulations, but could in the 
  // case of comparing forces.
  if (topo->molecules[topo->atoms[atom1].molecule].water ||
      topo->molecules[topo->atoms[atom2].molecule].water)
    return;
  
  //r_ij
  Real dist = sqrt(distSquared);
  
  //store pairwise distances in some data structure for future use
  topo->atoms[atom1].myGBSA_T->distij[atom2] = dist;
  topo->atoms[atom2].myGBSA_T->distij[atom1] = dist;
  
  Real radius_i = topo->atoms[atom1].myGBSA_T->vanDerWaalRadius;
  Real radius_j = topo->atoms[atom2].myGBSA_T->vanDerWaalRadius;
  
  //offset radii ({\tilde{\rho}_{j}}) Equation (4) in the writeup
  Real offsetRadius_i = radius_i - topo->atoms[atom1].myGBSA_T->offsetRadius;
  Real offsetRadius_j = radius_j - topo->atoms[atom2].myGBSA_T->offsetRadius;
  
  //Scaling factors
  Real S_i = topo->atoms[atom1].myGBSA_T->scalingFactor;
  Real S_j = topo->atoms[atom2].myGBSA_T->scalingFactor;
  
  Real Lij, Uij, Cij;
  
  //Equations (6-8)
  if (offsetRadius_i >=  dist + S_j*offsetRadius_j) {
    Lij = one;
    Uij = one;
    
  }else {
    Lij = (offsetRadius_i > fabs(dist - S_j*offsetRadius_j)) ? offsetRadius_i : fabs(dist - S_j*offsetRadius_j);
    
    Uij = dist + S_j*offsetRadius_j;
    
  }
  if (offsetRadius_i < offsetRadius_j*S_j- dist) {
    Cij = two*(one/offsetRadius_i - one/Lij);
  }else {
    Cij = 0;
  }
  
  //store Lvalues Uvalues in an array so that we can use later
  topo->atoms[atom1].myGBSA_T->Lvalues[atom2] = Lij;
  topo->atoms[atom1].myGBSA_T->Uvalues[atom2] = Uij;
  
  
  Real invLij = one/Lij;
  Real invUij = one/Uij;
  
  Real invLij2 = invLij*invLij;
  Real invUij2 = invUij*invUij;
  
  Real ratio_i = log(Lij/Uij);
  
  
  //add this to the burial term for atom i (atom1) (see Equation (5) in the writeup)
  Real term_i = (invLij - invUij) + (dist/four)*(invUij2 - invLij2) + (one/(two*dist))*(Real)ratio_i + ((S_j*S_j*offsetRadius_j*offsetRadius_j)/(four*dist))*(invLij2 - invUij2) + Cij;
  topo->atoms[atom1].myGBSA_T->burialTerm += term_i;
  
  // These are required for jith term. Note that Lij not necessary equal to Lji. Same is true for Uij and Cij.
  // These quantities may not be symmetric. I need these for calculating R_j.
  Real Lji, Uji, Cji;
  
  if (offsetRadius_j >=  dist + S_i*offsetRadius_i) {
    Lji = one;
    Uji = one;
  }else {
    Lji = (offsetRadius_j > fabs(dist - S_i*offsetRadius_i)) ? offsetRadius_j : fabs(dist - S_i*offsetRadius_i);         
    Uji = dist + S_i*offsetRadius_i;
  }
  
  //Store calculated L and U values
  topo->atoms[atom2].myGBSA_T->Lvalues[atom1] = Lji;
  topo->atoms[atom2].myGBSA_T->Uvalues[atom1] = Uji;
  
  if (offsetRadius_j < offsetRadius_i*S_i - dist) {
    Cji = two*(one/offsetRadius_j - one/Lji);
    
  }else {
    Cji = 0;
  }
  
  Real invLji = one/Lji;
  Real invUji = one/Uji;
  
  Real invLji2 = invLji*invLji;
  Real invUji2 = invUji*invUji;
  
  Real ratio_j = log(Lji/Uji);
  
  //add this to the burial term for atom j (atom2) (see Equation (5) in the writeup)
  Real term_j = (invLji - invUji) + (dist/four)*(invUji2 - invLji2) + (one/(two*dist))*(Real)ratio_j + ((S_i*S_i*offsetRadius_i*offsetRadius_i)/(four*dist))*(invLji2 - invUji2) + Cji;
  topo->atoms[atom2].myGBSA_T->burialTerm += term_j;
}

void GBBornRadii::accumulateEnergy(ScalarStructure *energies, Real energy) {
  
}

Real GBBornRadii::getEnergy(const ScalarStructure *energies) {
  return 0;
}

//now initialize GB structs here
void GBBornRadii::preProcess(const GenericTopology *apptopo, const Vector3DBlock *positions) {
  const unsigned int atoms = apptopo->atoms.size();
  for(unsigned int i=0;i<atoms;i++)
    apptopo->atoms[i].myGBSA_T->preForce();
}

void GBBornRadii::postProcess(const GenericTopology *topo, ScalarStructure *energies, Vector3DBlock *forces) {
  const unsigned int atomnumber = topo->atoms.size();
  
  //~~~~calculate born radius from burial term~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  for(unsigned int i=0; i<atomnumber; i++){
    Real radius_i = topo->atoms[i].myGBSA_T->vanDerWaalRadius;
    
    //offset radii ({\tilde{\rho}_{j}})
    Real offsetRadius_i = radius_i - topo->atoms[i].myGBSA_T->offsetRadius;
    
    Real psi_i;
    
    Real tanhparam_i;
    
    //calculate born radius
    //set flag to true
    
    //Equation (3) and (5)
    topo->atoms[i].myGBSA_T->PsiValue = 0.5*topo->atoms[i].myGBSA_T->burialTerm*offsetRadius_i;
    psi_i = topo->atoms[i].myGBSA_T->PsiValue;
    
    //part of Equation (1)
    tanhparam_i = topo->alphaObc*psi_i - topo->betaObc*psi_i*psi_i + topo->gammaObc*psi_i*psi_i*psi_i;
    //Second part of Equation (1)
    Real invBornRad_i = (1/offsetRadius_i) - (1/radius_i)*tanh(tanhparam_i);
    topo->atoms[i].myGBSA_T->bornRad = 1/invBornRad_i;
    topo->atoms[i].myGBSA_T->doneCalculateBornRadius = true;
  }
}

void GBBornRadii::parallelPostProcess(const GenericTopology *topo, ScalarStructure *energies) {
  const unsigned int atomnumber = topo->atoms.size();
  
  // Copy Radii
  Real *burial = new Real[ atomnumber ];
  
  //put radii into array
  for( unsigned int i = 0; i < atomnumber; i++ ){
    burial[i] = topo->atoms[i].myGBSA_T->burialTerm;
  }
  
  //sum accross nodes
  Parallel::reduce(burial, burial + atomnumber); //reduceSlaves only?
  
  //put radii back
  for( unsigned int i = 0; i < atomnumber; i++ ){
    topo->atoms[i].myGBSA_T->burialTerm = burial[i];
  }
  
  delete [] burial;
}

bool GBBornRadii::doParallelPostProcess() {
  return true;
}

// Parsing
std::string GBBornRadii::getId() {
  return keyword;
}

unsigned int GBBornRadii::getParameterSize() {
  return 0;
}

void GBBornRadii::getParameters(std::vector<Parameter> &) const {

}

GBBornRadii GBBornRadii::make(const std::vector<Value> &) {
  return GBBornRadii();
}
