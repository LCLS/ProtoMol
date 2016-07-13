#include <protomol/force/GB/GBPartialSum.h>
#include <protomol/parallel/Parallel.h>
#include <cmath>

#include "protomol/topology/Atom.h"
#include "protomol/topology/ExclusionTable.h"
#include "protomol/topology/GenericTopology.h"

namespace ProtoMol {
class ScalarStructure;
class Value;
class Vector3DBlock;
struct Parameter;
}  // namespace ProtoMol

using namespace ProtoMol;

const string GBPartialSum::keyword("GBPartialSum");

GBPartialSum::GBPartialSum(){}

void GBPartialSum::operator()(Real &energy, Real &force,
                Real distSquared, Real rDistSquared, const Vector3D &diff,
                const GenericTopology *topo,
                int atom1, int atom2, ExclusionClass excl) const {
  
  Real bornRad_i = topo->atoms[atom1].myGBSA_T->bornRad;
  Real bornRad_l = topo->atoms[atom2].myGBSA_T->bornRad;
  
  Real scaledCharge_i = topo->atoms[atom1].scaledCharge;
  Real scaledCharge_l = topo->atoms[atom2].scaledCharge;
  
  Real ril = std::sqrt(distSquared);
  
  Real expterm = std::exp( -(ril*ril)/(4.0*bornRad_i*bornRad_l) );
  topo->atoms[atom1].myGBSA_T->expTerm[atom2] = expterm;
  topo->atoms[atom2].myGBSA_T->expTerm[atom1] = expterm;
  
  Real filGB = std::sqrt(ril*ril + bornRad_i*bornRad_l*expterm);
  topo->atoms[atom1].myGBSA_T->filTerm[atom2] = filGB;
  topo->atoms[atom2].myGBSA_T->filTerm[atom1] = filGB;
  
  Real part = scaledCharge_i*scaledCharge_l*(1/(filGB*filGB))*0.5*(1/filGB)*expterm;
  
  Real aTerm = part*(bornRad_l + (ril*ril)/(4.0*bornRad_i));
  Real bTerm = part*(bornRad_i + (ril*ril)/(4.0*bornRad_l));
  
  topo->atoms[atom1].myGBSA_T->partialGBForceTerms += aTerm;
  topo->atoms[atom2].myGBSA_T->partialGBForceTerms += bTerm;
  
  topo->atoms[atom1].myGBSA_T->partialTerm[atom2] = aTerm;
  topo->atoms[atom2].myGBSA_T->partialTerm[atom1] = bTerm;
}

void GBPartialSum::accumulateEnergy(ScalarStructure *energies, Real energy) {
  
}

Real GBPartialSum::getEnergy(const ScalarStructure *energies) {
  return 0;
}

void GBPartialSum::postProcess(const GenericTopology *topo, ScalarStructure *energies, Vector3DBlock *forces) {
  const unsigned int atomnumber = topo->atoms.size();
  
  for( unsigned int i=0; i<atomnumber; i++)
    topo->atoms[i].myGBSA_T->havePartialGBForceTerms = true;
}

void GBPartialSum::parallelPostProcess(const GenericTopology *topo, ScalarStructure *energies) {
  const unsigned int atomnumber = topo->atoms.size();
  
  // Copy Radii
  Real *partial = new Real[ atomnumber ];
  
  //put radii into array
  for( unsigned int i = 0; i < atomnumber; i++ ){
    partial[i] = topo->atoms[i].myGBSA_T->partialGBForceTerms;
  }
  
  //sum accross nodes
  Parallel::reduce(partial, partial + atomnumber); //reduceSlaves only?
  
  //put radii back
  for( unsigned int i = 0; i < atomnumber; i++ ){
    topo->atoms[i].myGBSA_T->partialGBForceTerms = partial[i];
  }
  
  delete [] partial;
}

bool GBPartialSum::doParallelPostProcess() {
  return true;
}

// Parsing
std::string GBPartialSum::getId() {
  return keyword;
}

unsigned int GBPartialSum::getParameterSize() {
  return 0;
}

void GBPartialSum::getParameters(std::vector<Parameter> &) const {

}

GBPartialSum GBPartialSum::make(const std::vector<Value> &) {
  return GBPartialSum();
}
