/* -*- c++ -*- */
#ifndef GBBORNRADII_H
#define GBBORNRADII_H

#include <protomol/topology/ExclusionTable.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/type/ScalarStructure.h>
#include <string>
#include <vector>

#include "protomol/type/Real.h"
#include "protomol/type/Vector3D.h"

namespace ProtoMol {
class GenericTopology;
class ScalarStructure;
class Value;
class Vector3DBlock;
struct Parameter;

  class GBBornRadii {
    public:
      enum {DIST_R2 = 1};
      enum {CUTOFF = 0};
    public:
      GBBornRadii();
      void operator()(Real &energy, Real &force, Real distSquared, Real rDistSquared, 
                      const Vector3D &diff, const GenericTopology *topo, int atom1, 
                      int atom2, ExclusionClass excl) const;
      
      static void accumulateEnergy(ScalarStructure *energies, Real energy);
      static Real getEnergy(const ScalarStructure *energies);
      
      static void preProcess(const GenericTopology *apptopo, const Vector3DBlock *positions);
      static void postProcess(const GenericTopology *topo, ScalarStructure *energies, Vector3DBlock *forces);
      static void parallelPostProcess(const GenericTopology *topo, ScalarStructure *energies);
      static bool doParallelPostProcess();
      
      // Parsing
      static std::string getId();
      static unsigned int getParameterSize();
      void getParameters(std::vector<Parameter> &) const;
      
      static GBBornRadii make(const std::vector<Value> &);
    public:
      static const std::string keyword;
  };
}

#endif
