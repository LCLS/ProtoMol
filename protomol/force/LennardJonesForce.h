/* -*- c++ -*- */
#ifndef LENNARDJONESFORCE_H
#define LENNARDJONESFORCE_H

#include <protomol/base/Report.h>
#include <protomol/config/Parameter.h>
#include <protomol/topology/ExclusionTable.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/topology/LennardJonesParameters.h>
#include <protomol/type/ScalarStructure.h>
#include <string>
#include <vector>

#include "protomol/topology/Atom.h"
#include "protomol/topology/LennardJonesParameterTable.h"
#include "protomol/type/Real.h"
#include "protomol/type/Vector3D.h"

namespace ProtoMol {
class Value;
class Vector3DBlock;
struct Parameter;
}  // namespace ProtoMol

using namespace ProtoMol::Report;

namespace ProtoMol {
  //____ LennardJonesForce
  class LennardJonesForce {
  public:
    enum {DIST_R2 = 1};
    enum {CUTOFF = 0};

  public:
    virtual ~LennardJonesForce() {} // Compiler needs this

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class LennardJonesForce
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    void operator()(Real &energy, Real &force,
                    Real /*distSquared*/, Real rDistSquared, const Vector3D &,
                    const GenericTopology *topo, int atom1, int atom2,
                    ExclusionClass excl) const {
#ifdef DEBUG_LJ_DISTANCE_CHECK
      if (rDistSquared < 1.0 / 0.25)
        report << warning << "LennardJonesForce::operator(): atom " << atom1 <<
        " and " << atom2 << " get closer"
                << " than 0.5 angstroms (=" <<
        sqrt(1.0 / rDistSquared) << ")." << endr;
#endif
      const LennardJonesParameters &
      params = topo->lennardJonesParameters(topo->atoms[atom1].type,
                                            topo->atoms[atom2].type);

      Real A, B;
      if (excl != EXCLUSION_MODIFIED) {
        A = params.A;
        B = params.B;
      } else {
        A = params.A14;
        B = params.B14;
      }

      // Fast LJ computation
      //Real r2   = rDistSquared;
      Real r6 = rDistSquared * rDistSquared * rDistSquared;
      Real r12 = r6 * r6;
      Real r6B = B * r6;
      Real r12A = A * r12;
      energy = r12A - r6B;
      force = 12.0 * r12A * rDistSquared - 6.0 * r6B * rDistSquared;
    }

    static void accumulateEnergy(ScalarStructure *energies, Real energy) {
      (*energies)[ScalarStructure::LENNARDJONES] += energy;
    }

    static Real getEnergy(const ScalarStructure *energies) {
      return (*energies)[ScalarStructure::LENNARDJONES];
    }
    
    virtual void preProcess(const GenericTopology *apptopo, const Vector3DBlock *positions) {
		  
	  }
    
    static void postProcess(const GenericTopology *topo, ScalarStructure *energies, Vector3DBlock *forces){
		  
	  }

    static void parallelPostProcess(const GenericTopology *topo, ScalarStructure *energies) {
      
    }

    //do parallel post process?
    static bool doParallelPostProcess() {
      return false;
    }

    // Parsing
    static std::string getId() {return keyword;}
    static unsigned int getParameterSize() {return 0;}
    void getParameters(std::vector<Parameter> &) const {}

    static LennardJonesForce make(const std::vector<Value> &) {
      return LennardJonesForce();
    }

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    static const std::string keyword;
  private:
  };
  //____ INLINES
}

#endif /* LENNARDJONESFORCE_H */
