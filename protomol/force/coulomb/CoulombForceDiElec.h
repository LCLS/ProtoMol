/* -*- c++ -*- */
#ifndef COULOMBFORCEDIELEC_H
#define COULOMBFORCEDIELEC_H

#include <math.h>
#include <protomol/config/Parameter.h>
#include <protomol/topology/ExclusionTable.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/type/ScalarStructure.h>
#include <string>
#include <vector>

#include "protomol/topology/Atom.h"
#include "protomol/type/Real.h"
#include "protomol/type/Vector3D.h"

namespace ProtoMol {
  //____ CoulombForce
class Value;
class Vector3DBlock;
struct Parameter;

  class CoulombForceDiElec {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // This uses the weighted charges on each atom.
    // A DiElectric value is added to the electrostatic computation according
    // to Shen and Freed 2002 "All-atom fast protein folding simulations:The
    // Villin Headpiece" PROTEINS:Structure, Function, and Genetics 49:439-445
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  public:
    enum {DIST_R2 = 1};
    enum {CUTOFF = 0};
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    // Default constructor
    CoulombForceDiElec();
    // Constructor with parameters
    CoulombForceDiElec(Real EPSval, Real Dval, Real Sval);
    // Dielectric function constant terms D = 53.0 and S = 0.25
    // To change the value of eps(0)= 1 given with DiElecCom
    // Set X in DiElecCom = (D-X)*0.5 to your chosen value 2,3,4,etc.. 
    // which yields eps(0)=2,3,4,etc...
    // eps(0)=1 is insufficient to model solvent polarization effects between
    // near molecules
    // Try 2,3,or 4
    virtual ~CoulombForceDiElec() {} // Compiler needs this

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class CoulombForce
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    void operator()(Real &energy,
                    Real &force,
                    Real /*distSquared*/, Real rDistSquared, const Vector3D &,
                    const GenericTopology *topo,
                    int atom1, int atom2, ExclusionClass excl) const {
      Real rDist = sqrt(rDistSquared);
      Real Dist = 1.0 / rDist;

      Real DiElecCom = (D - EPS) * 0.5;

      // Dielectric term for energy calculation
      Real rDiElecEXP = exp(-S * Dist);
      Real rDiElecP1 = DiElecCom * (S * S * Dist * Dist + 2 * S * Dist + 2);
      Real rDiElec = 1.0 / (D - (rDiElecP1 * rDiElecEXP));

      energy = topo->atoms[atom1].scaledCharge *
               topo->atoms[atom2].scaledCharge *
               rDist * rDiElec;

      if (excl == EXCLUSION_MODIFIED)
        energy *= topo->coulombScalingFactor;

      // reciprocal force term calculation includes dielectric derivative
      Real drDistrDiElecPoly = S * S * Dist * Dist * Dist + 2 * S * Dist *
                               Dist + 2 * Dist;
      Real drDistrDiElecP1 = DiElecCom * rDiElecEXP *
                             (3 * S * S * Dist * Dist + 4 * S * Dist + 2);
      Real drDistrDiElecP2 = DiElecCom * rDiElecEXP * S * drDistrDiElecPoly;
      Real drDistrDiElecP3 = D * Dist - DiElecCom * rDiElecEXP *
                             drDistrDiElecPoly;

      Real drDistrDiElec =
        (D - drDistrDiElecP1 + drDistrDiElecP2) / 
        (drDistrDiElecP3 * drDistrDiElecP3);
      force = energy / rDiElec * drDistrDiElec;

      // no negative term - magnitude only
      // additional rDist term normalized force for later mult by directional
      // vec components
    }

    static void accumulateEnergy(ScalarStructure *energies, Real energy) {
      (*energies)[ScalarStructure::COULOMB] += energy;
    }

    static Real getEnergy(const ScalarStructure *energies) {
      return (*energies)[ScalarStructure::COULOMB];
    }
    
    virtual void preProcess(const GenericTopology *apptopo, const Vector3DBlock *positions) {
		  
	  }
    
    static void postProcess(const GenericTopology *topo, ScalarStructure *energies, Vector3DBlock *forces) {
      
    }

    static void parallelPostProcess(const GenericTopology *topo, ScalarStructure *energies) {
      
    }

    //do parallel post process?
    static bool doParallelPostProcess() {
      return false;
    }
    
    // Parsing
    static std::string getId() {return keyword;}
    void getParameters(std::vector<Parameter> &) const;
    static unsigned int getParameterSize() {return 3;}
    static CoulombForceDiElec make(const std::vector<Value> &);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    static const std::string keyword;
  private:
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    Real EPS;
    Real D;
    Real S;
  };
}
#endif /* COULOMBFORCEDIELEC_H */
