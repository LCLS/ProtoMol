/* -*- c++ -*- */
#ifndef COULOMBSCPISMFORCE_H
#define COULOMBSCPISMFORCE_H

#include <math.h>
#include <protomol/base/PMConstants.h>
#include <protomol/base/Report.h>
#include <protomol/config/Parameter.h>
#include <protomol/topology/ExclusionTable.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/type/ScalarStructure.h>
#include <iostream>
#include <string>
#include <vector>

#include "protomol/topology/Atom.h"
#include "protomol/topology/Molecule.h"
#include "protomol/type/Real.h"
#include "protomol/type/Vector3D.h"

namespace ProtoMol {
class Value;
class Vector3DBlock;
struct Parameter;
}  // namespace ProtoMol

using namespace std;

using namespace ProtoMol::Report;

namespace ProtoMol {
  //____ CoulombForce
  class CoulombSCPISMForce {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Electrostatic interaction term for SCP ISM for MD of
    // Hassan et al. PROTEINS (51):109 (2003)
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    enum {DIST_R2 = 1};
    enum {CUTOFF = 0};

  public:
    CoulombSCPISMForce();
    CoulombSCPISMForce(Real Dval);
    virtual ~CoulombSCPISMForce() {} // Compiler needs this

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class CoulombSCPISForce
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    void operator()(Real &energy, Real &force,
                    Real /*distSquared*/, Real rDistSquared, const Vector3D &,
                    const GenericTopology *topo,
                    int atom1, int atom2, ExclusionClass excl) const {

      //SCPISM valid?
      if(!topo->doSCPISM)
        report << error << "CoulombSCPISMForce requires SCPISM parameters." << endr;
      // If either molecule belongs to a water, do nothing.
      // Won't happen in most simulations, but could in the 
      // case of comparing forces.
      if (topo->molecules[topo->atoms[atom1].molecule].water ||
	      topo->molecules[topo->atoms[atom2].molecule].water)
	        return;

      //std::cout << "EPS: " << EPS << std::endl;
      //std::cout << "D: " << D << std::endl;
      //std::cout << "S: " << S << std::endl;

      Real rDist = sqrt(rDistSquared);
      Real Dist = 1.0 / rDist;

      // Screening term for interaction energy calculation
      Real alpha_ij = topo->atoms[atom1].mySCPISM_A->sqrtalphaSCPISM *
                      topo->atoms[atom2].mySCPISM_A->sqrtalphaSCPISM;
      Real Dp1 = D + 1.0;
      Real rDp1 = 1.0 / Dp1;
      Real Dm1 = D - 1.0;
      Real k = Dm1 * 0.5;
      Real rDiElecEXP = k * exp(-alpha_ij * Dist);
      Real DiElec = Dp1 / (1.0 + rDiElecEXP) - 1.0;
      Real rDiElec = 1.0 / DiElec;

      energy = topo->atoms[atom1].scaledCharge *
               topo->atoms[atom2].scaledCharge *
               rDist * rDiElec *
               (excl == EXCLUSION_MODIFIED ? topo->coulombScalingFactor : 1);

      // reciprocal force term calculation
      Real drDiElec = alpha_ij * rDp1 * (1.0 + DiElec) * (D - DiElec);
      force = (rDiElec * drDiElec + rDist) * energy * rDist;
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
    static unsigned int getParameterSize() {return 1;}
    static CoulombSCPISMForce make(const std::vector<Value> &);

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
    Real D;
  };

  //____ INLINES
}
#endif /* COULOMBSCPISMFORCE_H */
