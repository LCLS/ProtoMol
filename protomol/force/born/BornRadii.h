/* -*- c++ -*- */
#ifndef BORNRADII_H
#define BORNRADII_H

#include <math.h>
#include <protomol/base/Report.h>
#include <protomol/config/Parameter.h>
#include <protomol/force/born/BornSwitch.h>
#include <protomol/parallel/Parallel.h>
#include <protomol/topology/ExclusionTable.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/type/Real.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/type/Vector3D.h>
#include <string>
#include <vector>

#include "protomol/topology/Atom.h"
#include "protomol/topology/AtomType.h"
#include "protomol/topology/Molecule.h"

namespace ProtoMol {
class ScalarStructure;
class Value;
class Vector3DBlock;
struct Parameter;
}  // namespace ProtoMol

using namespace ProtoMol::Report;

namespace ProtoMol {

  //____ BornRadii
  class BornRadii {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Born radii calculation.
    // Equation numbers from "Notes on SCPISM" by C.R.Sweet, based on Hassan's
    // papers.
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  public:
    enum {DIST_R2 = 1};
    enum {CUTOFF = 0};

    BornRadii() : myBornSwitch(3){};
    BornRadii(int bsw) : myBornSwitch(bsw){};
    virtual ~BornRadii() {} // Compiler needs this

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // New methods of class BornRadii
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    void operator()(Real &energy, Real &force,
                    Real distSquared, Real rDistSquared, const Vector3D &diff,
                    const GenericTopology *topo,
                    int atom1, int atom2, ExclusionClass excl) const {
      
      //SCPISM valid?
      if(!topo->doSCPISM)
        report << error << "BornRadii requires SCPISM parameters." << endr;
      // If either molecule belongs to a water, do nothing.
      // Won't happen in most simulations, but could in the 
      // case of comparing forces.
      if (topo->molecules[topo->atoms[atom1].molecule].water ||
	      topo->molecules[topo->atoms[atom2].molecule].water)
	        return;

      // Distance and shift function
      Real dist = sqrt(distSquared);

      BornSwitch bSw(myBornSwitch);  // Always quartic now  6/1/2008

      Real f_ij = bSw.switchValue(dist);

      // Atom 1 variables
      int type1 = topo->atoms[atom1].type;
      int type2 = topo->atoms[atom2].type;

      //**********************************
      // Part of Eq. (1)
      Real C_i = topo->atomTypes[type1].mySCPISM_T->C_i;
      Real eta_i = topo->atoms[atom1].mySCPISM_A->eta;
      topo->atoms[atom1].mySCPISM_A->bornRadius += eta_i * f_ij * exp(-C_i * dist);
      //**********************************

      // If atom 1 is a polar H+, accumulate the derivative dR.
      if (topo->atomTypes[type1].mySCPISM_T->isHbonded == PH &&
          topo->atomTypes[type2].mySCPISM_T->isHbonded == PA &&
          (topo->atoms[atom1].residue_seq !=
           topo->atoms[atom2].residue_seq)) { // Polar
        Real E_i = 0.80; // Currently all polar H+ have this value
        Real g_i;
        if (topo->atoms[atom1].name == "HN" &&
            topo->atoms[atom2].name == "O")
          g_i = -0.378;
        else
          g_i = topo->atomTypes[type1].mySCPISM_T->g_i;
        Real g_j = topo->atomTypes[type2].mySCPISM_T->g_i;
        // Eq. (2)
        topo->atoms[atom1].mySCPISM_A->bornRadius += g_i * g_j * f_ij * exp(-E_i * dist);
      }

      //**********************************
      // Part of Eq. (1)      
      Real C_j = topo->atomTypes[type2].mySCPISM_T->C_i;      
      Real eta_j = topo->atoms[atom2].mySCPISM_A->eta;
      topo->atoms[atom2].mySCPISM_A->bornRadius += eta_j * f_ij * exp(-C_j * dist);
      //**********************************

      // If atom 2 is a polar H+, accumulate polar
      // fraction and derivative
      if (topo->atomTypes[type2].mySCPISM_T->isHbonded == PH &&
          topo->atomTypes[type1].mySCPISM_T->isHbonded == PA &&
          (topo->atoms[atom2].residue_seq !=
           topo->atoms[atom1].residue_seq)) {
        Real E_i = 0.80; // Currently all polar H+ have this value
        Real g_i;
        if (topo->atoms[atom2].name == "HN" && topo->atoms[atom1].name == "O")
          g_i = -0.378;
        else g_i = topo->atomTypes[type2].mySCPISM_T->g_i;
        Real g_j = topo->atomTypes[type1].mySCPISM_T->g_i;
        //Eq. (2)
        topo->atoms[atom2].mySCPISM_A->bornRadius += g_i * g_j * f_ij * exp(-E_i * dist);

      }
    }

    static void accumulateEnergy(ScalarStructure *energies, Real energy) {
      // Nothing to accumulate here
    }

    static Real getEnergy(const ScalarStructure *energies) {return 0;}

    //reset SCPISM struct before calculating
    virtual void preProcess(const GenericTopology *apptopo, const Vector3DBlock *positions) {
      const unsigned int atoms = apptopo->atoms.size();
      for(unsigned int i=0;i<atoms;i++)
        apptopo->atoms[i].mySCPISM_A->preForce();
	  }
    
    static void postProcess(const GenericTopology *topo, ScalarStructure *energies, Vector3DBlock *forces) {
      
    }
    
    static void parallelPostProcess(const GenericTopology *topo, ScalarStructure *energies) {
      
      const unsigned int atomnumber = topo->atoms.size();
      
      // Copy Radii
      Real *radii = new Real[ atomnumber ];
      
      //put radii (minus zeta) into array
      for( unsigned int i = 0; i < atomnumber; i++ ){
        radii[i] = topo->atoms[i].mySCPISM_A->bornRadius - topo->atoms[i].mySCPISM_A->zeta;
      }
      
      //sum accross nodes
      Parallel::reduce(radii, radii + atomnumber); //reduceSlaves only?
      
      //put radii back and add in zeta
      for( unsigned int i = 0; i < atomnumber; i++ ){
        topo->atoms[i].mySCPISM_A->bornRadius = radii[i] + topo->atoms[i].mySCPISM_A->zeta;
      }
      
      delete [] radii;
      
    }
    
    //do parallel post process?
    static bool doParallelPostProcess() {
      return true;
    }

    // Parsing
    static std::string getId() {return keyword;}
    void getParameters(std::vector<Parameter> &) const;
    static unsigned int getParameterSize() {return 1;}

    static BornRadii make(const std::vector<Value> &);

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // My data members
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    static const std::string keyword;

  private:
    int myBornSwitch;
  };

}
#endif /* BORNRADII_H */
