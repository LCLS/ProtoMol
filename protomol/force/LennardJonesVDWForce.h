/* -*- c++ -*- */
#ifndef LENNARDJONESVDWFORCE_H
#define LENNARDJONESVDWFORCE_H

#include <protomol/base/Report.h>
#include <protomol/config/Parameter.h>
#include <protomol/topology/ExclusionTable.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/topology/LennardJonesParameters.h>
#include <protomol/type/ScalarStructure.h>
#include <stddef.h>
#include <string>
#include <vector>

#include "protomol/base/MathUtilities.h"
#include "protomol/config/ConstraintValueType.h"
#include "protomol/config/Value.h"
#include "protomol/topology/Atom.h"
#include "protomol/topology/LennardJonesParameterTable.h"
#include "protomol/type/Real.h"
#include "protomol/type/Vector3D.h"

namespace ProtoMol {
class Vector3DBlock;
}  // namespace ProtoMol

using namespace ProtoMol::Report;

namespace ProtoMol {
  //____ LennardJonesVDWForce
  class LennardJonesVDWForce {
  public:
    enum {DIST_R2 = 1};
    enum {CUTOFF = 0};

  public:
    virtual ~LennardJonesVDWForce() {} // Compiler needs this

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class LennardJonesVDWForce
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    LennardJonesVDWForce() : onvdwratio(1.0), offvdwratio(1.5) {
    }

    LennardJonesVDWForce(Real on, Real off) : onvdwratio(on), offvdwratio(off) {
      if(onvdwratio >= offvdwratio){
        report << error << "LennardJonesVDWForce(): Switch on less than switch off." << endr; 
      }
    }
    
    void operator()(Real &energy, Real &force,
                    Real /*distSquared*/, Real rDistSquared, const Vector3D &,
                    const GenericTopology *topo, int atom1, int atom2,
                    ExclusionClass excl) const {
#ifdef DEBUG_LJ_DISTANCE_CHECK
      if (rDistSquared < 1.0 / 0.25)
        report << warning << "LennardJonesVDWForce::operator(): atom " << atom1 <<
        " and " << atom2 << " get closer"
                << " than 0.5 angstroms (=" <<
        sqrt(1.0 / rDistSquared) << ")." << endr;
#endif
        
      //Calculate switch
      Real deriv = 0.0, value = 1.0;
      
      //check VDW data available
      if(topo->atoms[atom1].myGBSA_T != NULL && topo->atoms[atom2].myGBSA_T != NULL ){
        const Real distSquared = 1.0 / rDistSquared;
        
        //get combined atom radii
        const Real switchrange = topo->atoms[atom1].myGBSA_T->vanDerWaalRadius + topo->atoms[atom2].myGBSA_T->vanDerWaalRadius;
//       	if(sqrt(distSquared) < switchrange * 1.5) report << debug(1) << "LennardJonesVDWForce::operator(): sr = " << switchrange << ", dist = " << sqrt(distSquared) << endr;
        
        if(switchrange == 0){
          report << error << "LennardJonesVDWForce::operator(): Switch range zero." << endr;  
        }
      
        //get square of switch on and off
        const Real myCutoff2 = offvdwratio * switchrange * offvdwratio * switchrange, mySwitchon2 = onvdwratio * switchrange * onvdwratio * switchrange;
        
        //apply switch?
        if (distSquared > myCutoff2)
          value = 0.0;
        else if (distSquared >= mySwitchon2){          
          const Real mySwitch1 = 1.0 / power<3>(myCutoff2 - mySwitchon2),
          mySwitch2 = myCutoff2 - 3.0 * mySwitchon2,
          mySwitch3 = 4.0 / power<3>(myCutoff2 - mySwitchon2);
          
          Real c2 = myCutoff2 - distSquared;
          Real c4 = c2 * (mySwitch2 + 2.0 * distSquared);
          value = mySwitch1 * (c2 * c4);
          deriv = mySwitch3 * (c2 * c2 - c4);
        }
      }else{
        report << error << "LennardJonesVDWForce::operator(): No Van Der Waals data (myGBSA_T)." << endr;
      }
        
      //get params
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
      energy = (r12A - r6B) * value;
      force = (12.0 * r12A * rDistSquared - 6.0 * r6B * rDistSquared) * value - (r12A - r6B) * deriv;
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
    static unsigned int getParameterSize() {return 2;}
    void getParameters(std::vector<Parameter> &parameters) const {
      parameters.push_back
      (Parameter("-onvdwratio",
                 Value(onvdwratio, ConstraintValueType::Positive()),
                 0.8));
      parameters.push_back
      (Parameter("-offvdwratio",
                 Value(offvdwratio, ConstraintValueType::Positive()),
                 1.0));
    }

    static LennardJonesVDWForce make(const std::vector<Value> &values) {
      return LennardJonesVDWForce(values[0], values[1]);
    }

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    static const std::string keyword;
  private:
        Real onvdwratio, offvdwratio;
  };
  //____ INLINES
}

#endif /* LENNARDJONESVDWFORCE_H */
