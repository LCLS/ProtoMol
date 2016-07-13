/* -*- c++ -*- */
#ifndef HARMONICRESTRAINTFORCE_H
#define HARMONICRESTRAINTFORCE_H

#include <protomol/base/Report.h>
#include <protomol/config/Parameter.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/topology/Topology.h>
#include <protomol/topology/TopologyUtilities.h>
#include <protomol/type/ScalarStructure.h>
#include <string>
#include <vector>

#include "protomol/topology/ExclusionTable.h"
#include "protomol/type/Real.h"
#include "protomol/type/Vector3D.h"
#include "protomol/type/Vector3DBlock.h"

namespace ProtoMol {
class GenericTopology;
class Value;
struct Parameter;
}  // namespace ProtoMol

using namespace ProtoMol::Report;

namespace ProtoMol {
  //____ HarmonicRestraintForce
  class HarmonicRestraintForce {
  public:
    enum {DIST_R2 = 1};
    enum {CUTOFF = 0};
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class LennardJonesForce
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    HarmonicRestraintForce();
    HarmonicRestraintForce(Real k, Real r, Real so, Real co);
    
    void operator()(Real &energy, Real &force,
                    Real /*distSquared*/, Real rDistSquared, const Vector3D &,
                    const GenericTopology *topo, int atom1, int atom2,
                    ExclusionClass excl)  {
      
      //test if the diagonal set (as unique)
      if( atom2 == atom1 + 1 ){
        energy = doAtomCalculation(atom1);
        
        const int size = pos->size();
        
        //last atom?, runs from i=1-N->1, j=i+1->N
        if( atom1 == size - 2 ){  //now don't need to test atom2
          energy += doAtomCalculation(atom2);
        }
      }

    }
    
    //actual calculation
    Real doAtomCalculation( int atom1 ){
              
      Vector3D diff = (*pos)[atom1] - centerofmass;

      Real distance = diff.norm();
      Real distSquared = distance * distance;
      
      //switch
      Real value;
      Real deriv = 0.0;
      if (distSquared > myCutoff2)
        value = 0.0;
      else if (distSquared >= mySwitchon2) {
        Real c2 = myCutoff2 - distSquared;
        Real c4 = c2 * (mySwitch2 + 2.0 * distSquared);
        value = mySwitch1 * (c2 * c4);
        deriv = mySwitch3 * (c2 * c2 - c4);
      } else
        value = 1.0;
      
      //compensate for compliment
      value = 1.0 - value;
      deriv *= -1.0;
      
      //calc energy
      Real energy = sphereK * ( distance - sphereradius ) * ( distance - sphereradius );
      
      //Real myForce = 
      const Real dpotdr = 2.0 * sphereK * (distance - sphereradius);  // Calculate dpot/dr
      
      //chain rule
      const Real combinedForce = dpotdr * value + energy * deriv;
      
      //fix energy for switch
      energy *= value;
      
      const Real scalarForce = -combinedForce / distance;
      
      // Calculate force on atom1 due to atom2.
      Vector3D force1(diff * scalarForce);//diff * (-dpotdr / distance));
      
      // Add to the total force.
      myForces[atom1] = force1;
      
      return energy;
      
    }

    static void accumulateEnergy(ScalarStructure *energies, Real energy) {
      (*energies)[ScalarStructure::OTHER] += energy;
    }

    static Real getEnergy(const ScalarStructure *energies) {
      return (*energies)[ScalarStructure::OTHER];
    }
    
    void preProcess(const GenericTopology *apptopo, const Vector3DBlock *positions) {
      const unsigned int size = positions->size();
      
      if( firstEvaluate ){
        
        firstEvaluate = false;
        
        
        //resize my forces
        myForces.resize(size);
        
        //get center of mass
        centerofmass = centerOfMass(positions, apptopo);
        
      }
      
      //reset myForces
      for(unsigned int i=0; i<size; i++){
        myForces[i] = Vector3D(); //initializes to zero
      }
      
      //save position pointer
      pos = positions;
		  
    }
    
    //add in forces
    void postProcess(const GenericTopology *topo, ScalarStructure *energies, Vector3DBlock *forces){
      const unsigned int size = forces->size();

      for( unsigned int i=0; i<size; i++ ){
        (*forces)[i] += myForces[i];
      }

    }

    static void parallelPostProcess(const GenericTopology *topo, ScalarStructure *energies) {
      
    }

    //do parallel post process?
    static bool doParallelPostProcess() {
      return false;
    }

    // Parsing
    static std::string getId() {return keyword;}
    static unsigned int getParameterSize();
    void getParameters(std::vector<Parameter> &) const;

    static HarmonicRestraintForce make(const std::vector<Value> &);
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    static const std::string keyword;
  private:
    Vector3DBlock myForces;
    Vector3D centerofmass;
    const Vector3DBlock *pos;
    Real sphereK, sphereradius;
    
    //"first" flag
    bool firstEvaluate;
    
    //switch
    Real mySwitchon, mySwitchon2, myCutoff, myCutoff2, mySwitch1, mySwitch2,
    mySwitch3;
  
  };
  //____ INLINES
}

#endif /* HARMONICRESTRAINTFORCE_H */
