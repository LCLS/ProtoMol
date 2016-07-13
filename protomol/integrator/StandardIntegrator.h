/*  -*- c++ -*-  */
#ifndef STANDARDINTEGRATOR_H
#define STANDARDINTEGRATOR_H

#include <protomol/integrator/Integrator.h>

namespace ProtoMol {
  class ForceGroup;
  //________________________________________ StandardIntegrator
  class GenericTopology;
  class ScalarStructure;
  class Vector3DBlock;

  class StandardIntegrator : public Integrator {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    friend class MTSIntegrator;
    friend class S2HMCIntegrator;
  public:
    StandardIntegrator();
    StandardIntegrator(ForceGroup *forceGroup);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class StandardIntegrator
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    virtual void doHalfKick();
    virtual void doKick();
    virtual void initializeForces();
    virtual void doDriftOrNextIntegrator() = 0;
    virtual void calculateForces();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Integrator
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual long run(const long numTimesteps);
    virtual Integrator *previous();
    virtual const Integrator *previous() const;

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  protected:
    StandardIntegrator *myPreviousIntegrator;
  };
  //________________________________________ INLINES
}
#endif
