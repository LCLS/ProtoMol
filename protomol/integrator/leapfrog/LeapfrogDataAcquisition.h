/*  -*- c++ -*-  */
#ifndef LEAPFROGDATAACQUISITION_H
#define LEAPFROGDATAACQUISITION_H


#include <protomol/integrator/leapfrog/LeapfrogIntegrator.h>
#include <string>
#include <vector>

#include "protomol/type/Real.h"


namespace ProtoMol {
  class DCDTrajectoryWriter;
  class ForceGroup;
  class ScalarStructure;
class ProtoMolApp;
class STSIntegrator;
class Value;
class Vector3DBlock;
struct Parameter;

  class LeapfrogDataAcquisition : public LeapfrogIntegrator {

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    LeapfrogDataAcquisition();
    LeapfrogDataAcquisition(Real timestep, string dcdfile, ForceGroup *overloadedForces);
    ~LeapfrogDataAcquisition();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Makeable
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual std::string getIdNoAlias() const{return keyword;}
    virtual unsigned int getParameterSize() const{return 2;}
    virtual void getParameters(std::vector<Parameter> &parameters) const;

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Integrator
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual void initialize(ProtoMolApp *app);
    virtual long run(const long numTimesteps);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class STSIntegrator
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    virtual STSIntegrator *doMake(const std::vector<Value> &values,
                                  ForceGroup *fg) const;

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    static const std::string keyword;
    int getNumNonWaters() {return numNonWaters;}

  private:
    void writeDCD();
    DCDTrajectoryWriter *myWriterF, *myWriterX, *myWriterV;
    string myDCDFile;
    int numNonWaters;
    Vector3DBlock *nonWaterF, *nonWaterX, *nonWaterV;
  };

}


#endif
