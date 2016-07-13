/*  -*- c++ -*-  */
#ifndef HESSIANINT_H
#define HESSIANINT_H

#include <protomol/force/Force.h>
#include <protomol/integrator/STSIntegrator.h>
#include <protomol/integrator/hessian/BlockHessian.h>
#include <protomol/integrator/hessian/BlockHessianDiagonalize.h>
#include <protomol/type/TypeSelection.h>
#include <protomol/type/Vector3DBlock.h>
#include <string>
#include <vector>

#include "protomol/type/Real.h"

namespace ProtoMol {
class ProtoMolApp;
class Value;
struct Parameter;
}  // namespace ProtoMol

using namespace std;

namespace ProtoMol {
  class ForceGroup;
  class ReducedHessAngle;
  class ScalarStructure;

  //____ HessianInt
  class HessianInt : public STSIntegrator {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    HessianInt();
    HessianInt(Real timestep, std::string evec_s, std::string eval_s,
               std::string hess_s, bool sorta, int fm, bool tef, 
               bool fdi, Real evt, int bvc, int rpb, Real bcd, bool masswt,
               bool bnm, bool aparm, bool geo, bool num, Real eps,
               ForceGroup *overloadedForces);
    ~HessianInt();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class HessianInt
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    void outputDiagHess(int numModes);
    Real calcQ();
    void numericalHessian();

  protected:
    void doKickdoDrift();
    void doHalfKickdoDrift();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Makeable
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual std::string getIdNoAlias() const {return keyword;}
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

  protected:

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    static const std::string keyword;
  private:
    typedef TypeSelection::Int<4>::type int32;
    double *eigVec;
    int totStep;
    unsigned int sz;
    std::string evecfile, evalfile, hessfile;
    bool sortOnAbs;
    unsigned int numberOfModes;
    bool textEig, fullDiag, massWeight, noseMass, autoParmeters;
    bool geometricfdof, numerichessians;
    BlockHessian hsn;     
    BlockHessianDiagonalize blockDiag;
    Real eigenValueThresh, blockCutoffDistance;
    int blockVectorCols, residuesPerBlock, residues_total_eigs;
    Real max_eigenvalue;
    Real epsilon;

  };
}
#endif

