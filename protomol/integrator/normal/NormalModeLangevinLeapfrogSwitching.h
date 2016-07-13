/*  -*- c++ -*-  */
#ifndef NORMALMODELANGLFSWITCHING_H
#define NORMALMODELANGLFSWITCHING_H

#include <protomol/integrator/MTSIntegrator.h>
#include <protomol/integrator/normal/NormalModeUtilities.h>
#include <protomol/type/Vector3DBlock.h>
#include <string>
#include <vector>

#include "protomol/type/Real.h"
#include "protomol/type/Vector3D.h"

namespace ProtoMol {

  class ForceGroup;
  class ScalarStructure;
class ProtoMolApp;
class StandardIntegrator;
class Value;
struct Parameter;

  //__________________________________________________ NormalModeLangevinLeapfrogSwitchingSwitching
  class NormalModeLangevinLeapfrogSwitching : public MTSIntegrator, public NormalModeUtilities {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    NormalModeLangevinLeapfrogSwitching();
    NormalModeLangevinLeapfrogSwitching(int cycles, int firstmode, int nummode, Real gammaInside, Real gammaOutside, int seed, Real temperature, bool gencn,
					Real switchon, Real cutoff, 
                            ForceGroup *overloadedForces, StandardIntegrator *nextIntegrator);
    ~NormalModeLangevinLeapfrogSwitching(); 

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class NormalModeLangevinLeapfrogSwitching
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  protected:
    void drift();  
    void doDrift();  
    void doHalfKick();

  public:

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Makeable
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual std::string getIdNoAlias() const{return keyword;}
    virtual unsigned int getParameterSize() const{return 10;}
    virtual void getParameters(std::vector<Parameter>& parameters) const;

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Integrator
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual void initialize(ProtoMolApp* app);
    virtual long run(const long numTimesteps);
  protected:
    virtual void addModifierAfterInitialize();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class STSIntegrator
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    virtual MTSIntegrator* doMake(const std::vector<Value>& values, ForceGroup* fg, StandardIntegrator *nextIntegrator)const;
  public:

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    static const std::string keyword;

  protected:
    Real myGammaInside;
    Real myGammaOutside;
    Vector3D myCenterOfMass;
    Real mySwitchOn;
    Real mySwitchOn2;
    Real myCutoff;
    Real myCutoff2;
    Real mySwitch1;
    Real mySwitch2;
    Real mySwitch3;

  private:
    bool genCompNoise;

  };

}

#endif


