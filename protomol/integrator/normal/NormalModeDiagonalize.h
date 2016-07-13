/*  -*- c++ -*-  */
#ifndef NORMALMODEDIAGONALIZE_H
#define NORMALMODEDIAGONALIZE_H

#include <protomol/base/Timer.h>
#include <protomol/integrator/MTSIntegrator.h>
#include <protomol/integrator/hessian/BlockHessian.h>
#include <protomol/integrator/hessian/BlockHessianDiagonalize.h>
#include <protomol/integrator/normal/NormalModeUtilities.h>
#include <protomol/type/BlockMatrix.h>
#include <protomol/type/Vector3DBlock.h>
#include <iosfwd>
#include <string>
#include <vector>

#include "protomol/type/Real.h"

namespace ProtoMol {
class ProtoMolApp;
class StandardIntegrator;
class Value;
struct Parameter;
}  // namespace ProtoMol

namespace ProtoMol
{

  class ForceGroup;
  class ScalarStructure;

  //__________________________________________________ NormalModeDiagonalize
  class NormalModeDiagonalize :
        public MTSIntegrator, public NormalModeUtilities
  {

    private:
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // Types and Enums
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      enum {MAX_ATOMS_PER_RES = 30};
      enum {REGRESSION_T = 0};

      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // Constructors, destructors, assignment
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    public:
      NormalModeDiagonalize();
      NormalModeDiagonalize( int cycles, int redi, bool fDiag,
                             bool rRand,
                             Real redhy, Real eTh, int bvc, int rpb, Real dTh, 
                             bool apar, bool adts, bool pdm, Real ml, int maxit,
                             bool geo, bool num, ForceGroup *overloadedForces,
                             StandardIntegrator *nextIntegrator );
      ~NormalModeDiagonalize();

      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // New methods of class NormalModeDiagonalize
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    protected:
      void utilityCalculateForces();

    private:
    public:

      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // From class Makeable
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    public:
      virtual std::string getIdNoAlias() const {return keyword;}
      virtual unsigned int getParameterSize() const {return 16;}
      virtual void getParameters( std::vector<Parameter>& parameters ) const;

      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // From class Integrator
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    public:
      virtual void initialize( ProtoMolApp *app );
      virtual long run( const long numTimesteps );
    protected:
      //virtual void addModifierAfterInitialize();

      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // From class STSIntegrator
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    private:
      virtual MTSIntegrator* doMake( const std::vector<Value>& values,
                                     ForceGroup *fg,
                                     StandardIntegrator *nextIntegrator ) const;
    public:

      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // New methods of class NormalModeUtilities
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    protected:
      virtual void streamRead( std::istream& inStream );
      virtual void streamWrite( std::ostream& outStream ) const;
    private:
      Real doDiagonalization();

      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // My data members
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    public:
      static const std::string keyword;

    private:
      //Hessian/Diag Hessian
      BlockHessian rHsn;
      BlockHessianDiagonalize blockDiag;

      bool firstDiag, fullDiag, removeRand;
      int rediagCount, nextRediag;
      bool validMaxEigv;
      NormalModeUtilities *myNextNormalMode, *myLastNormalMode;
      Real rediagHysteresis;

      //Diagnostic data
      int hessianCounter, rediagCounter, rediagUpdateCounter;

      //Residues
      Real eigenValueThresh, blockCutoffDistance;
      int blockVectorCols, residuesPerBlock;

      //Diagnostics
      unsigned int memory_Hessian, memory_eigenvector;

      //Checkpointing
      bool checkpointUpdate;    
      double origCEigVal, origTimestep;

      //auto-parameters?
      bool autoParmeters;
    
      //adaptive timestep
      bool adaptiveTimestep;

      //post diagonalize minimize
      bool postDiagonalizeMinimize;
      Real minLim;
      int maxMinSteps;

      //numerical and geometric Hessian
      bool geometricfdof, numerichessians;

  };
}

#endif


