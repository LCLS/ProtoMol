/*  -*- c++ -*-  */
#ifndef NORMALMODEQUADRATIC_H
#define NORMALMODEQUADRATIC_H

#include <protomol/integrator/STSIntegrator.h>
#include <protomol/integrator/normal/NormalModeUtilities.h>
#include <string>
#include <vector>

#include "protomol/type/Real.h"

namespace ProtoMol {
class ProtoMolApp;
class Value;
class Vector3DBlock;
struct Parameter;
}  // namespace ProtoMol

namespace ProtoMol
{

  class ForceGroup;
  class ScalarStructure;

  //__________________________________________________ NormalModeQuadratic
  class NormalModeQuadratic : public STSIntegrator, public NormalModeUtilities
  {
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // Constructors, destructors, assignment
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    public:
      NormalModeQuadratic();
      NormalModeQuadratic( Real timestep, int firstmode, int nummode,
                           Real temperature, int cyclest, bool scle, bool stepm,
                           ForceGroup *overloadedForces );

      ~NormalModeQuadratic();

      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // New methods of class NormalModeQuadratic
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    private:
      Vector3DBlock* subspaceProj( double *tmpC, Vector3DBlock * iPos );

      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // From class Makeable
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    public:
      virtual std::string getIdNoAlias() const {return keyword;}
      virtual unsigned int getParameterSize() const {return 7;}
      virtual void getParameters( std::vector<Parameter>& parameters ) const;

      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // From class Integrator
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    public:
      virtual void initialize( ProtoMolApp* app );
      virtual long run( const long numTimesteps );

      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // From class STSIntegrator
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    private:
      virtual STSIntegrator* doMake( const std::vector<Value>& values, ForceGroup* fg ) const;

      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // My data members
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    public:
      static const std::string keyword;

    private:
      NormalModeUtilities *myPreviousNormalMode;
      Vector3DBlock* ex0;
      double *cPos;
      Real total_time;
      int cycleSteps, currMode;
      bool fScale, stepModes;
  };
}

#endif


