#ifndef NORMALMODEOPENMM_H
#define NORMALMODEOPENMM_H

#include <protomol/integrator/normal/NormalModeUtilities.h>
#include <protomol/integrator/openMM/OpenMMIntegrator.h>
#include <string>
#include <vector>

#include "protomol/type/Real.h"

namespace ProtoMol {
	class ForceGroup;
	class ScalarStructure;
class ProtoMolApp;
class STSIntegrator;
class Value;
struct Parameter;

	class NormalModeOpenMM : public OpenMMIntegrator, public NormalModeUtilities {
		public:
			NormalModeOpenMM();
			NormalModeOpenMM( const std::vector<Value>& base, const std::vector<Value>& params, ForceGroup *forces );
			~NormalModeOpenMM();

			virtual std::string getIdNoAlias() const { return keyword; }
			virtual void getParameters( std::vector<Parameter>& parameters ) const;
			virtual unsigned int getParameterSize() const;

			virtual void initialize( ProtoMolApp *app );
			virtual long run( const long numTimesteps );
		private:
			virtual STSIntegrator *doMake( const std::vector<Value>& values, ForceGroup *fg )const;
		public:
			static const std::string keyword;
		private:
			bool shoudForceRediagOnMinFail;
			bool mRediagOnQuadratic;
			bool mAlwaysQuadratic;
			bool mRediagQuadraticLambda;
			bool mMetropolisMinimize;

			int mBlockSplit;

			int mRediagonalizationFrequency;

			Real mMinimizationLimit;
			int mResiduesPerBlock, mBlockDOF;
			Real mBlockDelta, mSDelta;
			int mModes, mBlockPlatform;
			bool mProtomolDiagonalize;
			NormalModeUtilities *myPreviousNormalMode;
	};
}

#endif
