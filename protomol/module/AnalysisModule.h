#ifndef ANALYSIS_MODULE_H
#define ANALYSIS_MODULE_H

#include <protomol/base/Module.h>
#include <protomol/config/InputValue.h>
#include <string>

#include "protomol/config/InputValueMacros.h"

namespace ProtoMol {
	class ProtoMolApp;

	declareInputValue(InputAnalysis, BOOL, NOCONSTRAINTS);

	class AnalysisModule : public Module {
		public:
			const std::string getName() const { return "Analysis"; }
			int getPriority() const { return 0;  }
			const std::string getHelp() const { return ""; }
			void init (ProtoMolApp *app);
	};
}

#endif	// OUTPUT_MODULE_H
