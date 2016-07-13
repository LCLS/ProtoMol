#include <protomol/ProtoMolApp.h>
#include <protomol/analysis/AnalysisDihedral.h>
#include <protomol/factory/AnalysisFactory.h>
#include <protomol/module/AnalysisModule.h>

using namespace ProtoMol;

defineInputValue(InputAnalysis,"analysis");

void AnalysisModule::init(ProtoMolApp *app) {
	InputAnalysis::registerConfiguration(&app->config, true);

	app->analysisFactory.registerExemplar(new AnalysisDihedral());
}
