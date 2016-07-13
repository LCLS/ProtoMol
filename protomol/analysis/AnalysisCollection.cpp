#include <protomol/ProtoMolApp.h>
#include <protomol/analysis/Analysis.h>
#include <protomol/analysis/AnalysisCollection.h>
#include <protomol/base/Exception.h>
#include <string>

#include "protomol/output/OutputCache.h"

using namespace ProtoMol;

AnalysisCollection::~AnalysisCollection() {
	for( iterator i = begin(); i != end(); i++ ) {
		delete ( *i );
	}
}

void AnalysisCollection::initialize(const ProtoMolApp *app) {
	this->app = app;
	for( iterator i = begin(); i != end(); i++ ) {
		( *i )->initialize(app);
	}
}

bool AnalysisCollection::run(long step) {
	bool outputRan = false;

	app->outputCache.uncache();
	for( iterator i = begin(); i != end(); ++i ) {
		if( ( *i )->isOnOutput() ) {
			outputRan |= ( *i )->run(step);
		}
	}

	return outputRan;
}

bool AnalysisCollection::shouldStop() {
	for( iterator i = begin(); i != end(); ++i ) {
		if(( *i )->shouldStop() ) {
			return true;
		}
	}
	return false;
}

void AnalysisCollection::finalize(long step) {
	app->outputCache.uncache();
	for( iterator i = begin(); i != end(); i++ ) {
		( *i )->finalize(step);
	}
}

void AnalysisCollection::adoptAnalysis(Analysis *output) {
	if( !output ) { THROW("null pointer"); }
	outputList.push_back(output);
}
