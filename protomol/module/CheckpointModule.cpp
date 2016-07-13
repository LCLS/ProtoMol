#include "CheckpointModule.h"

#include <protomol/ProtoMolApp.h>
#include <protomol/factory/OutputFactory.h>
#include <protomol/io/CheckpointConfigReader.h>
#include <protomol/output/OutputCheckpoint.h>
#include <iostream>

#include "protomol/base/Random.h"
#include "protomol/config/Configuration.h"
#include "protomol/config/Value.h"

using namespace std;
using namespace ProtoMol;

void CheckpointModule::init(ProtoMolApp *app) {
	OutputFactory &f = app->outputFactory;
	f.registerExemplar(new OutputCheckpoint());
}

void CheckpointModule::configure(ProtoMolApp *app) {
}

void CheckpointModule::read(ProtoMolApp *app) {
	if( !enabled ) { return; }

	Configuration &config = app->config;

	std::cout << "Reading Checkpoint Base Data" << std::endl;

	CheckpointConfigReader confReader;
	if( confReader.open(config["Checkpoint"], ios::in)) {
		confReader.readBase(config, Random::Instance());
	}
}

void CheckpointModule::postBuild(ProtoMolApp *app) {
	if( !enabled ) { return; }

	Configuration &config = app->config;

	std::cout << "Reading Checkpoint Integrator Data" << std::endl;
	// Load integrator data
	CheckpointConfigReader confReader;
	if( confReader.open(config["Checkpoint"], ios::in)) {
		confReader.readIntegrator(app->integrator);
	}
}

string CheckpointModule::WithoutExt(const string &path) {
	return path.substr(0, path.rfind(".") + 1);
}
