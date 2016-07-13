#include <protomol/ProtoMolApp.h>
#include <protomol/base/StringUtilities.h>
#include <protomol/base/SystemUtilities.h>
#include <protomol/io/CheckpointConfigWriter.h>
#include <protomol/io/XYZWriter.h>
#include <protomol/module/MainModule.h>
#include <protomol/output/OutputCheckpoint.h>
#include <iostream>

#include "protomol/base/Exception.h"
#include "protomol/base/Random.h"
#include "protomol/config/Configuration.h"
#include "protomol/config/ConstraintValueType.h"
#include "protomol/config/Parameter.h"
#include "protomol/config/Value.h"
#include "protomol/output/Output.h"
#include "protomol/output/OutputCache.h"
#include "protomol/topology/GenericTopology.h"
#include "protomol/type/SimpleTypes.h"

namespace ProtoMol {
class Vector3DBlock;
}  // namespace ProtoMol

#ifdef HAVE_LIBFAH
#include <cbang/os/File.h>

typedef cb::File fileStream;
#else
typedef std::fstream fileStream;
#endif

using namespace std;
using namespace ProtoMol::Report;
using namespace ProtoMol;

const string OutputCheckpoint::keyword("Checkpoint");

OutputCheckpoint::OutputCheckpoint(const string &name, int freq, int start, const string &posbase, const string &velbase)
	: Output(freq), current(start), name(name), posBase(posbase), velBase(velbase) {}

void OutputCheckpoint::doInitialize() {
	posBase = app->config["CheckpointPosBase"].get<std::string>();
	velBase = app->config["CheckpointVelBase"].get<std::string>();
}

void OutputCheckpoint::doIt(long step) {
	cout << "Checkpointing: Step " << step << ". . ." << flush;

	WritePositions(step);
	WriteVelocities(step);
	WriteConfig(step);

	//  Remove old checkpoint fil
	SystemUtilities::unlink(Append(Append(posBase, current - 1), ".pos"));
	SystemUtilities::unlink(Append(Append(velBase, current - 1), ".vel"));

	current += 1;

	cout << "done" << endl;
}

void OutputCheckpoint::doRun(long step) {
	const long firstStep = toLong(app->config["firststep"]);
	const long finalStep = firstStep + toLong(app->config["numsteps"]);

	if( step != firstStep && step != finalStep ) {
		if( getOutputFreq() > 0 && ( step % getOutputFreq()) == 0 ) {
			doIt(step);
		}
	}
}

Output *OutputCheckpoint::doMake(const vector<Value> &values) const {
	return new OutputCheckpoint(values[0], toInt(values[1]), values[2], values[3], values[4]);
}

bool OutputCheckpoint::isIdDefined(const Configuration *config) const {
	return config->valid(getId());
}

void OutputCheckpoint::getParameters(vector<Parameter> &parameter) const {
	parameter.push_back(Parameter(getId(), Value(name, ConstraintValueType::NotEmpty())));
	parameter.push_back(Parameter(getId() + "Freq", Value(outputFreq, ConstraintValueType::Positive()), Text("output frequency")));
	parameter.push_back(Parameter(getId() + "Start", Value(current, ConstraintValueType::NotNegative())));
	parameter.push_back(Parameter(getId() + "PosBase", Value(posBase, ConstraintValueType::NoConstraints())));
	parameter.push_back(Parameter(getId() + "VelBase", Value(velBase, ConstraintValueType::NoConstraints())));
}

bool OutputCheckpoint::adjustWithDefaultParameters(vector<Value> &values, const Configuration *config) const {
	if( !checkParameterTypes(values)) { return false; }

	if( config->valid(InputOutputfreq::keyword) && !values[1].valid()) {
		values[1] = ( *config )[InputOutputfreq::keyword];
	}

	if( !values[0].valid()) { values[0] = name; }
	if( !values[2].valid()) { values[2] = 0; }

	if( !values[3].valid()) {
		if( config->valid("gromacstprfile") ) {
			const std::string temp = (*config)["gromacstprfile"];
			values[3] = temp.substr(0, temp.rfind('.') + 1);
		} else {
			const std::string temp = (*config)["posfile"];
			values[3] = temp.substr(0, temp.rfind('.') + 1);
		}
	}

	if( !values[4].valid()) {
		if( config->valid("gromacstprfile") ) {
			const std::string temp = (*config)["gromacstprfile"];
			values[4] = temp.substr(0, temp.rfind('.') + 1);
		} else {
			if( !config->valid("velfile")) {
				values[4] = values[3];
			} else {
				const std::string temp = (*config)["velfile"];
				values[4] = temp.substr(0, temp.rfind('.') + 1);
			}
		}
	}

	return checkParameters(values);
}

void OutputCheckpoint::WritePositions(long step) {
	string posFile = Append(Append(posBase, current), ".pos");

	XYZWriter posWriter;
	if( !posWriter.open(posFile)) {
		THROWS("Can't open " << getId() << " '" << posFile + "'.");
	}

	const Vector3DBlock *pos = &app->positions;
	posWriter.setComment("Time : " + toString(app->outputCache.getTime()) +
						 ", step : " + toString(step) + ".");

	if( !posWriter.write(*pos, app->topology->atoms, app->topology->atomTypes)) {
		THROWS("Could not write " << getId() << " '" << posFile << "'.");
	}
}

void OutputCheckpoint::WriteVelocities(long step) {
	string velFile = Append(Append(velBase, current), ".vel");

	XYZWriter velWriter;
	if( !velWriter.open(velFile)) {
		THROWS("Can't open " << getId() << " '" << velFile << "'.");
	}

	velWriter.setComment("Time : " + toString(app->outputCache.getTime()) +
						 ", step : " + toString(step) + ".");

	if( !velWriter.write(*&app->velocities, app->topology->atoms,
						 app->topology->atomTypes)) {
		THROWS("Could not write " << getId() << " '" << velFile << "'.");
	}
}

void OutputCheckpoint::WriteConfig(long step) {
	string confFile = name + ".tmp";

	CheckpointConfigWriter confWriter;
	if( !confWriter.open(confFile)) {
		THROWS("Can't open " << getId() << " '" << confFile << "'.");
	}

	if( !confWriter.write(current, step, Random::Instance(), app->integrator)) {
		THROWS("Could not write " << getId() << " '" << confFile << "'.");
	}

	confWriter.close();

	SystemUtilities::rename(confFile, name);
}
