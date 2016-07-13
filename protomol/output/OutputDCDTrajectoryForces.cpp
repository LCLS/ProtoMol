#include <protomol/ProtoMolApp.h>
#include <protomol/base/Exception.h>
#include <protomol/config/Configuration.h>
#include <protomol/io/DCDTrajectoryWriter.h>
#include <protomol/module/MainModule.h>
#include <protomol/output/OutputDCDTrajectoryForces.h>
#include <ostream>

#include "protomol/config/ConstraintValueType.h"
#include "protomol/config/Parameter.h"
#include "protomol/config/Value.h"
#include "protomol/integrator/Integrator.h"
#include "protomol/output/Output.h"
#include "protomol/type/SimpleTypes.h"

using namespace std;
using namespace ProtoMol::Report;
using namespace ProtoMol;


const string OutputDCDTrajectoryForces::keyword("DCDForcesFile");


OutputDCDTrajectoryForces::OutputDCDTrajectoryForces() :
  dCD(0), minimalImage(false) {}


OutputDCDTrajectoryForces::OutputDCDTrajectoryForces(const string &filename,
                                               int freq, bool minimal) :
  Output(freq), dCD(new DCDTrajectoryWriter(filename)),
  minimalImage(minimal) {}


OutputDCDTrajectoryForces::~OutputDCDTrajectoryForces() {
  if (dCD) delete dCD;
}


void OutputDCDTrajectoryForces::doInitialize() {
  if (!dCD || !dCD->open())
    THROWS("Can not open '" << (dCD ? dCD->getFilename() : "")
           << "' for " << getId() << ".");
}

void OutputDCDTrajectoryForces::doRun(long) {
  if (!dCD->write(*(app->integrator->getForces())))
    THROWS("Could not write " << getId() << " '" << dCD->getFilename()
           << "'.");
}


void OutputDCDTrajectoryForces::doFinalize(long) {
  dCD->close();
}


Output *OutputDCDTrajectoryForces::doMake(const vector<Value> &values) const {
  return new OutputDCDTrajectoryForces(values[0], values[1], values[2]);
}


void OutputDCDTrajectoryForces::getParameters(vector<Parameter> &parameter)
const {
  parameter.push_back
    (Parameter(getId(), Value(dCD ? dCD->getFilename() : "",
                              ConstraintValueType::NotEmpty())));
  Output::getParameters(parameter);
  parameter.push_back
    (Parameter(keyword + "MinimalImage", Value(minimalImage),
               Text("whether the coordinates should be transformed to minimal "
                    "image or not (NA for Vel)")));
}

bool OutputDCDTrajectoryForces::adjustWithDefaultParameters(
  vector<Value> &values, const Configuration *config) const {
  if (!checkParameterTypes(values)) return false;

  if (config->valid(InputOutputfreq::keyword) && !values[1].valid())
    values[1] = (*config)[InputOutputfreq::keyword];

  if (config->valid(InputMinimalImage::keyword) && !values[2].valid())
    values[2] = (*config)[InputMinimalImage::keyword];

  return checkParameters(values);
}
