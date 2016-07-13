#include <protomol/ProtoMolApp.h>
#include <protomol/base/Exception.h>
#include <protomol/integrator/Integrator.h>
#include <protomol/io/XYZTrajectoryWriter.h>
#include <protomol/output/OutputXYZTrajectoryForce.h>
#include <protomol/topology/GenericTopology.h>
#include <stddef.h>
#include <ostream>

#include "protomol/config/ConstraintValueType.h"
#include "protomol/config/Parameter.h"
#include "protomol/config/Value.h"
#include "protomol/output/Output.h"

using namespace std;
using namespace ProtoMol::Report;
using namespace ProtoMol;


const string OutputXYZTrajectoryForce::keyword("XYZForceFile");


OutputXYZTrajectoryForce::OutputXYZTrajectoryForce() : xYZ() {}


OutputXYZTrajectoryForce::OutputXYZTrajectoryForce(const string &filename,
                                                   int freq) :
  Output(freq), xYZ(new XYZTrajectoryWriter(filename)) {}


OutputXYZTrajectoryForce::~OutputXYZTrajectoryForce() {
  if (xYZ) delete xYZ;
}


void OutputXYZTrajectoryForce::doInitialize() {
  if (xYZ == NULL || !xYZ->open())
    THROWS("Can not open '" << (xYZ ? xYZ->getFilename() : "")
           << "' for " << getId() << ".");
}


void OutputXYZTrajectoryForce::doRun(long) {
  if (!xYZ->write(*(app->integrator->getForces()), app->topology->atoms,
                    app->topology->atomTypes))
    THROWS("Could not write " << getId() << " '" << xYZ->getFilename()
           << "'.");
}


void OutputXYZTrajectoryForce::doFinalize(long) {
  xYZ->close();
}


Output *OutputXYZTrajectoryForce::doMake(const vector<Value> &values) const {
  return new OutputXYZTrajectoryForce(values[0], values[1]);
}


void OutputXYZTrajectoryForce::getParameters(vector<Parameter> &parameter)
const {
  parameter.push_back
    (Parameter(getId(), Value(xYZ ? xYZ->getFilename() : "",
                              ConstraintValueType::NotEmpty())));
  Output::getParameters(parameter);
}
