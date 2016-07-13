#include <protomol/ProtoMolApp.h>
#include <protomol/base/ModuleManager.h>
#include <protomol/config/CommandLine.h>
#include <protomol/module/CommandLineModule.h>
#include <protomol/module/TopologyModule.h>
#include <protomol/topology/PeriodicBoundaryConditions.h>
#include <protomol/topology/VacuumBoundaryConditions.h>
#include <iostream>

#include "protomol/analysis/Analysis.h"
#include "protomol/base/Factory.h"
#include "protomol/base/StringUtilities.h"
#include "protomol/config/CommandLineOption.h"
#include "protomol/config/Configuration.h"
#include "protomol/config/Value.h"
#include "protomol/factory/AnalysisFactory.h"
#include "protomol/factory/ForceFactory.h"
#include "protomol/factory/IntegratorFactory.h"
#include "protomol/factory/OutputFactory.h"
#include "protomol/factory/TopologyFactory.h"
#include "protomol/force/Force.h"
#include "protomol/integrator/Integrator.h"
#include "protomol/output/Output.h"
#include "protomol/topology/GenericTopology.h"

using namespace std;
using namespace ProtoMol;

void CommandLineModule::init(ProtoMolApp *app) {
  this->app = app;
  CommandLineOption::ActionBase *action;
  CommandLine &cmdLine = app->cmdLine;

  action = new CommandLineOption::Action<CommandLine>
    (&cmdLine, &CommandLine::usageAction);
  cmdLine.add('h', "help", action, "Print this help screen.");

  action = new CommandLineOption::Action<CommandLine>
    (&cmdLine, &CommandLine::splashAction);
  cmdLine.add(0, "splash", action, "Print splash screen and exit.");

  action = new CommandLineOption::Action<CommandLineModule>
    (this, &CommandLineModule::listTopologies);
  cmdLine.add(0, "topologies", action, "List topologies and exit.");

  action = new CommandLineOption::Action<CommandLineModule>
    (this, &CommandLineModule::listIntegrators);
  cmdLine.add(0, "integrators", action, "List integrators and exit.");

  action = new CommandLineOption::Action<CommandLineModule>
    (this, &CommandLineModule::listForces);
  cmdLine.add(0, "forces", action, "List forces and exit.");

  action = new CommandLineOption::Action<CommandLineModule>
    (this, &CommandLineModule::listOutputs);
  cmdLine.add(0, "outputs", action, "List outputs and exit.");

  action = new CommandLineOption::Action<CommandLineModule>
    (this, &CommandLineModule::listAnalysis);
  cmdLine.add(0, "analysis", action, "List analysis and exit.");

#ifdef DEBUG
  action = new CommandLineOption::Action<CommandLine>
    (&cmdLine, &CommandLine::enableStackTraceAction);

  cmdLine.add('X', "enable-stack-trace", action,
    "Enable stack trace output on exceptions.");
#endif
}

int CommandLineModule::listTopologies(const vector<string> &args) {
  cout << headerRow("Topologies") << endl << app->topologyFactory << endl;
  return -1;
}

int CommandLineModule::listIntegrators(const vector<string> &args) {
  cout << headerRow("Integrators") << endl << app->integratorFactory << endl;
  return -1;
}

int CommandLineModule::listForces(const vector<string> &args) {
  app->forceFactory.unregisterAllExemplars();

  // PeriodicBoundaryConditions
  app->config[InputBoundaryConditions::keyword] =
    PeriodicBoundaryConditions::keyword;

  app->modManager->registerForces(app);

  cout << headerRow("Forces periodic boundary conditions") << endl
       << app->forceFactory << endl;

  app->forceFactory.unregisterAllExemplars();

  // VacuumBoundaryConditions
  app->config[InputBoundaryConditions::keyword] =
    VacuumBoundaryConditions::keyword;

  app->modManager->registerForces(app);

  cout << headerRow("Forces vacuum boundary conditions") << endl
       << app->forceFactory << endl;

  app->forceFactory.unregisterAllExemplars();

  return -1;
}

int CommandLineModule::listOutputs(const vector<string> &args) {
  cout << headerRow("Outputs") << endl << app->outputFactory << endl;
  return -1;
}

int CommandLineModule::listAnalysis(const vector<string> &args) {
  cout << headerRow("Analysis") << endl << app->analysisFactory << endl;
  return -1;
}
