#include <protomol/ProtoMolApp.h>
#include <protomol/factory/TopologyFactory.h>
#include <protomol/module/TopologyModule.h>
#include <protomol/topology/CubicCellManager.h>
#include <protomol/topology/PeriodicBoundaryConditions.h>
#include <protomol/topology/Topology.h>
#include <protomol/topology/VacuumBoundaryConditions.h>
#include <vector>

#include "protomol/config/Configuration.h"
#include "protomol/config/Value.h"
#include "protomol/topology/GenericTopology.h"
#include "protomol/type/Vector.h"

using namespace std;
using namespace ProtoMol;

defineInputValue(InputBoundaryConditions, "boundaryConditions")
defineInputValue(InputCellManager, "cellManager")

void TopologyModule::init(ProtoMolApp *app) {
  Configuration *config = &app->config;
  GenericTopology *topo;
  TopologyFactory *factory = &app->topologyFactory;

  // Vacuum or normal boundary conditions
  topo = new Topology<VacuumBoundaryConditions, CubicCellManager>();
  factory->registerExemplar(topo, Vector<string>("NormalCubic"));

  // Periodic boundary conditions
  topo = new Topology<PeriodicBoundaryConditions, CubicCellManager>();
  factory->registerExemplar(topo);

  // Register input values
  InputBoundaryConditions::registerConfiguration(config);
  InputCellManager::registerConfiguration(config);
}

void TopologyModule::configure(ProtoMolApp *app) {
  Configuration &config = app->config;

  // Fix for old topology definition
  if (!config[GenericTopology::keyword].valid())
    config[GenericTopology::keyword] =
      config[InputBoundaryConditions::keyword].getString() +
      config[InputCellManager::keyword].getString();
}
