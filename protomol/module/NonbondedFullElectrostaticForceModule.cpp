#include <math.h>
#include <protomol/ProtoMolApp.h>
#include <protomol/base/StringUtilities.h>
#include <protomol/force/nonbonded/NonbondedFullEwaldSystemForce.h>
#include <protomol/module/NonbondedFullElectrostaticForceModule.h>
#include <protomol/module/TopologyModule.h>
#include <protomol/switch/CutoffSwitchingFunction.h>
#include <protomol/topology/CubicCellManager.h>
#include <protomol/topology/PeriodicBoundaryConditions.h>

#include "protomol/config/Configuration.h"
#include "protomol/config/Value.h"
#include "protomol/factory/ForceFactory.h"
#include "protomol/force/Force.h"
#include "protomol/topology/GenericTopology.h"
#include "protomol/type/Vector.h"
#include "protomol/type/Vector3D.h"

using namespace std;
using namespace ProtoMol;

void NonbondedFullElectrostaticForceModule::registerForces(ProtoMolApp *app) {
  ForceFactory &f = app->forceFactory;
  string boundConds = app->config[InputBoundaryConditions::keyword];

  typedef CubicCellManager CCM;
  typedef PeriodicBoundaryConditions PBC;
  typedef CutoffSwitchingFunction Cutoff;
#define FullEwald NonbondedFullEwaldSystemForce

  if (equalNocase(boundConds, PBC::keyword)) {
    // Full Ewald
    f.reg(new FullEwald<PBC,CCM,true,true,true,Cutoff>(),Vector<std::string>("CoulombEwald"));
    f.reg(new FullEwald<PBC,CCM,true,false,false,Cutoff>());
    f.reg(new FullEwald<PBC,CCM,true,false,true,Cutoff>());
    f.reg(new FullEwald<PBC,CCM,false,true,false,Cutoff>());
    f.reg(new FullEwald<PBC,CCM,false,true,true,Cutoff>());
    f.reg(new FullEwald<PBC,CCM,false,false,true,Cutoff>());

  } 

}
