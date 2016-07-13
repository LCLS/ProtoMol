#include <protomol/ProtoMolApp.h>
#include <protomol/base/StringUtilities.h>
//GB
#include <protomol/force/GB/GBBornRadii.h>
#include <protomol/force/GB/GBPartialSum.h>
#include <protomol/force/OneAtomPairNoExclusion.h>
//SCPISM
#include <protomol/force/born/BornRadii.h>
#include <protomol/force/born/BornSelfForce.h>
#include <protomol/force/nonbonded/NonbondedIntermittentFullSystemForce.h>
#include <protomol/module/NonbondedIntermittentFullForceModule.h>
#include <protomol/module/TopologyModule.h>
#include <protomol/switch/CutoffSwitchingFunction.h>
#include <protomol/switch/UniversalSwitchingFunction.h>
#include <protomol/topology/PeriodicBoundaryConditions.h>
#include <protomol/topology/VacuumBoundaryConditions.h>
#include <vector>

#include "protomol/config/Configuration.h"
#include "protomol/config/Value.h"
#include "protomol/factory/ForceFactory.h"
#include "protomol/force/Force.h"

using namespace std;
using namespace ProtoMol;

void NonbondedIntermittentFullForceModule::registerForces(ProtoMolApp *app) {
  ForceFactory &f = app->forceFactory;
  string boundConds = app->config[InputBoundaryConditions::keyword];

  typedef PeriodicBoundaryConditions PBC;
  typedef VacuumBoundaryConditions VBC;
  typedef CutoffSwitchingFunction Cutoff;
  typedef UniversalSwitchingFunction Universal;
#define Complement ComplementSwitchingFunction 
#define IntermittentFullSystem NonbondedIntermittentFullSystemForce

  if (equalNocase(boundConds, PeriodicBoundaryConditions::keyword)) {


  } else if (equalNocase(boundConds,  VacuumBoundaryConditions::keyword)) {
    
    // SCPISM
    // Born radius
    f.reg(new IntermittentFullSystem<OneAtomPairNoExclusion<VBC, Cutoff, BornRadii> >());
    f.reg(new IntermittentFullSystem<OneAtomPairNoExclusion<VBC, Cutoff, BornSelfForce> >());
    
    // GB
    f.reg(new IntermittentFullSystem<OneAtomPairNoExclusion<VBC, Universal, GBBornRadii> >());
    f.reg(new IntermittentFullSystem<OneAtomPairNoExclusion<VBC, Universal, GBPartialSum> >());
    
  }
}
