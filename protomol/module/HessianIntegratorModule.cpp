#include <protomol/ProtoMolApp.h>
#include <protomol/integrator/hessian/HessianInt.h>
#include <protomol/module/HessianIntegratorModule.h>

#include "protomol/factory/IntegratorFactory.h"

using namespace std;
using namespace ProtoMol;

void HessianIntegratorModule::init(ProtoMolApp *app) {
  app->integratorFactory.registerExemplar(new HessianInt());
}
