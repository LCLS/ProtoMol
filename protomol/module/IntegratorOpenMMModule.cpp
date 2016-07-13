#include <protomol/module/IntegratorOpenMMModule.h>

#ifdef HAVE_OPENMM
#include <protomol/integrator/openMM/OpenMMIntegrator.h>

#ifdef HAVE_OPENMM_LTMD
#include <protomol/integrator/openMM/NormalModeOpenMM.h>
#endif

#ifdef HAVE_OPENMM_FBM
#include <protomol/integrator/openMM/OpenMMFBMIntegrator.h>
#endif

#endif

#include <protomol/ProtoMolApp.h>

#include "protomol/factory/IntegratorFactory.h"

using namespace std;
using namespace ProtoMol;

void IntegratorOpenMMModule::init(ProtoMolApp *app) {
#ifdef HAVE_OPENMM
  app->integratorFactory.registerExemplar(new OpenMMIntegrator());

#ifdef HAVE_OPENMM_LTMD
  app->integratorFactory.registerExemplar(new NormalModeOpenMM());
#endif

#ifdef HAVE_OPENMM_FBM
  app->integratorFactory.registerExemplar(new OpenMMFBMIntegrator());
#endif

#endif
}
