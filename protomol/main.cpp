#include <protomol/ProtoMolApp.h>           // for ProtoMolApp
#include <protomol/base/Exception.h>        // for operator<<, Exception
#include <protomol/base/ModuleManager.h>    // for ModuleManager
#include <protomol/module/MainModule.h>     // for InputDebug
#include <iostream>                         // for cout, operator<<, ostream, basic_ostream, cerr, endl

#include "protomol/base/Report.h"
#include "protomol/base/Report.h"           // for cerr
#include "protomol/config/Configuration.h"
#include "protomol/config/Configuration.h"  // for Configuration
#include "protomol/config/Value.h"
#include "protomol/config/Value.h"          // for Value

using namespace std;
using namespace ProtoMol;

extern void moduleInitFunction(ModuleManager *);

int main(int argc, char *argv[]) {
  try {
    ModuleManager modManager;
    moduleInitFunction(&modManager);
    ProtoMolApp app(&modManager);

    if (!app.configure(argc, argv)) return 0;
    app.splash(cout);
    app.build();
    if ((int)app.config[InputDebug::keyword]) app.print(cout);

    while (app.step()) continue;
    app.finalize();

    return 0;

  } catch (const Exception &e) {
    cerr << "ERROR: " << e << endl;
  }

  return 1;
}
