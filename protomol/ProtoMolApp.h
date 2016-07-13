#ifndef PROTOMOLAPP_H
#define PROTOMOLAPP_H

#include <protomol/config/CommandLine.h>         // for CommandLine
#include <protomol/config/Configuration.h>       // for Configuration
#include <protomol/factory/AnalysisFactory.h>    // for AnalysisFactory
#include <protomol/factory/ForceFactory.h>       // for ForceFactory
#include <protomol/factory/IntegratorFactory.h>  // for IntegratorFactory
#include <protomol/factory/OutputFactory.h>      // for OutputFactory
#include <protomol/factory/TopologyFactory.h>    // for TopologyFactory
#include <protomol/output/OutputCache.h>         // for OutputCache
#include <protomol/type/EigenvectorInfo.h>       // for EigenvectorInfo
#include <protomol/type/PAR.h>                   // for PAR
#include <protomol/type/PSF.h>                   // for PSF
#include <protomol/type/ScalarStructure.h>       // for ScalarStructure
#include <protomol/type/Vector3DBlock.h>         // for Vector3DBlock
#include <ostream>                               // for ostream
#include <string>                                // for string
#include <vector>                                // for vector

namespace ProtoMol {
  class OutputCollection;
  class AnalysisCollection;
  class Integrator;
  class GenericTopology;
  class ModuleManager;
  struct CoulombSCPISMParameterTable;

  class ProtoMolApp {
  public:
    ModuleManager *modManager;

    // Data
    Vector3DBlock positions;
    Vector3DBlock velocities;
    EigenvectorInfo eigenInfo;
    PSF psf;
    PAR par;
    ScalarStructure energies;
    CoulombSCPISMParameterTable *SCPISMParameters;

    // Factories
    TopologyFactory topologyFactory;
    ForceFactory forceFactory;
    IntegratorFactory integratorFactory;
    OutputFactory outputFactory;
    AnalysisFactory analysisFactory;

    // Containers
    CommandLine cmdLine;
    Configuration config;
    mutable OutputCache outputCache;
    OutputCollection *outputs;
    AnalysisCollection *analysis;
    Integrator *integrator;
    GenericTopology *topology;

    // Run
    long currentStep;
    long lastStep;

    ProtoMolApp() {}
    ProtoMolApp(ModuleManager *modManager);
    ~ProtoMolApp();

    static void splash(std::ostream &stream);
    void load(const std::string &configfile);
    bool load(int argc = 0, char *argv[] = 0);
    bool load(const std::vector<std::string> &args);
    void configure();
    void configure(const std::string &configfile);
    bool configure(int argc, char *argv[]);
    bool configure(const std::vector<std::string> &args);
    void build();
    void print(std::ostream &stream);
    bool step(long inc = 0);
    void finalize();
  };
}

#endif // PROTOMOLAPP_H
