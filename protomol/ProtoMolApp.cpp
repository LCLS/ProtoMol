#include <protomol/ProtoMolApp.h>

#include <protomol/base/ModuleManager.h>
#include <protomol/base/SystemUtilities.h>
#include <protomol/base/StringUtilities.h>
#include <protomol/base/PMConstants.h>
#include <protomol/base/TimerStatistic.h>
#include <protomol/base/Zap.h>
#include <protomol/base/Report.h>

#include <protomol/module/MainModule.h>
#include <protomol/module/IOModule.h>
#include <protomol/module/ConfigurationModule.h>
#include <protomol/module/AnalysisModule.h>

#include <protomol/type/String.h>

#include <protomol/config/CommandLine.h>
#include <protomol/config/Configuration.h>

#include <protomol/io/ConfigurationReader.h>

#include <protomol/factory/TopologyFactory.h>
#include <protomol/factory/OutputFactory.h>

#include <protomol/topology/GenericTopology.h>
#include <protomol/topology/BuildTopology.h>
#include <protomol/topology/BuildTopologyFromTpr.h>
#include <protomol/topology/BuildTopologyFromXML.h>
#include <protomol/topology/TopologyUtilities.h>

#include <protomol/output/OutputCollection.h>
#include <protomol/analysis/AnalysisCollection.h>

#include <protomol/parallel/Parallel.h>

#include <iomanip>
#ifdef HAVE_PACKAGE_H
#include <protomol/package.h>
#endif

#ifdef HAVE_LIBFAH
#include <fah/core/Core.h>
#endif

using namespace std;
using namespace ProtoMol;
using namespace ProtoMol::Report;
using namespace ProtoMol::Constant;


ProtoMolApp::ProtoMolApp(ModuleManager *modManager) :
  modManager(modManager), SCPISMParameters(0), cmdLine(&config), outputs(0),
  integrator(0), topology(0) {

  TimerStatistic::timer[TimerStatistic::WALL].start();

  modManager->init(this);

  topologyFactory.registerAllExemplarsConfiguration(&config);
  outputFactory.registerAllExemplarsConfiguration(&config);
  analysisFactory.registerAllExemplarsConfiguration(&config);
}


ProtoMolApp::~ProtoMolApp() {}


void ProtoMolApp::splash(ostream &stream) {
	if( Parallel::iAmMaster() ){
		const int w = 16;
		stream
		<< headerRow("ProtoMol") << endl
		<< setw(w) << "Description: ";
		fillFormat(stream, "A rapid PROTOtyping MOLecular dynamics object-oriented "
				 "component based framework.", w, w);
		stream
#ifdef HAVE_PACKAGE_H
		<< setw(w) << "Version: " << PACKAGE_VERSION << endl
		<< setw(w) << "SVN revision: " << PACKAGE_REVISION << endl
		<< setw(w) << "Repository: " << PACKAGE_SOURCE_REPO << endl
		<< setw(w) << "Homepage: " << PACKAGE_HOMEPAGE << endl
		<< setw(w) << "Report bugs to: " << PACKAGE_BUGREPORT << endl
		<< setw(w) << "Compiler: " << PACKAGE_COMPILER << " "
		<< PACKAGE_COMPILER_VERSION << endl
		<< setw(w) << "Flags: " << PACKAGE_COMPILER_FLAGS << endl
		<< setw(w) << "Extra libs: " << PACKAGE_COMPILER_LIBS << endl
		<< setw(w) << "Built by: " << PACKAGE_BUILT_BY << endl
		<< setw(w) << "Build platform: " << PACKAGE_PLATFORM << endl
#endif // HAVE_PACKAGE_H
		<< setw(w) << "Build date: " <<  __DATE__ << ", " << __TIME__ << endl
#ifdef HAVE_LIBFAH
		<< setw(w) << "Checksumming: "
		<< "Enabled for Folding@Home file protection." << endl
#endif // HAVE_LIBFAH
    << setw(w) << "Long/Int size: " << sizeof(long) << "/" << sizeof(int) << endl
#ifdef HAVE_PACKAGE_H
		<< setw(w) << "Please cite: ";
		fillFormat(stream, PACKAGE_CITE, w, w);
		stream
#endif // HAVE_PACKAGE_H
		<< PROTOMOL_HR << endl;
	}
}


void ProtoMolApp::load(const string &configfile) {
  vector<string> args;

  args.push_back("ProtoMol");
  args.push_back(configfile);
  load(args);
}


bool ProtoMolApp::load(int argc, char *argv[]) {
	Parallel::init( argc, argv );
  return load(vector<string>(argv, argv + argc));
}


bool ProtoMolApp::load(const vector<string> &args) {
  // Parse command line
  if (cmdLine.parse(args)) return false;

  return true;
}


void ProtoMolApp::configure() {
  // Read Config file
  if (config.valid(InputConfig::keyword))
    SystemUtilities::chdir(SystemUtilities::dirname
                           (config[InputConfig::keyword]));
  else THROW("Configuration file not set.");

  modManager->configure(this);
}


void ProtoMolApp::configure(const string &configfile) {
  load(configfile);
  configure();
}


bool ProtoMolApp::configure(int argc, char *argv[]) {
  if (!load(argc, argv)) return false;
  configure();
  return true;
}


bool ProtoMolApp::configure(const vector<string> &args) {
  if (!load(args)) return false;

  configure();

  return true;
}


void ProtoMolApp::build() {

  //Set parallel parameters
  if( Parallel::isParallel() ){
    //mode
    if (config[InputParallelMode::keyword].valid())
      Parallel::setMode(config[InputParallelMode::keyword].getString());
    //barrier
    if (config[InputUseBarrier::keyword].valid())
      Parallel::setBarrier(config[InputUseBarrier::keyword]);
    //pipe size
    if (config[InputParallelPipe::keyword].valid())
      Parallel::setPipeSize(config[InputParallelPipe::keyword]);
    //max pacjkages
    if (config[InputMaxPackages::keyword].valid())
      Parallel::setMaxPackages(config[InputMaxPackages::keyword]);

  }

  // TPR input for topology, positions and velocities?
  // Then check for Gromacs support
#if !defined(HAVE_GROMACS)
  if (config.valid(InputGromacsTprFile::keyword) && config.valid(InputGromacsXMLFile::keyword))
    THROWS("GROMACS support not available for '" <<
           config[InputGromacsTprFile::keyword].getString() << "'.");
#endif

  //flags for GROMACS TPR/XML
  bool GROMACRTPR(false);
  bool GROMACSXML(false);
  bool GROMACSNEWPOSITIONS(false);
  bool GROMACSNEWVELOCITIES(false);

  //test TPR file
  if( config.valid(InputGromacsTprFile::keyword)){
    GROMACRTPR = true;
    if( config.valid("Checkpoint") ) {
        GROMACSNEWPOSITIONS = true;
        GROMACSNEWVELOCITIES = true;
    }
    if( config.valid(InputPositions::keyword) ) GROMACSNEWPOSITIONS = true;
    if( config.valid(InputVelocities::keyword) ) GROMACSNEWVELOCITIES = true;
  }

  //test TPR file
  if( config.valid(InputGromacsXMLFile::keyword)){
    GROMACSXML = true;
    GROMACSNEWPOSITIONS = true;
    GROMACSNEWVELOCITIES = true;
    if( !config.valid(InputPositions::keyword) ){
      THROWS("GROMACS PDB input file required for XML input!");
    }
  }
  
  // Read data if not TPR unless positions defined
  if( !GROMACRTPR || GROMACSNEWPOSITIONS || GROMACSNEWVELOCITIES) {
    modManager->read(this);
    GROMACSNEWPOSITIONS = config.valid(InputPositions::keyword);
    GROMACSNEWVELOCITIES = config.valid(InputVelocities::keyword);
  }

  // Build topology
  try {
    topology = topologyFactory.make(&config);
  } catch (const Exception &e) {

    if( !GROMACRTPR && !GROMACSXML){
      // Try to get some defaults with the postions known ...
      const GenericTopology *prototype =
        topologyFactory.find(config[GenericTopology::keyword].getString());

      if (prototype) {
        vector<Parameter> parameters = prototype->getDefaults(positions);

        for (unsigned int i = 0; i < parameters.size(); i++)
          if (!config.valid(parameters[i].keyword) &&
              parameters[i].value.valid()) {
            config.set(parameters[i].keyword, parameters[i].value);
            report << hint << parameters[i].keyword << " undefined, using "
                   << parameters[i].value.getString() << "." << endr;
          }

        topology = topologyFactory.make(&config);
      }
    }

    if (!topology) throw e;
  }

  if (!GROMACRTPR && !GROMACSXML) {
    // Using SCPISM parameter? Flag or filename
    if (config[InputDoSCPISM::keyword] || SCPISMParameters) {

      if(config[InputDoSCPISM::keyword])
        topology->doSCPISM = config[InputDoSCPISM::keyword];

      if ((topology->doSCPISM < 1 || topology->doSCPISM > 3) &&
          !SCPISMParameters)
        THROW("doscpism should be between 1 and 3 or an input file should be "
              "used.");

      if  (SCPISMParameters) {
        if (!config[InputDoSCPISM::keyword]) topology->doSCPISM = 4;
      } else {
        SCPISMParameters = new CoulombSCPISMParameterTable;
        SCPISMParameters->populateTable();
      }

      report
        << "SCPISM: doSCPISM set to " << topology->doSCPISM << "." << endr;

      if (topology->doSCPISM == 3) {
        // Quartic switch parameters
        SCPISMParameters->myData["H"].hbond_factor = 0.4695;
        SCPISMParameters->myData["HC"].hbond_factor = 7.2560;
      }

      SCPISMParameters->displayTable();

    }

    // Using GBSA with openMM
    if (config[InputDoGBSAObc::keyword]) {
      topology->doGBSAOpenMM = 1;
      topology->obcType = config[InputDoGBSAObc::keyword];

      if (topology->obcType == 1) {
        topology->alphaObc = 0.8;
        topology->betaObc = 0;
        topology->gammaObc = 2.91;

      } else if (topology->obcType == 2) {
        //obctype = 2
        topology->alphaObc = 1.0;
        topology->betaObc = 0.8;
        topology->gammaObc = 4.85;
      }
      topology->dielecOffset = 0.09;

    }

    // Find force field type before building topology
    if (config.valid(InputGromacsTopo::keyword) &&
        config.valid(InputGromacsParamPath::keyword)) {
      topology->forceFieldFlag = GROMACS;
    }

    // Build the topology
    buildTopology(topology, psf, par, config[InputDihedralMultPSF::keyword],
                  SCPISMParameters);

  } else {
    //TPR/GROMACS if here
    topology->forceFieldFlag = GROMACS;//TPR or XML

    if(!GROMACSXML){  //if not XML then build topology from TPR
      // Build the topology from the tpr file
      buildTopologyFromTpr( topology, positions, velocities,
                            config[InputGromacsTprFile::keyword].getString(),
                             GROMACSNEWPOSITIONS, GROMACSNEWVELOCITIES);
    }else{  //else from XML
      buildTopologyFromXML( topology, positions, velocities,
                           config[InputGromacsXMLFile::keyword].getString(),
                           config[InputPositions::keyword].getString(),
                           outputCache.getAtoms());
    }
  }

  // Register Forces
  modManager->registerForces(this);

  // Build the integrators and forces
  integrator =
    integratorFactory.make(config[InputIntegrator::keyword], &forceFactory);

  // Update number of steps based on a provided total simulation time
  if( config[InputSimulationTime::keyword].valid() ) {
    std::string sTime = config[InputSimulationTime::keyword];
	std::transform(sTime.begin(), sTime.end(), sTime.begin(), ::tolower);

    int rTime;
	char sUnit[2];
    sscanf( sTime.c_str(), "%d", &rTime );

	// Is it an integer
	if( sTime[sTime.length()-2] > 47 && sTime[sTime.length()-2] < 58 ){
		sUnit[0] = sTime[sTime.length()-1];
	}else{
		sUnit[0] = sTime[sTime.length()-2];
		sUnit[1] = sTime[sTime.length()-1];
	}

    Real SimTime = 0.0f;
    if( strncmp(sUnit, "fs", 2) == 0 ) SimTime = rTime;
    if( strncmp(sUnit, "ps", 2) == 0 ) SimTime = rTime * 1e3;
    if( strncmp(sUnit, "ns", 2) == 0 ) SimTime = rTime * 1e6;
    if( strncmp(sUnit, "us", 2) == 0 ) SimTime = rTime * 1e9;
    if( strncmp(sUnit, "ms", 2) == 0 ) SimTime = rTime * 1e12;
    if( strncmp(sUnit, "s", 1) == 0 ) SimTime = rTime * 1e15;

	// Update number of steps
    config[InputNumsteps::keyword] = int(SimTime / integrator->getTimestep());
  }

  // Update number of steps between outputs based on a provided output step time
  if( config[InputOutputTime::keyword].valid() ) {
    std::string sTime = config[InputOutputTime::keyword];
	std::transform(sTime.begin(), sTime.end(), sTime.begin(), ::tolower);

    int rTime;
	char sUnit[2];
    sscanf( sTime.c_str(), "%d", &rTime );

	// Is it an integer
	if( sTime[sTime.length()-2] > 47 && sTime[sTime.length()-2] < 58 ){
		sUnit[0] = sTime[sTime.length()-1];
	}else{
		sUnit[0] = sTime[sTime.length()-2];
		sUnit[1] = sTime[sTime.length()-1];
	}

    Real SimTime = 0.0f;
    if( strncmp(sUnit, "fs", 2) == 0 ) SimTime = rTime;
    if( strncmp(sUnit, "ps", 2) == 0 ) SimTime = rTime * 1e3;
    if( strncmp(sUnit, "ns", 2) == 0 ) SimTime = rTime * 1e6;
    if( strncmp(sUnit, "us", 2) == 0 ) SimTime = rTime * 1e9;
    if( strncmp(sUnit, "ms", 2) == 0 ) SimTime = rTime * 1e12;
    if( strncmp(sUnit, "s", 1) == 0 ) SimTime = rTime * 1e15;

	// Update number of steps
    config[InputOutputfreq::keyword] = int(SimTime / integrator->getTimestep());
  }

  // Setup run parameters (used for GUI so required here)
  currentStep = config[InputFirststep::keyword];
  lastStep = config[InputNumsteps::keyword];
  if( lastStep != -1 ){
    lastStep += currentStep;
  }

  // Create outputs
  // TODO if !Parallel::iAmMaster() turn off outputs
	if( !Parallel::iAmMaster() ){
		outputs = new OutputCollection;
        analysis = new AnalysisCollection;
	}else{
		if (config[InputOutput::keyword]){
			outputs = outputFactory.makeCollection(&config);
		}

        if (config[InputAnalysis::keyword]){
            analysis = analysisFactory.makeCollection(&config);
        }
	}

  // Post build processing
  modManager->postBuild(this);

  report << plain << "Actual start temperature : "
         << temperature(topology, &velocities) << "K" << endr;


  // Add Integrator Modifiers
  modManager->addModifiers(this);

  // Initialize
  energies.molecularVirial(config[InputMolVirialCalc::keyword]);
  energies.virial(config[InputVirialCalc::keyword]);
  report << plain << "Virial tensor : " << energies.virial() << endr;
  report << plain << "Molecular virial tensor : "
         << energies.molecularVirial() << endr;

  topology->time = (Real)(toInt(config[InputFirststep::keyword]) * integrator->getTimestep());

  integrator->initialize(this);
  outputs->initialize(this);
  outputCache.initialize(this);

  // Init cache
  //outputs->addToCache(pdbAtoms); // TODO fix this
  outputCache.add(psf);
  outputCache.add(par);

  // Initialize Analysis
  analysis->initialize(this);

  // Print Factories
  if ((int)config[InputDebug::keyword] >= 5 && (int)config[InputDebugLimit::keyword] <= 5)
    cout
      << headerRow("Factories")     << endl
      << headerRow("Configuration") << endl << config            << endl
      << headerRow("Topology")      << endl << topologyFactory   << endl
      << headerRow("Integrator")    << endl << integratorFactory << endl
      << headerRow("Force")         << endl << forceFactory      << endl
      << headerRow("Output")        << endl << outputFactory     << endl
      << headerRow("Analysis")      << endl << analysisFactory   << endl;

  // Clear all factories
  topologyFactory.unregisterAllExemplars();
  integratorFactory.unregisterAllExemplars();
  forceFactory.unregisterAllExemplars();
  outputFactory.unregisterAllExemplars();
  analysisFactory.unregisterAllExemplars();

  TimerStatistic::timer[TimerStatistic::RUN].reset();
  TimerStatistic::timer[TimerStatistic::INTEGRATOR].reset();
  TimerStatistic::timer[TimerStatistic::FORCES].reset();
  TimerStatistic::timer[TimerStatistic::COMMUNICATION].reset();
  TimerStatistic::timer[TimerStatistic::IDLE].reset();
}


bool ProtoMolApp::step(long inc) {
  if (currentStep >= lastStep && lastStep != -1) return false;

  TimerStatistic::timer[TimerStatistic::RUN].start();

  if (outputs->run(currentStep)) {
#ifdef HAVE_LIBFAH
    // Make sure we save the latest checksum information after writing data.
    if (FAH::Core::isActive()) FAH::Core::instance().checkpoint();
#endif
  }

  // Handle analysis on each output step
  analysis->run(currentStep);

  if (!inc) inc = outputs->getNext() - currentStep;

  //fix inc so do not overrun
  if( lastStep != -1 ){
    inc = std::min(lastStep, currentStep + inc) - currentStep;
  }

  TimerStatistic::timer[TimerStatistic::INTEGRATOR].start();

  const long completed = integrator->run(inc);

  TimerStatistic::timer[TimerStatistic::INTEGRATOR].stop();

  // moved here so that current step is valid in integrator
  if( completed != inc ) {
    report << plain << "Did not complete all steps: " << completed << " of " << inc << std::endl;
  }
  currentStep += completed;

  if( analysis->shouldStop() ) lastStep = currentStep;

  TimerStatistic::timer[TimerStatistic::RUN].stop();

  return true;
}


void ProtoMolApp::finalize() {
  float nanoSeconds = ( (currentStep - config[InputFirststep::keyword].operator long()) * integrator->getTimestep() ) / 1000000;
  float secondsPerNanosecond = TimerStatistic::timer[TimerStatistic::RUN].getTime().getRealTime() / nanoSeconds;
  float nanoSecondsPerDay = 86400 / secondsPerNanosecond;

  report.setf( std::ios::fixed );
  report << plain << "Performance: " << nanoSeconds << "ns in " << TimerStatistic::timer[TimerStatistic::RUN].getTime().getRealTime() << "s = " << nanoSecondsPerDay << "ns/day" << std::endl;

  outputs->finalize(currentStep);
  analysis->finalize(currentStep);

  // Clean up
  zap(topology);
  zap(integrator);
  zap(outputs);
  zap(analysis);
  zap(SCPISMParameters);

  TimerStatistic::timer[TimerStatistic::WALL].stop();

	if( Parallel::iAmMaster() ){
		report << plain << "Timing: " << TimerStatistic() << "." << endr;
	}

	Parallel::finalize();
}


void ProtoMolApp::print(ostream &stream) {
	if( Parallel::iAmMaster() ){
		// Output
		stream << headerRow("Outputs") << endl;

		for (OutputCollection::const_iterator itr =
			 const_cast<const OutputCollection *>(outputs)->begin();
		   itr != const_cast<const OutputCollection*>(outputs)->end(); itr++) {

		stream << "Output " << (*itr)->getId();

		vector<Parameter> parameters;
		(*itr)->getParameters(parameters);
		for (unsigned int i = 0; i < parameters.size(); i++)
		  stream << " " << parameters[i].value.getString();
		stream << "." << endl;
		}

		if (!((bool)config[InputOutput::keyword]))
		stream << "All output suppressed!" << endl;


		// Integrator
		stream << headerRow("Integrator") << endl;
		vector<IntegratorDefinition> inter = integrator->getIntegratorDefinitionAll();
		stream  << InputIntegrator::keyword << " {" << endl;

		for (int i = inter.size() - 1; i >= 0; i--)
		stream << Constant::PRINTINDENT << "Level "
			   << i << " " << inter[i].print() << endl;

		stream << "}" << endl;


		// Topology
		stream << headerRow("Topology") << endl;
		stream << topology->print(&positions) << endl;

		stream << PROTOMOL_HR << endl;
	}
}
