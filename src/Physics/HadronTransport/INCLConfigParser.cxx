#include "Framework/Conventions/GBuild.h"
#ifdef __GENIE_INCL_ENABLED__

/** \file ConfigParser.cxx
 * \brief A class for parsing input and producing a INCL++ Config
 *        Copied from INC++ main
 * see ConfigParser.h for more information
 *
 * \date 17th July 2014
 * \author Davide Mancusi
 */

#ifdef HAS_BOOST_PROGRAM_OPTIONS
//#include "ConfigParser.hh"  // this would be INCL copy of header
#include "INCLConfigParser.h"  // GENIE copy

#include "G4INCLLogger.hh"
#include "G4INCLParticleTable.hh"
#include "DatafilePaths.hh"

#include <fstream>
#include <cstdlib>
#include <sstream>
#include <cerrno>
#include <cstdlib>

const std::string ConfigParser::suggestHelpMsg = "You might want to run `INCLCascade --help' to get a help message.\n";

// Define the names of the de-excitation models
const std::string ConfigParser::theNoneName = "none";
#ifdef INCL_DEEXCITATION_SMM
const std::string ConfigParser::theSMMName = "SMM";
#endif
#ifdef INCL_DEEXCITATION_GEMINIXX
const std::string ConfigParser::theGEMINIXXName = "GEMINIXX";
#endif
#ifdef INCL_DEEXCITATION_ABLAXX
const std::string ConfigParser::theABLAv3pName = "ABLAv3p";
#endif
#ifdef INCL_DEEXCITATION_ABLA07
const std::string ConfigParser::theABLA07Name = "ABLA07";
#endif

const std::string ConfigParser::listSeparator = "\n  \t";

ConfigParser::ConfigParser() :
  runOptDesc("Run options"),
  hiddenOptDesc("Hidden options"),
  genericOptDesc("Generic options"),
  physicsOptDesc("Physics options"),
  deExcitationModelList(
                        listSeparator + theNoneName
#ifdef INCL_DEEXCITATION_ABLA07
                        + listSeparator + theABLA07Name
#endif
#ifdef INCL_DEEXCITATION_ABLAXX
                        + listSeparator + theABLAv3pName
#endif
#ifdef INCL_DEEXCITATION_GEMINIXX
                        + listSeparator + theGEMINIXXName
#endif
#ifdef INCL_DEEXCITATION_SMM
                        + listSeparator + theSMMName
#endif
                       )
{
  // Define the default de-excitation model, in decreasing order of priority
  defaultDeExcitationModel = theNoneName;
#ifdef INCL_DEEXCITATION_SMM
  defaultDeExcitationModel = theSMMName;
#endif
#ifdef INCL_DEEXCITATION_GEMINIXX
  defaultDeExcitationModel = theGEMINIXXName;
#endif
#ifdef INCL_DEEXCITATION_ABLAXX
  defaultDeExcitationModel = theABLAv3pName;
#endif
#ifdef INCL_DEEXCITATION_ABLA07
  defaultDeExcitationModel = theABLA07Name;
#endif

  // Append " (default)" to the name of the default model
  size_t defaultModelIndex = deExcitationModelList.find(defaultDeExcitationModel);
  if(defaultModelIndex!=std::string::npos) {
    deExcitationModelList = deExcitationModelList.substr(0, defaultModelIndex+defaultDeExcitationModel.size())
      + " (default)"
      + deExcitationModelList.substr(defaultModelIndex+defaultDeExcitationModel.size(), std::string::npos);
  }

  // Hidden options
  hiddenOptDesc.add_options()
    ("input-file", po::value<std::string>(&config.inputFileName), "input file")
    ("impact-parameter", po::value<double>(&config.impactParameter)->default_value(-1.), "impact parameter")
    ("cascade-action", po::value<std::string>(&config.cascadeAction)->default_value("default"), "cascade action:\n  \tdefault (default)\n  \tavatar-dump")
    ;

  // Generic options
  std::stringstream verbosityDescription;
  verbosityDescription << "set verbosity level:\n"
    << " 0: \tquiet, suppress all output messages\n"
    << " " << G4INCL::InfoMsg << ": \tminimal logging\n"
    << " " << G4INCL::FatalMsg << ": \tlog fatal error messages as well\n"
    << " " << G4INCL::ErrorMsg << ": \tlog error messages as well\n"
    << " " << G4INCL::WarningMsg << ": \tlog warning messages as well\n"
    << " " << G4INCL::DebugMsg << ": \tlog debug messages as well\n"
    << " " << G4INCL::DataBlockMsg << ": \tlog data-block messages as well";

  genericOptDesc.add_options()
    ("help,h", "produce this help message")
    ("version", "print version string and exit")
    ;

  // Run-specific options
  std::stringstream randomSeedsDescription;
  randomSeedsDescription << "comma-separated list of seeds for the random-number generator. Allowed seed range: "
    << randomSeedMin << "-" << randomSeedMax << ".";

  runOptDesc.add_options()
    ("title", po::value<std::string>(&config.title)->default_value("INCL default run title"), "run title")
    ("output,o", po::value<std::string>(&config.outputFileRoot), "root for generating output file names. File-specific suffixes (.root, .out, etc.) will be appended to this root. Defaults to the input file name, if given; otherwise, defaults to a string composed of the explicitly specified options and of a customisable suffix, if provided using the -s option")
    ("suffix,s", po::value<std::string>(&config.fileSuffix), "suffix to be appended to generated output file names")
    ("logfile,l", po::value<std::string>(&config.logFileName), "log file name. Defaults to `<output_root>.log'. Use `-' if you want to redirect logging to stdout")
    ("number-shots,N", po::value<int>(&config.nShots), "* number of shots")
    ("target,t", po::value<std::string>(&config.targetString), "* target nuclide. Can be specified as Fe56, 56Fe, Fe-56, 56-Fe, Fe_56, 56_Fe or Fe. If the mass number is omitted, natural target composition is assumed.")
    ("projectile,p", po::value<std::string>(&config.projectileString), "* projectile name:\n"
     "  \tproton, p\n"
     "  \tneutron, n\n"
     "  \tpi+, piplus, pion+, pionplus\n"
     "  \tpi0, pizero, pion0, pionzero\n"
     "  \tpi-, piminus, pion-, pionminus\n"
     "  \td, t, a, deuteron, triton, alpha\n"
     "  \tHe-4, He4, 4He (and so on)")
    ("energy,E", po::value<double>(&config.projectileKineticEnergy), "* total kinetic energy of the projectile, in MeV")
    ("verbose-event", po::value<int>(&config.verboseEvent)->default_value(-1), "request verbose logging for the specified event only")
    ("random-number-generator", po::value<std::string>(&config.randomNumberGenerator)->default_value("Ranecu"), "Random number generator to use:\n  \tRanecu (2 seeds, default)\n  \tRanecu3 (3 seeds)")
    ("random-seeds", po::value<std::string>(&config.randomSeeds)->default_value("666,777,1234"), randomSeedsDescription.str().c_str())
    ("autosave-frequency", po::value<unsigned int>(&config.autosaveFrequency)->default_value(10000), "frequency between automatic saves of the output/log files")
#ifdef INCL_ROOT_USE
    ("root-selection", po::value<std::string>(&config.rootSelectionString)->default_value(""), "ROOT selection for abridged output ROOT tree. For example: \"A==1 && Z==0 && theta<3\" selects only events where a neutron is scattered in the forward direction.")
    ("concise-root-tree", po::value<bool>(&config.conciseROOTTree)->default_value(false), "whether INCL++ should output a concise ROOT event tree:\n  \ttrue, 1\n  \tfalse, 0 (default)")
#endif
    ("inverse-kinematics", po::value<bool>(&config.inverseKinematics)->default_value(false), "whether INCL++ should output variables describing the reaction in inverse kinematics:\n  \ttrue, 1\n  \tfalse, 0 (default)")
    ("inclxx-datafile-path", po::value<std::string>(&config.INCLXXDataFilePath)->default_value(defaultINCLXXDatafilePath),
     "path to the INCL++ data files")
#ifdef INCL_DEEXCITATION_ABLA07
    ("abla07-datafile-path", po::value<std::string>(&config.abla07DataFilePath)->default_value(defaultABLA07DatafilePath),
     "path to the ABLA07 data files")
#endif
#ifdef INCL_DEEXCITATION_ABLAXX
    ("ablav3p-cxx-datafile-path", po::value<std::string>(&config.ablav3pCxxDataFilePath)->default_value(defaultABLAXXDatafilePath),
     "path to the ABLAv3p data files")
#endif
#ifdef INCL_DEEXCITATION_GEMINIXX
    ("geminixx-datafile-path", po::value<std::string>(&config.geminixxDataFilePath)->default_value(defaultGEMINIXXDatafilePath),
     "path to the GEMINI++ data files")
#endif
    ("verbosity,v", po::value<int>(&config.verbosity)->default_value(4), verbosityDescription.str().c_str())
    ;

  // Physics options
  physicsOptDesc.add_options()
    ("de-excitation,d", po::value<std::string>(&config.deExcitationString)->default_value(defaultDeExcitationModel.c_str()), ("which de-excitation model to use:" + deExcitationModelList).c_str())
#ifdef INCL_DEEXCITATION_FERMI_BREAKUP
    ("max-mass-fermi-breakup", po::value<int>(&config.maxMassFermiBreakUp)->default_value(16), "Maximum remnant mass for Fermi break-up (0-16). Default: 16.")
    ("max-charge-fermi-breakup", po::value<int>(&config.maxChargeFermiBreakUp)->default_value(8), "Maximum remnant mass for Fermi break-up (0-8). Default: 8.")
#endif
    ("pauli", po::value<std::string>(&config.pauliString)->default_value("strict-statistical"), "Pauli-blocking algorithm:\n"
     "  \tstrict-statistical (default)\n"
     "  \tstrict\n"
     "  \tstatistical\n"
     "  \tglobal\n"
     "  \tnone")
    ("cdpp", po::value<bool>(&config.CDPP)->default_value(true), "whether to apply CDPP after collisions:\n  \ttrue, 1 (default)\n  \tfalse, 0")
    ("coulomb", po::value<std::string>(&config.coulombString)->default_value("non-relativistic"), "Coulomb-distortion algorithm:\n  \tnon-relativistic (default)\n  \tnone")
    ("potential", po::value<std::string>(&config.potentialString)->default_value("isospin-energy"), "nucleon potential:\n  \tisospin-energy-smooth\n  \tisospin-energy (default)\n  \tisospin\n  \tconstant")
    ("pion-potential", po::value<bool>(&config.pionPotential)->default_value("true"), "whether to use a pion potential:\n  \ttrue, 1 (default)\n  \tfalse, 0")
    ("local-energy-BB", po::value<std::string>(&config.localEnergyBBString)->default_value("first-collision"), "local energy in baryon-baryon collisions:\n  \talways\n  \tfirst-collision (default)\n  \tnever")
    ("local-energy-pi", po::value<std::string>(&config.localEnergyPiString)->default_value("first-collision"), "local energy in pi-N collisions and in delta decays:\n  \talways\n  \tfirst-collision (default)\n  \tnever")
    ("cluster-algorithm", po::value<std::string>(&config.clusterAlgorithmString)->default_value("intercomparison"), "clustering algorithm for production of composites:\n  \tintercomparison (default)\n  \tnone")
    ("cluster-max-mass", po::value<int>(&config.clusterMaxMass)->default_value(8), "maximum mass of produced composites:\n  \tminimum 2\n  \tmaximum 12")
    ("back-to-spectator", po::value<bool>(&config.backToSpectator)->default_value("true"), "whether to use back-to-spectator:\n  \ttrue, 1 (default)\n  \tfalse, 0")
    ("use-real-masses", po::value<bool>(&config.useRealMasses)->default_value("true"), "whether to use real masses for the outgoing particle energies:\n  \ttrue, 1 (default)\n  \tfalse, 0")
    ("separation-energies", po::value<std::string>(&config.separationEnergyString)->default_value("INCL"), "how to assign the separation energies of the INCL nucleus:\n  \tINCL (default)\n  \treal\n  \treal-light")
    ("fermi-momentum", po::value<std::string>(&config.fermiMomentumString)->default_value("constant"), "how to assign the Fermi momentum of the INCL nucleus:\n  \tconstant (default)\n  \tconstant-light\n  \tmass-dependent\n  \t[a positive value]")
    ("cutNN", po::value<double>(&config.cutNN)->default_value(1910.), "minimum CM energy for nucleon-nucleon collisions, in MeV. Default: 1910.")
    ("rp-correlation", po::value<double>(&config.rpCorrelationCoefficient)->default_value(1.), "correlation coefficient for the r-p correlation. Default: 1 (full correlation).")
    ("rp-correlation-p", po::value<double>(&config.rpCorrelationCoefficientProton)->default_value(1.), "correlation coefficient for the proton r-p correlation. Overrides the value specified using the rp-correlation option. Default: 1 (full correlation).")
    ("rp-correlation-n", po::value<double>(&config.rpCorrelationCoefficientNeutron)->default_value(1.), "correlation coefficient for the neutron r-p correlation. Overrides the value specified using the rp-correlation option. Default: 1 (full correlation).")
    ("neutron-skin", po::value<double>(&config.neutronSkin)->default_value(0.), "thickness of the neutron skin, in fm. Default: 0.")
    ("neutron-halo", po::value<double>(&config.neutronHalo)->default_value(0.), "thickness of the neutron halo, in fm. Default: 0.")
    ("refraction", po::value<bool>(&config.refraction)->default_value(false), "whether to use refraction when particles are transmitted. Default: false.")
    ("cross-sections", po::value<std::string>(&config.crossSectionsString)->default_value("multipions"), "cross-section parametrizations:\n"
     "  \tmultipions (default)\n"
     "  \ttruncated-multipions\n"
     "  \tincl46")
    ("max-number-multipions", po::value<int>(&config.maxNumberMultipions)->default_value(-1), "maximum number of pions that can be produced in multipion collisions. When specified with arg>0, it enforces cross-sections=truncated-multipions. Default: -1 (no limit).")
    ("phase-space-generator", po::value<std::string>(&config.phaseSpaceGenerator)->default_value("raubold-lynch"), "algorithm to generate phase-space decays:\n  \tRaubold-Lynch (default)\n  \tKopylov")
    ("hadronization-time", po::value<double>(&config.hadronizationTime)->default_value(0.), "Hadronization time, in fm/c. Particles emerging from collisions are forbidden to recollide within this time. Default: 0.")
    ;

  // Select options allowed on the command line
  cmdLineOptions.add(hiddenOptDesc).add(genericOptDesc).add(runOptDesc).add(physicsOptDesc);

  // Select options allowed in config files
  configFileOptions.add(runOptDesc).add(physicsOptDesc);

  // Select visible options
  visibleOptions.add(genericOptDesc).add(runOptDesc).add(physicsOptDesc);

  // Declare input-file as a positional option (if we just provide a file
  // name on the command line, it should be interpreted as an input-file
  // option).
  p.add("input-file", 1);

}

ConfigParser::~ConfigParser() {
}

G4INCL::Config *ConfigParser::parse(int argc, char *argv[]) {

  config.init();

  try {

    // Disable guessing of option names
    const int cmdstyle =
      po::command_line_style::default_style &
      ~po::command_line_style::allow_guessing;
    variablesMap.clear();

    // Result of the option processing
    po::store(po::command_line_parser(argc, argv).
              style(cmdstyle).
              options(cmdLineOptions).positional(p).run(), variablesMap);
    po::notify(variablesMap);

    // If an input file was specified, merge the options with the command-line
    // options.
    if(variablesMap.count("input-file")) {
      std::ifstream inputFileStream(config.inputFileName.c_str());
      if(!inputFileStream) {
        std::cerr << "Cannot open input file: " << config.inputFileName << '\n';
        return NULL;
      } else {
        // Merge options from the input file
        po::parsed_options parsedOptions = po::parse_config_file(inputFileStream, configFileOptions, true);

        // Make sure that the unhandled options are all "*-datafile-path"
        std::vector<std::string> unhandledOptions =
          po::collect_unrecognized(parsedOptions.options, po::exclude_positional);
        bool ignoreNext = false;
        const std::string match = "-datafile-path";
        for(std::vector<std::string>::const_iterator i=unhandledOptions.begin(), e=unhandledOptions.end(); i!=e; ++i) {
          if(ignoreNext) {
            ignoreNext=false;
            continue;
          }
          if(i->rfind(match) == i->length()-match.length()) {
            std::cout << "# Ignoring unrecognized option " << *i << '\n';
            ignoreNext = true;
          } else {
            std::cerr << "Error: unrecognized option " << *i << '\n';
            std::cerr << suggestHelpMsg;
            return NULL;
          }
        }

        // Store the option values in the variablesMap
        po::store(parsedOptions, variablesMap);
        po::notify(variablesMap);
      }
      inputFileStream.close();
    }

    // Process the options from the user-specific config file ~/.inclxxrc
    std::string configFileName;
    const char * const configFileVar = getenv("INCLXXRC");
    if(configFileVar)
      configFileName = configFileVar;
    else {
      const char * const homeDirectoryPointer = getenv("HOME");
      if(homeDirectoryPointer) { // Check if we can find the home directory
        std::string homeDirectory(homeDirectoryPointer);
        configFileName = homeDirectory + "/.inclxxrc";
      } else {
        std::cerr << "Could not determine the user's home directory. "
          << "Are you running Linux, Unix or BSD?"<< std::endl;
        return NULL;
      }
    }

    std::ifstream configFileStream(configFileName.c_str());
    if(configFileStream) {
      std::cout << "# Reading config file " << configFileName << std::endl;

      // Merge options from the input file
      po::parsed_options parsedOptions = po::parse_config_file(configFileStream, configFileOptions, true);
      po::store(parsedOptions, variablesMap);

      // Make sure that the unhandled options are all "*-datafile-path"
      std::vector<std::string> unhandledOptions =
        po::collect_unrecognized(parsedOptions.options, po::exclude_positional);
      bool ignoreNext = false;
      const std::string match = "-datafile-path";
      for(std::vector<std::string>::const_iterator i=unhandledOptions.begin(), e=unhandledOptions.end(); i!=e; ++i) {
        if(ignoreNext) {
          ignoreNext=false;
          continue;
        }
        if(i->rfind(match) == i->length()-match.length()) {
          std::cout << "Ignoring unrecognized option " << *i << std::endl;
          ignoreNext = true;
        } else {
          std::cerr << "Error: unrecognized option " << *i << std::endl;
          std::cerr << suggestHelpMsg;
          return NULL;
        }
      }

      // Store the option values in the variablesMap
      po::store(parsedOptions, variablesMap);
      po::notify(variablesMap);
    }
    configFileStream.close();



    /* *******************
     * Process the options
     * *******************/

    // -h/--help: print the help message and exit successfully
    if(variablesMap.count("help")) {
      std::cout
        << "Usage: INCLCascade [options] <input_file>" << std::endl
        << std::endl << "Options marked with a * are compulsory, i.e. they must be provided either on\nthe command line or in the input file." << std::endl
        << visibleOptions << std::endl;
      return NULL;
    }

    // --version: print the version string and exit successfully
    if(variablesMap.count("version")) {
      std::cout <<"INCL++ version " << config.getVersionString() << std::endl;
      return NULL;
    }

    // Check if the required options are present
    std::string missingOption("");
    if(!variablesMap.count("number-shots"))
      missingOption = "number-shots";
    else if(!variablesMap.count("target"))
      missingOption = "target";
    else if(!variablesMap.count("projectile"))
      missingOption = "projectile";
    else if(!variablesMap.count("energy"))
      missingOption = "energy";
    if(!missingOption.empty()) {
      std::cerr << "Required option " << missingOption << " is missing." << std::endl;
      std::cerr << suggestHelpMsg;
      return NULL;
    }

    // -p/--projectile: projectile species
    config.projectileSpecies = G4INCL::ParticleSpecies(config.projectileString);
    if(config.projectileSpecies.theType == G4INCL::UnknownParticle) {
      std::cerr << "Error: unrecognized particle type " << config.projectileString << std::endl;
      std::cerr << suggestHelpMsg;
      return NULL;
    }

    // -t/--target: target species
    if(variablesMap.count("target")) {
      config.targetSpecies = G4INCL::ParticleSpecies(config.targetString);
      if(config.targetSpecies.theType!=G4INCL::Composite) {
        std::cerr << "Unrecognized target. You specified: " << config.targetString << std::endl
          << "  The target nuclide must be specified in one of the following forms:" << std::endl
          << "    Fe56, 56Fe, Fe-56, 56-Fe, Fe_56, 56_Fe, Fe" << std::endl
          << "  You can also use IUPAC element names (such as Uuh)." << std::endl;
        std::cerr << suggestHelpMsg;
        return NULL;
      }
      if(config.targetSpecies.theA==0)
        config.naturalTarget = true;
    }

    // --pauli
    if(variablesMap.count("pauli")) {
      std::string pauliNorm = config.pauliString;
      std::transform(pauliNorm.begin(), pauliNorm.end(), pauliNorm.begin(), ::tolower);
      if(pauliNorm=="statistical")
        config.pauliType = G4INCL::StatisticalPauli;
      else if(pauliNorm=="strict")
        config.pauliType = G4INCL::StrictPauli;
      else if(pauliNorm=="strict-statistical")
        config.pauliType = G4INCL::StrictStatisticalPauli;
      else if(pauliNorm=="global")
        config.pauliType = G4INCL::GlobalPauli;
      else if(pauliNorm=="none")
        config.pauliType = G4INCL::NoPauli;
      else {
        std::cerr << "Unrecognized Pauli-blocking algorithm. Must be one of:" << std::endl
          << "  strict-statistical (default)" << std::endl
          << "  strict" << std::endl
          << "  statistical" << std::endl
          << "  global" << std::endl
          << "  none" << std::endl;
        std::cerr << suggestHelpMsg;
        return NULL;
      }
    }

    // --coulomb
    if(variablesMap.count("coulomb")) {
      std::string coulombNorm = config.coulombString;
      std::transform(coulombNorm.begin(), coulombNorm.end(), coulombNorm.begin(), ::tolower);
      if(coulombNorm=="non-relativistic")
        config.coulombType = G4INCL::NonRelativisticCoulomb;
      else if(coulombNorm=="none")
        config.coulombType = G4INCL::NoCoulomb;
      else {
        std::cerr << "Unrecognized Coulomb-distortion algorithm. Must be one of:" << std::endl
          << "  non-relativistic (default)" << std::endl
          << "  none" << std::endl;
        std::cerr << suggestHelpMsg;
        return NULL;
      }
    }

    // --potential
    if(variablesMap.count("potential")) {
      std::string potentialNorm = config.potentialString;
      std::transform(potentialNorm.begin(), potentialNorm.end(), potentialNorm.begin(), ::tolower);
      if(potentialNorm=="isospin-energy-smooth") {
        config.potentialType = G4INCL::IsospinEnergySmoothPotential;
      } else if(potentialNorm=="isospin-energy") {
        config.potentialType = G4INCL::IsospinEnergyPotential;
      } else if(potentialNorm=="isospin")
        config.potentialType = G4INCL::IsospinPotential;
      else if(potentialNorm=="constant")
        config.potentialType = G4INCL::ConstantPotential;
      else {
        std::cerr << "Unrecognized potential type. Must be one of:" << std::endl
          << "  isospin-energy-smooth" << std::endl
          << "  isospin-energy (default)" << std::endl
          << "  isospin" << std::endl
          << "  constant" << std::endl;
        std::cerr << suggestHelpMsg;
        return NULL;
      }
    }

    // --local-energy-BB
    if(variablesMap.count("local-energy-BB")) {
      std::string localEnergyBBNorm = config.localEnergyBBString;
      std::transform(localEnergyBBNorm.begin(), localEnergyBBNorm.end(), localEnergyBBNorm.begin(), ::tolower);
      if(localEnergyBBNorm=="always") {
        config.localEnergyBBType = G4INCL::AlwaysLocalEnergy;
      } else if(localEnergyBBNorm=="first-collision")
        config.localEnergyBBType = G4INCL::FirstCollisionLocalEnergy;
      else if(localEnergyBBNorm=="never")
        config.localEnergyBBType = G4INCL::NeverLocalEnergy;
      else {
        std::cerr << "Unrecognized local-energy-BB type. Must be one of:" << std::endl
          << "  always" << std::endl
          << "  first-collision (default)" << std::endl
          << "  never" << std::endl;
        std::cerr << suggestHelpMsg;
        return NULL;
      }
    }

    // --local-energy-pi
    if(variablesMap.count("local-energy-pi")) {
      std::string localEnergyPiNorm = config.localEnergyPiString;
      std::transform(localEnergyPiNorm.begin(), localEnergyPiNorm.end(), localEnergyPiNorm.begin(), ::tolower);
      if(localEnergyPiNorm=="always") {
        config.localEnergyPiType = G4INCL::AlwaysLocalEnergy;
      } else if(localEnergyPiNorm=="first-collision")
        config.localEnergyPiType = G4INCL::FirstCollisionLocalEnergy;
      else if(localEnergyPiNorm=="never")
        config.localEnergyPiType = G4INCL::NeverLocalEnergy;
      else {
        std::cerr << "Unrecognized local-energy-pi type. Must be one of:" << std::endl
          << "  always" << std::endl
          << "  first-collision" << std::endl
          << "  never (default)" << std::endl;
        std::cerr << suggestHelpMsg;
        return NULL;
      }
    }

    // -d/--de-excitation
    if(variablesMap.count("de-excitation")) {
      std::string deExcitationNorm = config.deExcitationString;
      std::transform(deExcitationNorm.begin(),
                     deExcitationNorm.end(),
                     deExcitationNorm.begin(), ::tolower);
      if(deExcitationNorm=="none")
        config.deExcitationType = G4INCL::DeExcitationNone;
#ifdef INCL_DEEXCITATION_ABLAXX
      else if(deExcitationNorm=="ablav3p")
        config.deExcitationType = G4INCL::DeExcitationABLAv3p;
#endif
#ifdef INCL_DEEXCITATION_ABLA07
      else if(deExcitationNorm=="abla07")
        config.deExcitationType = G4INCL::DeExcitationABLA07;
#endif
#ifdef INCL_DEEXCITATION_SMM
      else if(deExcitationNorm=="smm")
        config.deExcitationType = G4INCL::DeExcitationSMM;
#endif
#ifdef INCL_DEEXCITATION_GEMINIXX
      else if(deExcitationNorm=="geminixx")
        config.deExcitationType = G4INCL::DeExcitationGEMINIXX;
#endif
      else {
        std::cerr << "Unrecognized de-excitation model. "
          << "Must be one of:" << std::endl
          << deExcitationModelList << std::endl;
        std::cerr << suggestHelpMsg;
        return NULL;
      }
    }

#ifdef INCL_DEEXCITATION_FERMI_BREAKUP
    // --max-mass-fermi-breakup and --max-charge-fermi-breakup
    if(variablesMap.count("max-mass-fermi-breakup")) {
      if(config.maxMassFermiBreakUp<0 || config.maxMassFermiBreakUp>16) {
        std::cerr << "The maximum mass for Fermi breakup must belong to the [0,16] interval. " << std::endl;
        std::cerr << suggestHelpMsg;
        return NULL;
      }
    }
    if(variablesMap.count("max-charge-fermi-breakup")) {
      if(config.maxChargeFermiBreakUp<0 || config.maxChargeFermiBreakUp>8) {
        std::cerr << "The maximum charge for Fermi breakup must belong to the [0,8] interval. " << std::endl;
        std::cerr << suggestHelpMsg;
        return NULL;
      }
    }
#endif

    // --cluster-algorithm
    if(variablesMap.count("cluster-algorithm")) {
      std::string clusterAlgorithmNorm = config.clusterAlgorithmString;
      std::transform(clusterAlgorithmNorm.begin(),
                     clusterAlgorithmNorm.end(),
                     clusterAlgorithmNorm.begin(), ::tolower);
      if(clusterAlgorithmNorm=="none")
        config.clusterAlgorithmType = G4INCL::NoClusterAlgorithm;
      else if(clusterAlgorithmNorm=="intercomparison")
        config.clusterAlgorithmType = G4INCL::IntercomparisonClusterAlgorithm;
      else {
        std::cerr << "Unrecognized cluster algorithm. "
          << "Must be one of:" << std::endl
          << "  intercomparison (default)" << std::endl
          << "  none" << std::endl;
        std::cerr << suggestHelpMsg;
        return NULL;
      }
    }

    // --cluster-max-mass
    if(variablesMap.count("cluster-max-mass") && config.clusterMaxMass < 2 && config.clusterMaxMass > 12) {
      std::cerr << "Maximum cluster mass outside the allowed range. Must be between 2 and 12 (included)"
        << std::endl
        << suggestHelpMsg;
      return NULL;
    }

    // --separation-energies
    if(variablesMap.count("separation-energies")) {
      std::string separationEnergyNorm = config.separationEnergyString;
      std::transform(separationEnergyNorm.begin(),
                     separationEnergyNorm.end(),
                     separationEnergyNorm.begin(), ::tolower);
      if(separationEnergyNorm=="incl")
        config.separationEnergyType = G4INCL::INCLSeparationEnergy;
      else if(separationEnergyNorm=="real")
        config.separationEnergyType = G4INCL::RealSeparationEnergy;
      else if(separationEnergyNorm=="real-light")
        config.separationEnergyType = G4INCL::RealForLightSeparationEnergy;
      else {
        std::cerr << "Unrecognized separation-energies option. "
          << "Must be one of:" << std::endl
          << "  INCL (default)" << std::endl
          << "  real" << std::endl
          << "  real-light" << std::endl;
        std::cerr << suggestHelpMsg;
        return NULL;
      }
    } else {
      config.separationEnergyType = G4INCL::INCLSeparationEnergy;
    }

    // --fermi-momentum
    if(variablesMap.count("fermi-momentum")) {
      std::string fermiMomentumNorm = config.fermiMomentumString;
      std::transform(fermiMomentumNorm.begin(),
                     fermiMomentumNorm.end(),
                     fermiMomentumNorm.begin(), ::tolower);
      if(fermiMomentumNorm=="constant")
        config.fermiMomentumType = G4INCL::ConstantFermiMomentum;
      else if(fermiMomentumNorm=="constant-light")
        config.fermiMomentumType = G4INCL::ConstantLightFermiMomentum;
      else if(fermiMomentumNorm=="mass-dependent")
        config.fermiMomentumType = G4INCL::MassDependentFermiMomentum;
      else {
        // Try to convert the option value to a float, and bomb out on failure
        errno = 0;
        char *tail;
        config.fermiMomentum = strtod(fermiMomentumNorm.c_str(), &tail);
        if(errno || *tail!='\0') {
          std::cerr << "Unrecognized fermi-momentum option. "
            << "Must be one of:" << std::endl
            << "  constant (default)" << std::endl
            << "  constant-light" << std::endl
            << "  mass-dependent" << std::endl
            << "  [a postiive value]" << std::endl;
          std::cerr << suggestHelpMsg;
          return NULL;
        }
        if(config.fermiMomentum<=0.) {
          std::cerr << "Values passed to fermi-momentum must be positive." << std::endl;
          std::cerr << suggestHelpMsg;
          return NULL;
        }
        config.fermiMomentumType = G4INCL::ConstantFermiMomentum;
      }
    } else {
      config.fermiMomentumType = G4INCL::ConstantFermiMomentum;
    }

    // --rp-correlation / --rp-correlation-p / --rp-correlation-n
    if(variablesMap.count("rp-correlation")) {
      if(!variablesMap.count("rp-correlation-p") || variablesMap.find("rp-correlation-p")->second.defaulted())
        config.rpCorrelationCoefficientProton = config.rpCorrelationCoefficient;
      if(!variablesMap.count("rp-correlation-n") || variablesMap.find("rp-correlation-n")->second.defaulted())
        config.rpCorrelationCoefficientNeutron = config.rpCorrelationCoefficient;
    }

    // --cross-sections
    if(variablesMap.count("cross-sections")) {
      std::string crossSectionsNorm = config.crossSectionsString;
      std::transform(crossSectionsNorm.begin(), crossSectionsNorm.end(), crossSectionsNorm.begin(), ::tolower);
      if(crossSectionsNorm=="incl46")
        config.crossSectionsType = G4INCL::INCL46CrossSections;
      else if(crossSectionsNorm=="multipions")
        config.crossSectionsType = G4INCL::MultiPionsCrossSections;
      else if(crossSectionsNorm=="truncated-multipions")
        config.crossSectionsType = G4INCL::TruncatedMultiPionsCrossSections;
      else {
        std::cerr << "Unrecognized cross section parametrization. Must be one of:" << std::endl
          << "  multipions (default)" << std::endl
          << "  truncated-multipions" << std::endl
          << "  incl46" << std::endl;
        std::cerr << suggestHelpMsg;
        return NULL;
      }
    }

    // --max-number-multipions
    if(!variablesMap["max-number-multipions"].defaulted() && !variablesMap["cross-sections"].defaulted()) {
      if(config.crossSectionsType!=G4INCL::TruncatedMultiPionsCrossSections) {
        // enforce cross-sections=truncated-multipions
        config.crossSectionsString = "truncated-multipions";
        config.crossSectionsType = G4INCL::TruncatedMultiPionsCrossSections;
      }
    }

    // --phase-space-generator
    if(variablesMap.count("phase-space-generator")) {
      std::string phaseSpaceGeneratorNorm = config.phaseSpaceGenerator;
      std::transform(phaseSpaceGeneratorNorm.begin(),
                     phaseSpaceGeneratorNorm.end(),
                     phaseSpaceGeneratorNorm.begin(), ::tolower);
      if(phaseSpaceGeneratorNorm=="raubold-lynch")
        config.phaseSpaceGeneratorType = G4INCL::RauboldLynchType;
      else if(phaseSpaceGeneratorNorm=="kopylov")
        config.phaseSpaceGeneratorType = G4INCL::KopylovType;
      else {
        std::cerr << "Unrecognized phase-space-generator option. "
          << "Must be one of:" << std::endl
          << "  Raubold-Lynch (default)" << std::endl
          << "  Kopylov" << std::endl;
        std::cerr << suggestHelpMsg;
        return NULL;
      }
    } else {
      config.phaseSpaceGeneratorType = G4INCL::RauboldLynchType;
    }

    // --cascade-action
    if(variablesMap.count("cascade-action")) {
      std::string cascadeActionNorm = config.cascadeAction;
      std::transform(cascadeActionNorm.begin(),
                     cascadeActionNorm.end(),
                     cascadeActionNorm.begin(), ::tolower);
      if(cascadeActionNorm=="default")
        config.cascadeActionType = G4INCL::DefaultActionType;
      else if(cascadeActionNorm=="avatar-dump")
        config.cascadeActionType = G4INCL::AvatarDumpActionType;
      else {
        std::cerr << "Unrecognized cascade-action option. "
          << "Must be one of:" << std::endl
          << "  default (default)" << std::endl
          << "  avatar-dump" << std::endl;
        std::cerr << suggestHelpMsg;
        return NULL;
      }
    } else {
      config.cascadeActionType = G4INCL::DefaultActionType;
    }

    // -s/--suffix
    if(!variablesMap.count("suffix")) {
      // update the value in the variables_map
      variablesMap.insert(std::make_pair("suffix", po::variable_value(boost::any(config.fileSuffix), false)));
    }

    // --*-path: perform tilde expansion on the datafile paths
    if(variablesMap.count("inclxx-datafile-path"))
      config.INCLXXDataFilePath = G4INCL::String::expandPath(config.INCLXXDataFilePath);
#ifdef INCL_DEEXCITATION_ABLAXX
    if(variablesMap.count("ablav3p-cxx-datafile-path"))
      config.ablav3pCxxDataFilePath = G4INCL::String::expandPath(config.ablav3pCxxDataFilePath);
#endif
#ifdef INCL_DEEXCITATION_ABLA07
    if(variablesMap.count("abla07-datafile-path"))
      config.abla07DataFilePath = G4INCL::String::expandPath(config.abla07DataFilePath);
#endif
#ifdef INCL_DEEXCITATION_GEMINIXX
    if(variablesMap.count("geminixx-datafile-path"))
      config.geminixxDataFilePath = G4INCL::String::expandPath(config.geminixxDataFilePath);
#endif

    // --output: path expansion is applied
    if(variablesMap.count("output")) {
      config.outputFileRoot = G4INCL::String::expandPath(config.outputFileRoot);
    } else {
      // construct a reasonable output file root if not specified
      std::stringstream outputFileRootStream;
      // If an input file was specified, use its name as the output file root
      if(variablesMap.count("input-file"))
        outputFileRootStream << config.inputFileName << config.fileSuffix;
      else {
        outputFileRootStream.precision(0);
        outputFileRootStream.setf(std::ios::fixed, std::ios::floatfield);
        outputFileRootStream <<
          G4INCL::ParticleTable::getShortName(config.projectileSpecies) << "_" <<
          G4INCL::ParticleTable::getShortName(config.targetSpecies) << "_" <<
          config.projectileKineticEnergy;
        outputFileRootStream.precision(2);

        // Append suffixes to the output file root for each explicitly specified CLI option
        typedef po::variables_map::const_iterator BPOVMIter;
        for(BPOVMIter i=variablesMap.begin(), e=variablesMap.end(); i!=e; ++i) {
          std::string const &name = i->first;
          // Only process CLI options
          if(name!="projectile"
             && name!="target"
             && name!="energy"
             && name!="number-shots"
             && name!="random-seeds"
             && name!="random-number-generator"
             && name!="verbosity"
             && name!="verbose-event"
             && name!="suffix"
#ifdef INCL_ROOT_USE
             && name!="root-selection"
             && name!="concise-root-tree"
#endif
             && name!="inverse-kinematics"
             && name!="inclxx-datafile-path"
#ifdef INCL_DEEXCITATION_ABLA07
             && name!="abla07-datafile-path"
#endif
#ifdef INCL_DEEXCITATION_ABLAXX
             && name!="ablav3p-cxx-datafile-path"
#endif
#ifdef INCL_DEEXCITATION_GEMINIXX
             && name!="geminixx-datafile-path"
#endif
             ) {
               po::variable_value v = i->second;
               if(!v.defaulted()) {
                 const std::type_info &type = v.value().type();
                 if(type==typeid(std::string))
                   outputFileRootStream << "_" << name << "=" << v.as<std::string>();
                 else if(type==typeid(float))
                   outputFileRootStream << "_" << name << "=" << v.as<float>();
                 else if(type==typeid(double))
                   outputFileRootStream << "_" << name << "=" << v.as<double>();
                 else if(type==typeid(int))
                   outputFileRootStream << "_" << name << "=" << v.as<int>();
                 else if(type==typeid(bool))
                   outputFileRootStream << "_" << name << "=" << v.as<bool>();
               }
             }
        }

        outputFileRootStream << config.fileSuffix;
      }

      // update the variable
      config.outputFileRoot = outputFileRootStream.str();
      // update the value in the variables_map
      variablesMap.insert(std::make_pair("output", po::variable_value(boost::any(config.outputFileRoot), false)));
    }

    // -l/--logfile
    if(!variablesMap.count("logfile")) {
      // update the variable
      config.logFileName = config.outputFileRoot + ".log";
      // update the value in the variables_map
      variablesMap.insert(std::make_pair("logfile", po::variable_value(boost::any(config.logFileName), false)));
    } else {
      // perform path expansion
      config.logFileName = G4INCL::String::expandPath(config.logFileName);
    }

    // --random-number-generator
    if(variablesMap.count("random-number-generator")) {
      std::string randomNumberGeneratorNorm = config.randomNumberGenerator;
      std::transform(randomNumberGeneratorNorm.begin(),
                     randomNumberGeneratorNorm.end(),
                     randomNumberGeneratorNorm.begin(), ::tolower);
      if(randomNumberGeneratorNorm=="ranecu")
        config.rngType = G4INCL::RanecuType;
      else if(randomNumberGeneratorNorm=="ranecu3")
        config.rngType = G4INCL::Ranecu3Type;
      else {
        std::cerr << "Unrecognized random-number-generator option. "
          << "Must be one of:" << std::endl
          << "  Ranecu (2 seeds, default)" << std::endl
          << "  Ranecu3 (3 seeds)" << std::endl;
        std::cerr << suggestHelpMsg;
        return NULL;
      }
    } else {
      config.rngType = G4INCL::RanecuType;
    }

    // --random-seeds
    if(variablesMap.count("random-seeds")) {
      config.randomSeedVector.clear();

      std::vector<std::string> tokens = G4INCL::String::tokenize(config.randomSeeds, ", \t");
      for(std::vector<std::string>::const_iterator i=tokens.begin(), e=tokens.end(); i!=e; ++i) {
        if(!G4INCL::String::isInteger(*i)) {
          std::cerr << "Invalid random seed, must be an integer. Parsed token: "
            << *i << std::endl;
          std::cerr << suggestHelpMsg;
          return NULL;
        }

        std::stringstream ss(*i);
        int seed;
        ss >> seed;

        if(seed<randomSeedMin || seed>randomSeedMax) {
          std::cerr << "Invalid value for random-seed-1. "
            << "Allowed range: [" << randomSeedMin << ", " << randomSeedMax << "]." << std::endl;
          std::cerr << suggestHelpMsg;
          return NULL;
        }

        config.randomSeedVector.push_back(seed);
      }
    }

  }
  catch(std::exception& e)
  {
    std::cerr << e.what() << "\n";
    std::cerr << suggestHelpMsg;
    return NULL;
  }

  // Copy the config into a new object and return it
  G4INCL::Config *aConfig = new G4INCL::Config(config);
  return aConfig;
}

std::string ConfigParser::echo(G4INCL::Config const * const aConfig) {
  config = *aConfig;
  std::stringstream ss;
  ss << "###########################\n"
    << "### Start of input echo ###\n"
    << "###########################\n\n"
    << "# You may re-use this snippet of the log file as an input file!\n"
    << "# Options marked with a * are compulsory.\n"
    << "\n### Run options\n" << echoOptionsDescription(runOptDesc)
    << "\n### Physics options\n" << echoOptionsDescription(physicsOptDesc)
    << "\n# the projectile nuclide was parsed as Z=" << config.projectileSpecies.theZ
    << ", A=" << config.projectileSpecies.theA
    << "\n# the target nuclide was parsed as Z=" << config.targetSpecies.theZ;
  if(config.targetSpecies.theA>0)
    ss << ", A=" << config.targetSpecies.theA;
  else
    ss << ", natural target";
  ss << "\n\n#########################\n"
    << "### End of input echo ###\n"
    << "#########################" << std::endl;

  return ss.str();
}

std::string ConfigParser::echoOptionsDescription(const po::options_description &aDesc) {
  typedef std::vector< boost::shared_ptr< po::option_description > > OptVector;
  typedef std::vector< boost::shared_ptr< po::option_description > >::const_iterator OptIter;

  std::stringstream ss;
  ss << std::boolalpha;
  OptVector const &anOptVect = aDesc.options();
  for(OptIter opt=anOptVect.begin(), e=anOptVect.end(); opt!=e; ++opt) {
    std::string optDescription = (*opt)->description();
    G4INCL::String::wrap(optDescription);
    G4INCL::String::replaceAll(optDescription, "\n", "\n# ");
    ss << "\n# " << optDescription << '\n';
    const std::string &name = (*opt)->long_name();
    ss << name << " = ";
    po::variable_value const &value = variablesMap.find(name)->second;
    std::type_info const &type = value.value().type();
    if(type == typeid(std::string)) {
      const std::string s = value.as<std::string>();
      if(s.empty())
        ss << "\"\"";
      else
        ss << s;
    } else if(type == typeid(int))
      ss << value.as<int>();
    else if(type == typeid(unsigned int))
      ss << value.as<unsigned int>();
    else if(type == typeid(float))
      ss << value.as<float>();
    else if(type == typeid(double))
      ss << value.as<double>();
    else if(type == typeid(bool))
      ss << value.as<bool>();
    ss << '\n';
  }
  return ss.str();
}
#endif

#endif // __GENIE_INCL_ENABLED__
