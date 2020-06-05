#include "Framework/Conventions/GBuild.h"
#ifdef __GENIE_INCL_ENABLED__

/** \file INCLConfigParser.h
 * \brief A class for parsing input and producing a INCL++ Config
 *        Copied from INC++ main
 *  class name can't be changed from "ConfigParser" to "INCLConfigParser"
 *  due to code in $INCL_SRC_DIR/utils/include/G4INCLConfig.hh:
 *       friend class ::ConfigParser;
 *  and this code needs the access to the private data members
 *  By default copy of the class isn't compiled into any available INC++ library
 *  so we're stuck with the name for the class, but we'll leave the
 *  filename changed to indicate that it references INCL++
 *
 * \date 17th July 2014
 * \author Davide Mancusi
 */

#ifndef INCLCONFIGPARSER_HH_
#define INCLCONFIGPARSER_HH_

#ifdef HAS_BOOST_PROGRAM_OPTIONS

#include <string>
#include "G4INCLConfig.hh"

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/positional_options.hpp>
#include <boost/program_options/variables_map.hpp>

namespace po = boost::program_options;

class ConfigParser {
  public:
     ConfigParser();
    ~ConfigParser();

    G4INCL::Config *parse(int argc, char *argv[]);
    std::string echo(G4INCL::Config const * const aConfig);

  private:
    G4INCL::Config config;
    po::options_description runOptDesc;
    po::options_description hiddenOptDesc;
    po::options_description genericOptDesc;
    po::options_description physicsOptDesc;
    po::options_description cmdLineOptions;
    po::options_description configFileOptions;
    po::options_description visibleOptions;
    po::positional_options_description p;
    po::variables_map variablesMap;

    static const int randomSeedMin = 1;
    static const int randomSeedMax = ((1<<30)-1)+(1<<30);

    static const std::string suggestHelpMsg;

    // Define the names of the de-excitation models
    static const std::string theNoneName;
#ifdef INCL_DEEXCITATION_SMM
    static const std::string theSMMName;
#endif
#ifdef INCL_DEEXCITATION_GEMINIXX
    static const std::string theGEMINIXXName;
#endif
#ifdef INCL_DEEXCITATION_ABLAXX
    static const std::string theABLAv3pName;
#endif
#ifdef INCL_DEEXCITATION_ABLA07
    static const std::string theABLA07Name;
#endif

    static const std::string listSeparator;

    std::string deExcitationModelList;
    std::string defaultDeExcitationModel;

    std::string echoOptionsDescription(const po::options_description &aDesc);
};

// alias for GENIE
typedef ConfigParser INCLConfigParser;

#endif // HAS_BOOST_PROGRAM_OPTIONS

#endif // INCLCONFIGPARSER_HH_

#endif // __GENIE_INCL_ENABLED__
