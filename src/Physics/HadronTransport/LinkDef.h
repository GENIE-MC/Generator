#ifdef __CINT__
#include "Framework/Conventions/GBuild.h"

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ namespace genie;
#pragma link C++ namespace genie::utils::intranuke;

#pragma link C++ class genie::INukeHadroData;
#pragma link C++ class genie::INukeHadroData2018;
#pragma link C++ class genie::INukeDeltaPropg;
//#pragma link C++ class genie::INukePhotoPropg;

#pragma link C++ class genie::HadronTransporter;
#pragma link C++ class genie::NucBindEnergyAggregator;

#pragma link C++ class genie::Intranuke;
#pragma link C++ class genie::HAIntranuke;
#pragma link C++ class genie::Intranuke2018;
#pragma link C++ class genie::HAIntranuke2018;
#pragma link C++ class genie::HNIntranuke2018;

#ifdef __GENIE_INCL_ENABLED__
#pragma link C++ class genie::HINCLCascadeIntranuke;
#endif

#ifdef __GENIE_GEANT4_INTERFACE_ENABLED__
#pragma link C++ class genie::HG4BertCascIntranuke;
#endif

#endif
