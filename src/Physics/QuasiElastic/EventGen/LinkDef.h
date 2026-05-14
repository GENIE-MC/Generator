#ifdef __CINT__
#include "Framework/Conventions/GBuild.h"

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ namespace genie;

#pragma link C++ class genie::QELPrimaryLeptonGenerator;
#pragma link C++ class genie::QELHadronicSystemGenerator;
#pragma link C++ class genie::QELInteractionListGenerator;
#pragma link C++ class genie::QELKinematicsGenerator;
#pragma link C++ class genie::QELEventGeneratorSuSA;
#pragma link C++ class genie::QELEventGenerator;
#pragma link C++ class genie::QELEventGeneratorSM;
#ifdef __GENIE_INCL_ENABLED__
#pragma link C++ class genie::QELEventGeneratorINCL;
#endif

#endif
