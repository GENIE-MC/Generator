#ifdef __CINT__
#include "Framework/Conventions/GBuild.h"

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ namespace genie;
//#pragma link C++ namespace genie::utils::mec;

#pragma link C++ class genie::MECInteractionListGenerator;
#pragma link C++ class genie::MECGenerator;
#ifdef __GENIE_INCL_ENABLED__
#pragma link C++ class genie::MECGeneratorINCL;
#endif

#endif
