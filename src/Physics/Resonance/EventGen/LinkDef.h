#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ namespace genie;

#pragma link C++ class genie::RESPrimaryLeptonGenerator;
#pragma link C++ class genie::RESKinematicsGenerator;
#pragma link C++ class genie::RESInteractionListGenerator;

#pragma link C++ class genie::RSPPHadronicSystemGenerator;
#pragma link C++ class genie::SPPEventGenerator;
#pragma link C++ class genie::RESHadronicSystemGenerator;
#pragma link C++ class genie::RSPPResonanceSelector;
#pragma link C++ class genie::RSPPInteractionListGenerator;

// Wrappers for GSL/MathMore lib
#pragma link C++ class genie::utils::gsl::d4XSecMK_dWQ2CosThetaPhi_E;

#endif
