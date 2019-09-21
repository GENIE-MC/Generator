#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ namespace genie;
#pragma link C++ namespace genie::utils::frgmfunc;
#pragma link C++ namespace genie::utils::fragmrec;

#pragma link C++ class Pythia8Singleton;

#pragma link C++ class genie::Pythia6Hadronization;
#pragma link C++ class genie::Pythia8Hadronization;
#pragma link C++ class genie::KNOHadronization;
#pragma link C++ class genie::KNOPythiaHadronization;
#pragma link C++ class genie::CharmP6Hadronization;
#pragma link C++ class genie::CharmP8Hadronization;

#pragma link C++ class genie::FragmentationFunctionI;
#pragma link C++ class genie::PetersonFragm;
#pragma link C++ class genie::CollinsSpillerFragm;
#pragma link C++ function genie::utils::frgmfunc::peterson_func;
#pragma link C++ function genie::utils::frgmfunc::collins_spiller_func;

#endif
