#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ namespace genie;

#pragma link C++ class genie::FragmentationFunctionI;
#pragma link C++ class genie::PetersonFragm;
#pragma link C++ class genie::CollinsSpillerFragm;

#pragma link C++ class genie::HadronizationModelI;
#pragma link C++ class genie::PythiaHadronization;
#pragma link C++ class genie::KNOHadronization;
#pragma link C++ class genie::KNOHadronization2;

#pragma link C++ class genie::MultiplicityProbModelI;
#pragma link C++ class genie::SchmitzMultiplicityModel;
//#pragma link C++ class genie::MultiplicityProb;

#pragma link C++ class genie::CharmFractionTablePool;
#pragma link C++ class genie::CharmFractionTable;

#pragma link C++ function genie::peterson_fragmentation_function;
#pragma link C++ function genie::collins_spiller_fragmentation_function;

#endif
