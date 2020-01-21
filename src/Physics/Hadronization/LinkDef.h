#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ namespace genie;
#pragma link C++ namespace genie::utils::frgmfunc;
//#pragma link C++ namespace genie::utils::fragmrec;

#pragma link C++ class genie::PythiaBaseHadro2019;
#pragma link C++ class genie::Pythia6Hadro2019;
#pragma link C++ class genie::Pythia8Hadro2019;

#pragma link C++ class genie::AGKYLowW2019;
#pragma link C++ class genie::AGKY2019;
#pragma link C++ class genie::AGCharm2019;

#pragma link C++ class genie::FragmentationFunctionI;
#pragma link C++ class genie::PetersonFragm;
#pragma link C++ class genie::CollinsSpillerFragm;
#pragma link C++ function genie::utils::frgmfunc::peterson_func;
#pragma link C++ function genie::utils::frgmfunc::collins_spiller_func;

#endif
