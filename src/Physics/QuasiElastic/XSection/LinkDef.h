#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ namespace genie;

#pragma link C++ class genie::ELFormFactorsModelI;
#pragma link C++ class genie::ELFormFactors;

#pragma link C++ class genie::BBA03ELFormFactorsModel;
#pragma link C++ class genie::BBA05ELFormFactorsModel;
#pragma link C++ class genie::BBA07ELFormFactorsModel;
#pragma link C++ class genie::DipoleELFormFactorsModel;
#pragma link C++ class genie::KellyELFormFactorsModel;
#pragma link C++ class genie::TransverseEnhancementFFModel;

#pragma link C++ class genie::AhrensNCELPXSec;
#pragma link C++ class genie::RosenbluthPXSec;

#pragma link C++ class genie::LwlynSmithQELCCPXSec;
#pragma link C++ class genie::LwlynSmithFFCC;
#pragma link C++ class genie::LwlynSmithFFNC;
#pragma link C++ class genie::LwlynSmithFF;
#pragma link C++ class genie::LwlynSmithFFEM;
#pragma link C++ class genie::LwlynSmithFFDeltaS;
#pragma link C++ class genie::AxialFormFactorModelI;
#pragma link C++ class genie::AxialFormFactor;
#pragma link C++ class genie::DipoleAxialFormFactorModel;
#pragma link C++ class genie::ZExpAxialFormFactorModel;
#pragma link C++ class genie::KuzminNaumov2016AxialFormFactorModel;
#pragma link C++ class genie::NievesQELCCPXSec;
#pragma link C++ class genie::SuSAv2QELPXSec;
#pragma link C++ class genie::SmithMonizQELCCPXSec;
#pragma link C++ class genie::SmithMonizQELCCXSec;
#pragma link C++ class genie::SmithMonizUtils;

#pragma link C++ class genie::QELFormFactors;
#pragma link C++ class genie::QELFormFactorsModelI;

#pragma link C++ class genie::QELXSec;
#pragma link C++ class genie::NewQELXSec;

#pragma link C++ class genie::Rank2LorentzTensorI;
#pragma link C++ class genie::LeptonTensor;
#pragma link C++ class genie::ManualResponseTensor;

// Wrappers for GSL/MathMore lib
#pragma link C++ class genie::utils::gsl::d2Xsec_dQ2dv;

#endif
