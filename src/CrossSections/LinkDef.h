#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ namespace genie;
#pragma link C++ namespace genie::utils;
#pragma link C++ namespace genie::utils::gsl;

#pragma link C++ class genie::QELXSec;
#pragma link C++ class genie::RESXSec;
#pragma link C++ class genie::DISXSec;
#pragma link C++ class genie::IMDXSec;
#pragma link C++ class genie::COHXSec;
#pragma link C++ class genie::NuElectronXSec;

// Wrappers for GSL/MathMore lib
#pragma link C++ class genie::utils::gsl::dXSec_dQ2_E;
#pragma link C++ class genie::utils::gsl::dXSec_dy_E;
#pragma link C++ class genie::utils::gsl::d2XSec_dxdy_E;
#pragma link C++ class genie::utils::gsl::d2XSec_dWdQ2_E;
#pragma link C++ class genie::utils::gsl::d2XSec_dxdy_Ex;
#pragma link C++ class genie::utils::gsl::d2XSec_dxdy_Ey;
#pragma link C++ class genie::utils::gsl::d2XSec_dWdQ2_EW;
#pragma link C++ class genie::utils::gsl::d2XSec_dWdQ2_EQ2;
#pragma link C++ class genie::utils::gsl::d5XSecAR;
#pragma link C++ class genie::utils::gsl::d5Xsec_dEldOmegaldOmegapi;
#pragma link C++ class genie::utils::gsl::d4Xsec_dEldThetaldOmegapi;
#pragma link C++ class genie::utils::gsl::d3Xsec_dOmegaldThetapi;
#pragma link C++ class genie::utils::gsl::dXSec_dElep_AR;
#pragma link C++ class genie::utils::gsl::dXSec_Log_Wrapper;

#endif
