#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ namespace genie;
#pragma link C++ namespace genie::utils;
#pragma link C++ namespace genie::utils::gsl;
#pragma link C++ namespace genie::utils::gsl::wrap;

#pragma link C++ class genie::QELXSec;
#pragma link C++ class genie::RESXSec;
#pragma link C++ class genie::DISXSec;
#pragma link C++ class genie::IMDXSec;
#pragma link C++ class genie::COHXSec;
#pragma link C++ class genie::NuElectronXSec;

// Wrappers for GSL/MathMore lib
#pragma link C++ class genie::utils::gsl::wrap::dXSec_dQ2_E;
#pragma link C++ class genie::utils::gsl::wrap::dXSec_dy_E;
#pragma link C++ class genie::utils::gsl::wrap::d2XSec_dxdy_E;
#pragma link C++ class genie::utils::gsl::wrap::d2XSec_dWdQ2_E;
#pragma link C++ class genie::utils::gsl::wrap::d2XSec_dxdy_Ex;
#pragma link C++ class genie::utils::gsl::wrap::d2XSec_dxdy_Ey;
#pragma link C++ class genie::utils::gsl::wrap::d2XSec_dWdQ2_EW;
#pragma link C++ class genie::utils::gsl::wrap::d2XSec_dWdQ2_EQ2;

// Depreciated wrappers
#pragma link C++ class genie::GXSecFunc;
#pragma link C++ class genie::Integrand_DXSec_DQ2_E;
#pragma link C++ class genie::Integrand_DXSec_Dy_E;
#pragma link C++ class genie::Integrand_D2XSec_DxDy_E;
#pragma link C++ class genie::Integrand_D2XSec_DWDQ2_E;
#pragma link C++ class genie::Integrand_D2XSec_DxDy_Ex;
#pragma link C++ class genie::Integrand_D2XSec_DxDy_Ey;
#pragma link C++ class genie::Integrand_D2XSec_DWDQ2_EW;
#pragma link C++ class genie::Integrand_D2XSec_DWDQ2_EQ2;

#endif
