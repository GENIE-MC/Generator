//____________________________________________________________________________
/*!

\namespace  genie::xsec_utils

\brief      Simple non-neutrino cross section functions

\author     Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
            CCLRC, Rutherford Appleton Laboratory

\created    March 11, 2004

*/
//____________________________________________________________________________

#ifndef _XSEC_UTILS_H_
#define _XSEC_UTILS_H_

namespace genie {

namespace xsec_utils
{
  double InelasticPionNucleonXSec (double Epion);
  double TotalPionNucleonXSec     (double Epion);

}      // xsec_utils namespace
}      // genie namespace

#endif // _XSEC_UTILS_H_
