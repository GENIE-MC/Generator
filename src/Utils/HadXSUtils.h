//____________________________________________________________________________
/*!

\namespace  genie::utils::hadxs

\brief      Simple functions and data for computing hadron interaction xsecs

\author     Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
            CCLRC, Rutherford Appleton Laboratory

\created    March 11, 2004

*/
//____________________________________________________________________________

#ifndef _HADXS_UTILS_H_
#define _HADXS_UTILS_H_

namespace genie {
namespace utils {

namespace hadxs
{
  double InelasticPionNucleonXSec (double Epion);
  double TotalPionNucleonXSec     (double Epion);

}      // hadxs namespace
}      // utils namespace
}      // genie namespace

#endif // _HADXS_UTILS_H_
