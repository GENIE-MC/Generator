//____________________________________________________________________________
/*!

\file     breit_wigner_func

\brief    Breit Wigner functions

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  November 22, 2004
 
*/
//____________________________________________________________________________

#ifndef _BREIT_WIGNER_FUNCTIONS_H_
#define _BREIT_WIGNER_FUNCTIONS_H_

#include "Config/Registry.h"

namespace genie {

  //-- A realistic Breit-Wigner distribution with L-dependent width.

  double breit_wigner_L (double W, const Registry & res_config);

  //-- A simple Breit-Wigner distribution.

  double breit_wigner   (double W, const Registry & res_config);

}        // genie namespace

#endif   // _BREIT_WIGNER_FUNCTIONS_H_
