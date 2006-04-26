//____________________________________________________________________________
/*!

\file     fragmentation_functions.h

\brief    Fragmentation functions

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  June 15, 2004

\cpright  Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _FRAGMENTATION_FUNCTIONS_H_
#define _FRAGMENTATION_FUNCTIONS_H_

namespace genie {

/*!
 \fn    double collins_spiller_fragmentation_function(double * x, double * par)
 \brief The Collins-Spiller fragmentation function
*/
  double collins_spiller_fragmentation_function  (double * x, double * par);

/*!
  \fn    double peterson_fragmentation_function(double * x, double * par)
  \brief The Peterson fragmentation function
*/

  double peterson_fragmentation_function         (double * x, double * par);

}      // genie namespace

#endif // _FRAGMENTATION_FUNCTIONS_H_
