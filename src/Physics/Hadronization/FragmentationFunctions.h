//____________________________________________________________________________
/*!

\namespace genie::utils::frgmfunc

\brief     Fragmentation functions

\author    Costas Andreopoulos <c.andreopoulos \at cern.ch>
           University of Liverpool

\created   June 15, 2004

\cpright   Copyright (c) 2003-2025, The GENIE Collaboration
           For the full text of the license visit http://copyright.genie-mc.org           
*/
//____________________________________________________________________________

#ifndef _FRAGMENTATION_FUNCTIONS_H_
#define _FRAGMENTATION_FUNCTIONS_H_

namespace genie    {
namespace utils    {
namespace frgmfunc {

/*!
 \fn    double collins_spiller_func(double * x, double * par)
 \brief The Collins-Spiller fragmentation function
*/
  double collins_spiller_func(double * x, double * par);

/*!
  \fn    double peterson_func(double * x, double * par)
  \brief The Peterson fragmentation function
*/
  double peterson_func(double * x, double * par);

} // frgmfunc namespace
} // utils    namespace
} // genie    namespace

#endif // _FRAGMENTATION_FUNCTIONS_H_
