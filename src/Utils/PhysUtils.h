//____________________________________________________________________________
/*!

\namespace  genie::utils::phys

\brief      Various physics formulas

\author     Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
            STFC, Rutherford Appleton Laboratory

\created    January 22, 2008

\cpright    Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
            For the full text of the license visit http://copyright.genie-mc.org
            or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _PHYS_UTILS_H_
#define _PHYS_UTILS_H_

namespace genie {
namespace utils {

namespace phys
{
  //
  // parameterizations of the longitudinal to transverse cross section ratio R 
  //
  double R99118   (double x, double Q2); ///< PRL 98, 142301, 2007
  double RWhitlow (double x, double Q2);

} // phys  namespace
} // utils namespace
} // genie namespace

#endif // _PHYS_UTILS_H_
