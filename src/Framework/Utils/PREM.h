//____________________________________________________________________________
/*!

\namespace  genie::utils::prem

\brief      Preliminary Earth Model

\author     Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
            University of Liverpool & STFC Rutherford Appleton Lab

\created    August 07, 2009

\cpright    Copyright (c) 2003-2019, The GENIE Collaboration
            For the full text of the license visit http://copyright.genie-mc.org
            or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _PREM_H_
#define _PREM_H_

namespace genie {
namespace utils {

namespace prem
{
  //
  // the earth density profile as given by the Preliminary Earth Model
  // by Adam Dziewonski, Earth Structure, Global,
  // in The Encyclopedia of Solid Earth Geophysics, David E. James ed.
  // (Van Nostrand Reinhold, New York, 1989) p 331.
  //
  double Density(double r);

} // prem  namespace
} // utils namespace
} // genie namespace

#endif // _PREM_H_
