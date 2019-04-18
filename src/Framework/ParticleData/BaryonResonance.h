//____________________________________________________________________________
/*!

\class    genie::BaryonResonance

\brief    An enumeration of Baryon Resonances more frequently used in 
          resonance neutrino-nucleon/nucleus models.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  May 06, 2004

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _BARYON_RESONANCE_H_
#define _BARYON_RESONANCE_H_

namespace genie {

  typedef enum EResonance {
       kNoResonance = -1,
       kP33_1232    =  0,
       kS11_1535    =  1,
       kD13_1520    =  2,
       kS11_1650    =  3,
       kD13_1700    =  4,
       kD15_1675    =  5,
       kS31_1620    =  6,
       kD33_1700    =  7,
       kP11_1440    =  8,
       kP33_1600    =  9,
       kP13_1720    = 10,
       kF15_1680    = 11,
       kP31_1910    = 12,
       kP33_1920    = 13,
       kF35_1905    = 14,
       kF37_1950    = 15,
       kP11_1710    = 16,
       kF17_1970    = 17
  } Resonance_t;

}        // genie namespace

#endif   // _BARYON_RESONANCE_H_
