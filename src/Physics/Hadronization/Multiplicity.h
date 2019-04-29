//____________________________________________________________________________
/*!

\class    genie::Multiplicity

\brief    An enumeration of hadronic final state multiplicities.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  November 26, 2004

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _MULTIPLICITY_H_
#define _MULTIPLICITY_H_

namespace genie {

  typedef enum EMultiplicity {

       kMltUndefined    = -1,
       kMltTotal,              
       kMltTotalCharged,      
       kMltTotalNegative,
       kMltTotalPositive,
       kMltTotalNeutral,
       kMltFwdTotal,          /* Fwd: in Forward Hemisphere: xF>0 */
       kMltFwdCharged,
       kMltFwdNegative,
       kMltFwdPositive,
       kMltFwdNeutral,
       kMltBkwTotal,          /* Bkw: in Backward Hemisphere xF<0 */
       kMltBkwCharged,
       kMltBkwNegative,
       kMltBkwPositive,
       kMltBkwNeutral

  } Multiplicity_t;

}        // genie namespace

#endif   // _MULTIPLICITY_H_
