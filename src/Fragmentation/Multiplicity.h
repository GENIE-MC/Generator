//____________________________________________________________________________
/*!

\class    genie::Multiplicity

\brief    An enumeration of hadronic final state multiplicities.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  November 26, 2004

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
