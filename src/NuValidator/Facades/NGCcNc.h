//_____________________________________________________________________________
/*!

\class    genie::nuvld::facades::NGCcNc

\brief    NeuGEN's CC/NC enumeration

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>
          Hugh Gallagher      (Tufts University) <gallag@minos.phy.tufts.edu>

\created  August 2004
*/
//_____________________________________________________________________________

#ifndef _CC_NC_H_
#define _CC_NC_H_

#ifndef ROOT_Rtypes
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "Rtypes.h"
#endif
#endif

namespace genie   {
namespace nuvld   {
namespace facades {

typedef enum ENGCcNc {

  e_cc = 1,
  e_nc,
  e_both_ccnc,
  e_undefined_ccnc

} NGCcNc_t;

class NGCcNc {

  public:

     static const char * AsString(NGCcNc_t ccnc) 
     {
       switch(ccnc) {
         case e_cc:          return "CC";     break;
         case e_nc:          return "NC";     break;
         case e_both_ccnc:   return "CC+NC";  break;

         case e_undefined_ccnc:
         default:            
                      return "unknown CCNC"; break;
       }
     }

     static NGCcNc_t GetFromCode(int code) 
     {
        if      (code == 1) return e_cc;
        else if (code == 2) return e_nc;
        else if (code == 3) return e_both_ccnc;
        else                return e_undefined_ccnc;
     }
     
ClassDef(NGCcNc, 0)
};

} // facades namespace
} // nuvld   namespace
} // genie   namespace

#endif

