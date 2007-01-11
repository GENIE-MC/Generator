//_____________________________________________________________________________
/*!

\class    genie::nuvld::facades::NGFlavor

\brief    NeuGEN's Flavor enumeration

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>
          Hugh Gallagher      (Tufts University) <gallag@minos.phy.tufts.edu>

\created  August 2004
*/
//_____________________________________________________________________________

#ifndef _FLAVOR_H_
#define _FLAVOR_H_

#ifndef ROOT_Rtypes
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "Rtypes.h"
#endif
#endif

namespace genie   {
namespace nuvld   {
namespace facades {

typedef enum ENGFlavor {

  e_e = 1,
  e_mu,
  e_tau,
  e_undefined_flavor

} NGFlavor_t;

class NGFlavor {

  public:

     virtual ~NGFlavor() { }

     static const char * AsString(NGFlavor_t flavor) 
     {
       switch(flavor) {
         case e_e:                return "Electron Flavor "; break;
         case e_mu:               return "Muon Flavor";      break;
         case e_tau:              return "Tau Flavor";       break;
         case e_undefined_flavor:
         default:            
              return "Unknown Flavor";   break;
       }
       return "Unknown Flavor"; 
     }

     static NGFlavor_t GetFromCode(int pdgc) 
     {
       if      (pdgc == 5 || pdgc ==  6) return e_e;
       else if (pdgc == 7 || pdgc ==  8) return e_mu;
       else if (pdgc == 9 || pdgc == 10) return e_tau;
       else                              return e_undefined_flavor;
     }

ClassDef(NGFlavor, 0)
};

} // facades namespace
} // nuvld   namespace
} // genie   namespace

#endif

