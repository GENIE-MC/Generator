//_____________________________________________________________________________
/*!

\class    genie::nuvld::facades::NGNucleus

\brief    NeuGEN's Nucleus enumeration

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>
          Hugh Gallagher      (Tufts University) <gallag@minos.phy.tufts.edu>

\created  August 2004
*/
//_____________________________________________________________________________

#ifndef _NUCLEUS_H_
#define _NUCLEUS_H_

#ifndef ROOT_Rtypes
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "Rtypes.h"
#endif
#endif

namespace genie   {
namespace nuvld   {
namespace facades {

typedef enum ENGNucleus {

  e_free              = 0,
  e_C12               = 274,
  e_O16               = 284,
  e_Fe56              = 372,
  e_undefined_nucleus

} NGNucleus_t;

class NGNucleus {

  public:

     static const char * AsString(NGNucleus_t nucleus) {
       
       switch(nucleus) {
         case e_free:  return "Free nucleon";       break;
         case e_C12:   return "Carbon";       break;
         case e_O16:   return "Oxygen";       break;
         case e_Fe56:  return "Iron";   break;

         case e_undefined_nucleus:
         default:            
                      return "unknown Nucleus"; break;
       }
     }
     
ClassDef(NGNucleus, 0)
};

} // facades namespace
} // nuvld   namespace
} // genie   namespace

#endif

