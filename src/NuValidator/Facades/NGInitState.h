//_____________________________________________________________________________
/*!

\class    genie::nuvld::facades::NGInitState

\brief    NeuGEN's Initial State information

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>
          Hugh Gallagher      (Tufts University) <gallag@minos.phy.tufts.edu>

\created  August 2004
*/
//_____________________________________________________________________________

#ifndef _INIT_STATE_H_
#define _INIT_STATE_H_

#ifndef ROOT_Rtypes
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "Rtypes.h"
#endif
#endif

namespace genie   {
namespace nuvld   {
namespace facades {

typedef enum ENGInitState {

  e_vp = 1,
  e_vn,
  e_vbp,
  e_vbn,
  e_vN,
  e_vbN,
  e_vA,
  e_vbA,
  e_lp,
  e_ln,
  e_lN,
  e_undefined_init_state

} NGInitState_t;

class NGInitState {

  public:

     virtual ~NGInitState(){}

     static const char * AsString(NGInitState_t initial_state) 
     {  
       switch(initial_state) {
         case e_vp:   return "v + p";       break;
         case e_vn:   return "v + n";       break;
         case e_vN:   return "v + N";       break;
         case e_vA:   return "v + A";       break;
         case e_vbp:  return "v_bar + p";   break;
         case e_vbn:  return "v_bar + n";   break;
         case e_vbN:  return "v_bar + N";   break;
         case e_vbA:  return "v_bar + A";   break;
         case e_lp:   return "l + p";       break;
         case e_ln:   return "l + n";       break;
         case e_lN:   return "l + N";       break;

         case e_undefined_init_state:
         default:            
            return "Unknown initial state"; break;
       }
       return "Unknown initial state"; 
     }
     
     static NGInitState_t GetInitStateFromCode(int code)
     {
        if      (code ==  1) return e_vp;
        else if (code ==  2) return e_vn;
        else if (code ==  5) return e_vN;
        else if (code ==  7) return e_vA;
        else if (code ==  3) return e_vbp;
        else if (code ==  4) return e_vbn;
        else if (code ==  6) return e_vbN;
        else if (code ==  8) return e_vbA;
        else if (code ==  9) return e_lp;
        else if (code == 10) return e_ln;
        else if (code == 11) return e_lN;
        else                 return e_undefined_init_state;
     }

ClassDef(NGInitState, 0)
};

} // facades namespace
} // nuvld   namespace
} // genie   namespace


#endif

