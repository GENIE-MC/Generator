//_____________________________________________________________________________
/*!

\class    genie::nuvld::facades::NGKineVar

\brief    NeuGEN's Kinematic Variables enumeration

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>
          Hugh Gallagher      (Tufts University) <gallag@minos.phy.tufts.edu>

\created  August 2004
*/
//_____________________________________________________________________________

#ifndef _KINEMATIC_VARIABLE_H_
#define _KINEMATIC_VARIABLE_H_

#ifndef ROOT_Rtypes
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "Rtypes.h"
#endif
#endif

namespace genie   {
namespace nuvld   {
namespace facades {

typedef enum ENGKineVar {

  e_nokv = 0,
  e_qqs = 1,
  e_w,
  e_x,
  e_y,
  e_logqqs,
  e_costh,    /* cosine of scattering angle in lab */
  e_ep,       /* scattered electron energy */
  e_undefined_kinematic_variable

} NGKineVar_t;

class NGKineVar {

  public:

     virtual ~NGKineVar() { }

     static const char * AsString(NGKineVar_t kid) 
     {
       switch(kid) {
         case e_qqs:      return "qqs";        break;
         case e_w:        return "W";          break;
         case e_x:        return "x";          break;
         case e_y:        return "y";          break;
         case e_logqqs:   return "log(qqs)";   break;
         case e_costh:    return "cos(theta)"; break;
         case e_ep:       return "Ep";         break;

         case e_undefined_kinematic_variable:
         default:            
            return "Unknown kinematic variable"; break;
       }
       return "Unknown kinematic variable";
     }
     
     static NGKineVar_t GetKineVarFromCode(int code)
     {
        if      (code == 0) return e_nokv;
        else if (code == 1) return e_qqs;
        else if (code == 2) return e_w;
        else if (code == 3) return e_x;
        else if (code == 4) return e_y;
        else if (code == 5) return e_logqqs;
        else if (code == 6) return e_costh;
        else if (code == 7) return e_ep;
        else                return e_undefined_kinematic_variable;
     }

ClassDef(NGKineVar, 0)
};

} // facades namespace
} // nuvld   namespace
} // genie   namespace

#endif

