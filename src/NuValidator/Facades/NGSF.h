//_____________________________________________________________________________
/*!

\class    genie::nuvld::facades::NGSF

\brief    NeuGEN's structure function enumeration

\author   Costas Andreopoulos (Rutherford Lab.)  <costas.andreopoulos \at stfc.ac.uk>
          Hugh Gallagher      (Tufts University) <gallag@minos.phy.tufts.edu>

\created  August 2004
*/
//_____________________________________________________________________________

#ifndef _NGSF_H_
#define _NGSF_H_

#ifndef ROOT_Rtypes
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "Rtypes.h"
#endif
#endif

namespace genie   {
namespace nuvld   {
namespace facades {

typedef enum ENGSF {

  e_F1 = 0,
  e_F2,
  e_F3,
  e_F4,
  e_F5,
  e_F6,
  e_xF3,
  e_UndefinedSF

} NGSF_t;

class NGSF {

  public:

     virtual ~NGSF() { }

     static const char * AsString(NGSF_t sf) 
     {  
       switch(sf) {
         case e_F1:   return "F1";   break;
         case e_F2:   return "F2";   break;
         case e_F3:   return "F3";   break;
         case e_F4:   return "F4";   break;
         case e_F5:  return  "F5";   break;
         case e_F6:  return  "F6";   break;
         case e_xF3:  return "xF3";  break;

         case e_UndefinedSF:
         default:            
             return "Unknown structure function"; break;
       }
       return "Unknown structure function"; 
     }
     
     static NGSF_t GetSFFromCode(int code)
     {
        if      (code ==  0) return e_F1;
        else if (code ==  1) return e_F2;
        else if (code ==  2) return e_F3;
        else if (code ==  3) return e_F4;
        else if (code ==  4) return e_F5;
        else if (code ==  5) return e_F6;
        else if (code ==  6) return e_xF3;
        else                 return e_UndefinedSF;
     }

ClassDef(NGSF, 0)
};

} // facades namespace
} // nuvld   namespace
} // genie   namespace


#endif

