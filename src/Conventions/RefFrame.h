//____________________________________________________________________________
/*!

\class    genie::RefFrame

\brief    An enumeration of reference frames.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 04, 2004
 
*/
//____________________________________________________________________________

#ifndef _REF_FRAME_H_
#define _REF_FRAME_H_

namespace genie {

typedef enum ERefFrame {

  kNullRefFrame = 0,
  kRfCenterOfMass,
  kRfTargetAtRest,
  kRfStruckNucAtRest,
  kRfLab  

} RefFrame_t;

class RefFrame {

public:

  static char * AsString(RefFrame_t rf) 
  {
    switch (rf) {
       case (kRfCenterOfMass) :     return "CENTER-OF-MASS-FRAME";      break;
       case (kRfTargetAtRest) :     return "TARGET-REST-FRAME";         break;
       case (kRfStruckNucAtRest) :  return "STRUCK-NUCLEON-REST-FRAME"; break;
       case (kRfLab) :              return "LAB-FRAME";                 break;       
       case (kNullRefFrame) :
       default :  
                return "Undefined REF-FRAME";
    }
  }

};

}        // genie namespace

#endif   // _REF_FRAME_H_
