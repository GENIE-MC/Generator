//____________________________________________________________________________
/*!

\struct   genie::PDF_t

\brief    A struct to hold PDF set data

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 04, 2004

*/ 
//____________________________________________________________________________

#ifndef _PDF_T_H_
#define _PDF_T_H_

namespace genie {

typedef struct EPDF {

  double uval;
  double dval;
  double usea;
  double dsea;
  double str; 
  double chm; 
  double bot;
  double top; 
  double gl;

} PDF_t;

}        // genie namespace

#endif   // _PDF_T_H_
