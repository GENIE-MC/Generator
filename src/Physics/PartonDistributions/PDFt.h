//____________________________________________________________________________
/*!

\struct   genie::PDF_t

\brief    A struct to hold PDF set data

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  May 04, 2004

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
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
