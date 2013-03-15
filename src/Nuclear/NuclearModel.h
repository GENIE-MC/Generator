//____________________________________________________________________________
/*!

\class    genie::NuclearModel

\brief    Encapsulates an enumeration of nuclear model types

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  Jun 09, 2009

\cpright  Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE 
*/
//____________________________________________________________________________

#ifndef _NUCL_MODEL_TYPE_H_
#define _NUCL_MODEL_TYPE_H_

#ifndef ROOT_Rtypes
#include "Rtypes.h"
#endif

namespace genie {

typedef enum ENuclearModel {

   kNucmUndefined    = -1, 
   kNucmFermiGas,
   kNucmSpectralFunc

} NuclearModel_t; 
  

class NuclearModel {

public:
  static const char * AsString(NuclearModel_t nucmod) {
     switch (nucmod) {
     case kNucmUndefined:     return "Undefined nuclear model";  break;
     case kNucmFermiGas:      return "Fermi gas model";          break;
     case kNucmSpectralFunc:  return "Spectral function model";  break;
     default:                 break;
     }
     return " ";
  }

};

}
#endif
