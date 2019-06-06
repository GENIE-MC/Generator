//____________________________________________________________________________
/*!

\class    genie::NuclearModel

\brief    Encapsulates an enumeration of nuclear model types

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  Jun 09, 2009

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
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
   kNucmLocalFermiGas,
   kNucmSpectralFunc,
   kNucmEffSpectralFunc

} NuclearModel_t; 

typedef enum EFermiMoverInteractionType {
  kFermiMoveDefault = 0,
  kFermiMoveEffectiveSF1p1h,
  kFermiMoveEffectiveSF2p2h_eject,
  kFermiMoveEffectiveSF2p2h_noeject,
} FermiMoverInteractionType_t;  

class NuclearModel {

public:
  static const char * AsString(NuclearModel_t nucmod) {
     switch (nucmod) {
     case kNucmUndefined:       return "Undefined nuclear model";  break;
     case kNucmFermiGas:        return "Fermi gas model";          break;
     case kNucmLocalFermiGas:   return "Local Fermi gas model";    break;
     case kNucmSpectralFunc:    return "Spectral function model";  break;
     case kNucmEffSpectralFunc: return "Effective spectral function model"; break;
     default:                 break;
     }
     return " ";
  }

};

}
#endif
