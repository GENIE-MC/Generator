//____________________________________________________________________________
/*!

\class    genie::geant::GPrimaryGeneratorType

\brief    Encapsulates an enumeration of possible GEANT/GEANT PrimaryGenerator
          types

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  September 12, 2005
 
*/
//____________________________________________________________________________

#ifndef _G_PRIMARY_GENERATOR_TYPE_H_
#define _G_PRIMARY_GENERATOR_TYPE_H_

#ifndef ROOT_Rtypes
#include "Rtypes.h"
#endif

namespace genie {
namespace geant {

typedef enum EGPrimaryGeneratorType {

   kPGUndefined   = -1, 
   kPGNtpRd,
   kPGEVGDrv,
   kPGMCJDrv

} GPrimaryGeneratorType_t; 
  

class GPrimaryGeneratorType {

 public:

  static char * AsString(GPrimaryGeneratorType_t pgt) {
     switch (pgt) {
     case kPGUndefined: 
         return "Undefined PrimaryGenerator";                                
         break;
     case kPGNtpRd:     
         return "NTP PrimaryGenerator [Reads events from GENIE's ER ROOT Trees]"; 
         break;
     case kPGEVGDrv:    
         return "EVG PrimaryGenerator [Invokes GENIE's GEVGDriver]";              
         break;
     case kPGMCJDrv:    
         return "MCJ PrimaryGenerator [Invokes GENIE's GMCJDriver]";              
         break;
     default:           
         break;
     }
     return "---";
  }

};

} // geant namespace
} // genie namespace

#endif
