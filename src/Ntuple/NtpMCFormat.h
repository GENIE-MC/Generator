//____________________________________________________________________________
/*!

\class    genie::NtpMCFormat

\brief    Encapsulates an enumeration of possible GENIE output TTree formats

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  September 02, 2005

\cpright  Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _NTP_MC_FORMAT_H_
#define _NTP_MC_FORMAT_H_

#ifndef ROOT_Rtypes
#include "Rtypes.h"
#endif

namespace genie {

typedef enum ENtpMCFormat {

   kNFUndefined = -1,
   kNFGHEP   /* each mc tree leaf contains the full GHEP EventRecord */

} NtpMCFormat_t;

class NtpMCFormat {
 public:
  static char * AsString(NtpMCFormat_t fmt) {
     switch (fmt) {
     case kNFUndefined:
              return "Undefined";
              break;
     case kNFGHEP:
              return "[NtpMCEventRecord]";
              break;
     default:
              break;
     }
     return " ";
  }

  static char * FilenameTag(NtpMCFormat_t fmt) {

     // The output ROOT files containing GENIE ntuple are typically named as
     // gntp.[tag].root where TAG describes the tree format

     switch (fmt) {
     case kNFUndefined:
              return "undef";
              break;
     case kNFGHEP:
              return "ghep";
              break;
     default:
              break;
     }
     return "undef";
  }
};

}
#endif
