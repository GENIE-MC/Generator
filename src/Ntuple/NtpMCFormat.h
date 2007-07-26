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
   kNFEventRecord

} NtpMCFormat_t;


class NtpMCFormat {

 public:

  static char * AsString(NtpMCFormat_t fmt) {
     switch (fmt) {
     case kNFUndefined:
              return "Undefined";
              break;
     case kNFEventRecord:
              return "[NtpMCEventRecord]";
              break;
     default:
              break;
     }
     return " ";
  }

  static char * FilenameTag(NtpMCFormat_t fmt) {

     // The output ROOT files containing GENIE ntuple are typically named as
     // GNtp[TAG].root where TAG describes the tree format (This is just
     // a naming convention for helping out the user with his book-keeping.
     // GENIE understands the TTree format not by checking the filename but
     // by reading the format NtpMCFormat_t variable from the TTree header)

     switch (fmt) {
     case kNFUndefined:
              return "UNDEFINED";
              break;
     case kNFEventRecord:
              return "ER";
              break;
     default:
              break;
     }
     return "UNDEFINED";
  }

};

}
#endif
