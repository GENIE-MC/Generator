//____________________________________________________________________________
/*!

\class    genie::NtpMCFormat

\brief    Encapsulates an enumeration of possible GENIE output TTree formats

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  September 02, 2005

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
   kNFPlainRecord,
   kNFEventRecord,

} NtpMCFormat_t;


class NtpMCFormat {

 public:

  static char * AsString(NtpMCFormat_t fmt) {
     switch (fmt) {
     case kNFUndefined:
              return "Undefined";
              break;
     case kNFPlainRecord:
              return "TTree with [TBranch=NtpMCPlainRecord]";
              break;
     case kNFEventRecord:
              return "TTree with [TBranch=NtpMCEventRecord]";
              break;
     default:
              break;
     }
     return " ";
  }

  static char * SuggestedUsage(NtpMCFormat_t fmt) {
     switch (fmt) {
     case kNFUndefined:
              return "for wasting your time";
              break;
     case kNFPlainRecord:
              return "for analyzing data in bare ROOT session (without GENIE libs)";
              break;
     case kNFEventRecord:
              return "for passing GENIE data between applications (eg GENIE->GEANT)";
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
     case kNFPlainRecord:
              return "PR";
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
