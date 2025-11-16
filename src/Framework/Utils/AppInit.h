//____________________________________________________________________________
/*!

\namespace  genie::utils::app_init

\brief      Initialization code commonly occuring in GENIE apps,
            factored out from existing apps for convenience.
            Not generic GENIE initialization code.

\author     Costas Andreopoulos <c.andreopoulos \at cern.ch>
            University of Liverpool

\created    January 31, 2013

\cpright    Copyright (c) 2003-2025, The GENIE Collaboration
            For the full text of the license visit http://copyright.genie-mc.org          
*/
//____________________________________________________________________________

#ifndef _APP_INIT_UTILS_H_
#define _APP_INIT_UTILS_H_

#include <string>
using std::string;

namespace genie {
namespace utils {

namespace app_init
{
  void RandGen        (long int seed);
  void XSecTable      (string inpfile, bool require_table);
  void MesgThresholds (string inpfile);
  void CacheFile      (string inpfile);

} // app_init namespace
} // utils namespace
} // genie namespace

#endif // _APP_INIT_UTILS_H_
