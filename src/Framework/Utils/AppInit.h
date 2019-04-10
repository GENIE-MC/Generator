//____________________________________________________________________________
/*!

\namespace  genie::utils::app_init

\brief      Initialization code commonly occuring in GENIE apps,
            factored out from existing apps for convenience.
            Not generic GENIE initialization code.

\author     Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
            University of Liverpool & STFC Rutherford Appleton Lab

\created    January 31, 2013

\cpright    Copyright (c) 2003-2019, The GENIE Collaboration
            For the full text of the license visit http://copyright.genie-mc.org
            or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _APP_INIT_UTILS_H_
#define _APP_INIT_UTILS_H_

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
