//____________________________________________________________________________
/*!

\namespace  genie::utils::system

\brief      System utilities

\author     Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
            University of Liverpool & STFC Rutherford Appleton Lab

\created    Oct 08, 2009

\cpright    Copyright (c) 2003-2019, The GENIE Collaboration
            For the full text of the license visit http://copyright.genie-mc.org
            or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _SYST_UTILS_H_
#define _SYST_UTILS_H_

#include <vector>
#include <string>

using std::vector;
using std::string;

namespace genie {
namespace utils {

namespace system
{

  vector<string> GetAllFilesInPath  (string path, string extension="");

  int GenieMajorVrsNum (string tag);
  int GenieMinorVrsNum (string tag);
  int GenieRevisVrsNum (string tag);

  bool FileExists(string filename);
  bool DirectoryExists( const char * path ) ;

  string LocalTimeAsString(string format);

} // system namespace
} // utils  namespace
} // genie  namespace

#endif // _SYST_UTILS_H_
