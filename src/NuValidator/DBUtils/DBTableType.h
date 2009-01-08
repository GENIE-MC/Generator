//_____________________________________________________________________________
/*!

\class    genie::nuvld::DBTableType

\brief    An enumeration of possible DBTable<T> types.

\author   Costas Andreopoulos (Rutherford Lab.)  <costas.andreopoulos \at stfc.ac.uk>

\created  January 28, 2005
*/
//_____________________________________________________________________________


#ifndef _DB_TABLE_TYPE_H_
#define _DB_TABLE_TYPE_H_

#include <string>

using std::string;

namespace genie {
namespace nuvld {

typedef enum EDBTableType {

   eDbt_NuXSec            = 0,
   eDbt_NuDiffXSec        = 1,
   eDbt_ElXSec            = 2,
   eDbt_ElDiffXSec        = 3,
   eDbt_HadrMult          = 4,
   eDbt_SF                = 5,
   eDbt_Unknown           = 6

} DBTableType_t;

class DBTableType {

public:

  //___________________________________________________________________________
  static string AsString(DBTableType_t table_type) 
  {
    switch(table_type) {

      case eDbt_NuXSec:            return "vN-XSec";       break;
      case eDbt_NuDiffXSec:        return "vN-Diff-XSec";  break;
      case eDbt_ElXSec:            return "eN-XSec";       break;
      case eDbt_ElDiffXSec:        return "eN-Diff-XSec";  break;
      case eDbt_HadrMult:          return "Hadron-Mult";   break;
      case eDbt_SF:                return "SF";            break;
      case eDbt_Unknown:           
      default:                     return "Unknown";
    }
  }
  //___________________________________________________________________________
  static DBTableType_t FromString(string ttype) 
  {
    if     ( strcmp( ttype.c_str(),"vN-XSec" )      == 0 ) return eDbt_NuXSec;
    else if( strcmp( ttype.c_str(),"vN-Diff-XSec" ) == 0 ) return eDbt_NuDiffXSec;
    else if( strcmp( ttype.c_str(),"eN-XSec" )      == 0 ) return eDbt_ElXSec;
    else if( strcmp( ttype.c_str(),"eN-Diff-XSec" ) == 0 ) return eDbt_ElDiffXSec;
    else if( strcmp( ttype.c_str(),"Hadron-Mult" )  == 0 ) return eDbt_HadrMult;
    else if( strcmp( ttype.c_str(),"SF" )           == 0 ) return eDbt_SF;
    else
         return eDbt_Unknown;
  }
  //___________________________________________________________________________

};

} // nuvld namespace
} // genie namespace

#endif 

