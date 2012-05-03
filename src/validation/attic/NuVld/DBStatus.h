//_____________________________________________________________________________
/*!

\class    genie::nuvld::DBStatus

\brief

\author   Costas Andreopoulos (Rutherford Lab.)  <costas.andreopoulos \at stfc.ac.uk>

\created  August 2003
*/
//_____________________________________________________________________________


#ifndef _DB_STATUS_H_
#define _DB_STATUS_H_

namespace genie {
namespace nuvld {

typedef enum EDBStatus {

   eDbu_OK                = 0,
   eDbu_SQL_ERROR         = 1,
   eDbu_CONNECTION_LOST   = 2,
   eDbu_CONNECTION_FAILED = 3,
   eDbu_EMPTY_TABLE       = 4

} DBStatus_t;

class DBStatus {

public:

     static const char * AsString(DBStatus_t status) {
       switch(status) {
         case eDbu_OK:                 return "Successful SQL query";            break;
         case eDbu_SQL_ERROR:          return "SQL syntax error";                break;
         case eDbu_CONNECTION_LOST:    return "Lost connection with RDBMS";      break;
         case eDbu_CONNECTION_FAILED:  return "Connection to RDBMS failed";      break;
         case eDbu_EMPTY_TABLE:        return "SQL query returned empty table";  break;
         default:                      return "Unknown DBStatus_t";              break;
       }
     }
};

} // nuvld namespace
} // genie namespace

#endif 

