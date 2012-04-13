//_____________________________________________________________________________
/*!

\class    genie::nuvld::SqlTag

\brief

\author   Costas Andreopoulos (Rutherford Lab.)  <costas.andreopoulos \at stfc.ac.uk>

\created  January 2004
*/
//_____________________________________________________________________________

#ifndef _SQL_TAG_H_
#define _SQL_TAG_H_

#ifndef ROOT_Rtypes
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "Rtypes.h"
#endif
#endif

namespace genie {
namespace nuvld {

class SqlTag {

  public:

     typedef enum SqlEnum {

       kCreateDatabase_NuScat = 0,
       kUseDatabase_NuScat,
       kDropDatabase_NuScat,
       kDropTable_all,

       kCreateTable_BeamFlux,
       kCreateTable_ExpInfo,
       kCreateTable_Reference,
       kCreateTable_XmlMeasurements,
       kCreateTable_CrossSection,
       kCreateTable_PartialCrossSection,

       kGetXSec_nuDIS,
       kGetXSec_nubarDIS

     } SqlEnum_t;

     static const char * as_string(SqlEnum_t sql) {
       switch(sql) {
         case kCreateDatabase_NuScat:             return "create database NuScat"; break;
         case kUseDatabase_NuScat:                return "use database NuScat"; break;
         case kDropDatabase_NuScat:               return "drop database NuScat"; break;
         case kDropTable_all:                     return "drop all tables"; break;
         case kCreateTable_BeamFlux:              return "create table BEAM_FLUX"; break;
         case kCreateTable_ExpInfo:               return "create table EXP_INFO"; break;
         case kCreateTable_Reference:             return "create table REFERENCE"; break;
         case kCreateTable_CrossSection:          return "create table CROSS_SECTION"; break;
         case kCreateTable_PartialCrossSection:   return "create table PARTIAL_CROSS_SECTION"; break;
         case kCreateTable_XmlMeasurements:          return "create table MEASUREMENTS"; break;
         case kGetXSec_nuDIS:                     return "querying for v+N DIS XSec"; break;
         case kGetXSec_nubarDIS:                  return "querying for vbar+N DIS XSec"; break;
         default:                                 return "unrecognized default sql query"; break;
       }
       return "unrecognized default sql query";
     }
     
     static const char * filename(SqlEnum_t sql) {     
       switch(sql) {
         case kCreateDatabase_NuScat:             return "createDatabase_NuScat.sql"; break;
         case kUseDatabase_NuScat:                return "useDatabase_NuScat.sql"; break;
         case kDropDatabase_NuScat:               return "dropDatabase_NuScat.sql"; break;
         case kDropTable_all:                     return "dropTable_all.sql"; break;
         case kCreateTable_BeamFlux:              return "createTable_BeamFlux.sql"; break;
         case kCreateTable_ExpInfo:               return "createTable_ExpInfo.sql"; break;
         case kCreateTable_Reference:             return "createTable_Reference.sql"; break;
         case kCreateTable_CrossSection:          return "createTable_CrossSection.sql"; break;
         case kCreateTable_PartialCrossSection:   return "createTable_PartialCrossSection.sql"; break;
         case kCreateTable_XmlMeasurements:          return "createTable_XmlMeasurements.sql"; break;
         case kGetXSec_nuDIS:                     return "getXSec_nuDIS.sql"; break;
         case kGetXSec_nubarDIS:                  return "getXSec_nubarDIS.sql"; break;
         default:                                 return "??"; break;
       }
       return "??"; 
     }
};

} // nuvld namespace
} // genie namespace

#endif // SqlTag

