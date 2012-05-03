//_____________________________________________________________________________
/*!

\class    genie::nuvld::SqlUtils

\brief

\author   Costas Andreopoulos (Rutherford Lab.)  <costas.andreopoulos \at stfc.ac.uk>

\created  January 2004          
*/
//_____________________________________________________________________________

#ifndef _SQL_UTILS_H_
#define _SQL_UTILS_H_

#include <sstream>
#include <string>
#include <iostream>

#include <TSQLServer.h>
#include <TSQLResult.h>
#include <TSQLRow.h>

#include "Utils/StringUtils.h"

using std::ostringstream;
using std::string;
using std::vector;
using std::cout;
using std::endl;

namespace genie {
namespace nuvld {

class SqlUtils {

public:

 //__________________________________________________________________________
 static string build_v_conditional(
                string experiments, string xsecs, string nus, string targets)
 {
   unsigned int i = 0;

   ostringstream conditional;

   vector<string>::iterator str_iter;

   vector<string> vec_exp_names, vec_observables, vec_nu, vec_tgt;

   vec_exp_names   = utils::str::Split(experiments.c_str(),  ",");
   vec_observables = utils::str::Split(xsecs.c_str(),        ",");
   vec_nu          = utils::str::Split(nus.c_str(),          ",");
   vec_tgt         = utils::str::Split(targets.c_str(),      ",");

   //-- add "experiments"-selection

   if(strcmp(experiments.c_str(), "*") != 0
                                        && experiments.size() > 0) {
     i=0;
     conditional << "MEASUREMENT_HEADER.name IN (";
     for(str_iter = vec_exp_names.begin();
                              str_iter != vec_exp_names.end(); ++str_iter) {

      conditional << "\"" << utils::str::FilterString(" ", *str_iter) << "\"";

      if(i < vec_exp_names.size()-1 ) conditional << ", ";
      i++;
     }
     conditional << ")";
   } else if ( strcmp(experiments.c_str(), "*") == 0)  conditional << " 1 ";
   else conditional << " 0 ";

   //-- add "observables"-selection

   if( strcmp(xsecs.c_str(), "*") != 0  && xsecs.size() > 0) {
     i=0;
     conditional << " AND MEASUREMENT_HEADER.observable IN (";

     for(str_iter = vec_observables.begin();
            str_iter != vec_observables.end(); ++str_iter) {

       conditional << "\"" << utils::str::FilterString(" ", *str_iter) << "\"";

       if(i < vec_observables.size()-1 ) conditional << ", ";
       i++;
     }
     conditional << ")";
   } else if (strcmp(xsecs.c_str(), "*") == 0) conditional << " AND 1 ";
   else conditional << " AND 0 ";

   //-- add "neutrino"-selection

   if( strcmp(nus.c_str(),"*") != 0  && nus.size() > 0) {
    i=0;
    conditional << " AND (";

    for(str_iter = vec_nu.begin(); str_iter != vec_nu.end(); ++str_iter) {

       conditional << " MEASUREMENT_HEADER.reaction LIKE \"%"
                   << utils::str::FilterString(" ", *str_iter) << " %\" ";

       if(i < vec_nu.size()-1 ) conditional << " OR ";
       i++;
    }
    conditional << ") ";
   } else if (strcmp(nus.c_str(), "*") == 0) conditional << " AND 1 ";
   else conditional << " AND 0 ";

   //-- add "target"-selection

   if(strcmp(targets.c_str(),"*") != 0 && targets.size() > 0)  {
    i=0;
    conditional << " AND (";

    for(str_iter = vec_tgt.begin(); str_iter != vec_tgt.end(); ++str_iter) {

         conditional << "MEASUREMENT_HEADER.target LIKE \"%"
                     << utils::str::FilterString(" ", *str_iter) << "%\"";

         if(i < vec_tgt.size()-1 ) conditional << " OR ";
         i++;
     }
     conditional << ") ";

   } else if (strcmp(targets.c_str(), "*") == 0) conditional << " AND 1 ";
   else conditional << " AND 0 ";

   return conditional.str();
 }
 //__________________________________________________________________________
 static string build_v_key_list(TSQLServer * db,
                string experiments, string xsecs, string nus, string targets)
 {
   string conditional = build_v_conditional(experiments, xsecs, nus, targets);

   ostringstream query;

   query << "SELECT MEASUREMENT_HEADER.name, MEASUREMENT_HEADER.measurement_tag"
         << " FROM MEASUREMENT_HEADER WHERE " << conditional << ";";

   //cout << query.str() << endl;
         
   TSQLResult * result = db->Query(query.str().c_str());

   const int nrows = result->GetRowCount();

   TSQLRow * row = 0;

   ostringstream key_list;

   for (int i = 0; i < nrows; i++) {

      row = result->Next();

      key_list << row->GetField(0) << "," << row->GetField(1) << ";";

      delete row;
  }
  delete result;

  cout << key_list.str() << endl;

  return key_list.str();
 }
 //__________________________________________________________________________ 
 static string build_e_conditional(string experiments, string targets)
 {
   unsigned int i = 0;

   ostringstream conditional;

   vector<string>::iterator str_iter;

   vector<string> vec_exp_names, vec_tgt;

   vec_exp_names   = utils::str::Split(experiments.c_str(),  ",");
   vec_tgt         = utils::str::Split(targets.c_str(),      ",");

   //-- add "experiments"-selection

   if(strcmp(experiments.c_str(), "*") != 0
                                        && experiments.size() > 0) {
     i=0;
     conditional << "MEASUREMENT_HEADER.name IN (";
     for(str_iter = vec_exp_names.begin();
           str_iter != vec_exp_names.end(); ++str_iter) {

      conditional << "\"" << utils::str::FilterString(" ", *str_iter) << "\"";

      if(i < vec_exp_names.size()-1 ) conditional << ", ";
      i++;
     }
     conditional << ")";
   } else if ( strcmp(experiments.c_str(), "*") == 0)  conditional << " 1 ";
   else conditional << " 0 ";

   //-- add "target"-selection

   if(strcmp(targets.c_str(),"*") != 0 && targets.size() > 0)  {
    i=0;
    conditional << " AND (";

    for(str_iter = vec_tgt.begin(); str_iter != vec_tgt.end(); ++str_iter) {

         conditional << "MEASUREMENT_HEADER.target LIKE \"%"
                     << utils::str::FilterString(" ", *str_iter) << "%\"";

         if(i < vec_tgt.size()-1 ) conditional << " OR ";
         i++;
     }
     conditional << ") ";

   } else if (strcmp(targets.c_str(), "*") == 0) conditional << " AND 1 ";
   else conditional << " AND 0 ";

   return conditional.str();
 }
 //__________________________________________________________________________
 static string build_e_key_list(TSQLServer * db,
                                 string experiments, string targets)
 {
   string conditional = build_e_conditional(experiments, targets);

   ostringstream query;

   query << "SELECT MEASUREMENT_HEADER.name, MEASUREMENT_HEADER.measurement_tag"
         << " FROM MEASUREMENT_HEADER WHERE " << conditional << ";";

   //cout << query.str() << endl;

   TSQLResult * result = db->Query(query.str().c_str());

   const int nrows = result->GetRowCount();

   TSQLRow * row = 0;

   ostringstream key_list;

   for (int i = 0; i < nrows; i++) {

      row = result->Next();

      key_list << row->GetField(0) << "," << row->GetField(1) << ";";

      delete row;
  }
  delete result;

  cout << key_list.str() << endl;

  return key_list.str();
 }
 //__________________________________________________________________________
 static string build_sf_conditional(
                string experiments, string sf, string probes, string targets)
 {
   unsigned int i = 0;

   ostringstream conditional;

   vector<string>::iterator str_iter;

   vector<string> vec_exp    = utils::str::Split(experiments.c_str(),  ",");
   vector<string> vec_probes = utils::str::Split(probes.c_str(),       ",");
   vector<string> vec_tgt    = utils::str::Split(targets.c_str(),      ",");

   //-- add "experiments"-selection

   if(strcmp(experiments.c_str(), "*") != 0 && experiments.size() > 0) {

     i=0;
     conditional << "MEASUREMENT_HEADER.name IN (";
     for(str_iter = vec_exp.begin(); str_iter != vec_exp.end(); ++str_iter) {

      conditional << "\"" << utils::str::FilterString(" ", *str_iter) << "\"";

      if(i < vec_exp.size()-1 ) conditional << ", ";
      i++;
     }
     conditional << ")";
   } else if ( strcmp(experiments.c_str(), "*") == 0)  conditional << " 1 ";
   else conditional << " 0 ";

   //-- add "probes"-selection

   if( strcmp(probes.c_str(),"*") != 0  && probes.size() > 0) {
    i=0;
    conditional << " AND (";

    for(str_iter = vec_probes.begin(); str_iter != vec_probes.end(); ++str_iter) {

       conditional << " MEASUREMENT_HEADER.reaction LIKE \"%"
                   << utils::str::FilterString(" ", *str_iter) << " %\" ";

       if(i < vec_probes.size()-1 ) conditional << " OR ";
       i++;
    }
    conditional << ") ";
   } else if (strcmp(probes.c_str(), "*") == 0) conditional << " AND 1 ";
   else conditional << " AND 0 ";

   //-- add "target"-selection

   if(strcmp(targets.c_str(),"*") != 0 && targets.size() > 0)  {
    i=0;
    conditional << " AND (";

    for(str_iter = vec_tgt.begin(); str_iter != vec_tgt.end(); ++str_iter) {

         conditional << "MEASUREMENT_HEADER.target LIKE \"%"
                    << utils::str::FilterString(" ", *str_iter) << "%\"";

         if(i < vec_tgt.size()-1 ) conditional << " OR ";
         i++;
     }
     conditional << ") ";

   } else if (strcmp(targets.c_str(), "*") == 0) conditional << " AND 1 ";
   else conditional << " AND 0 ";

   //-- add "observables"-selection

   assert(sf.size() > 0);   
   conditional << " AND MEASUREMENT_HEADER.observable = \"" << sf << "\"";

   cout << conditional.str() << endl;

   return conditional.str();   
 }
 //__________________________________________________________________________
 static string build_sf_key_list(TSQLServer * db,
                 string experiments, string sf, string probe, string targets)
 {
   string conditional = build_sf_conditional(experiments, sf, probe, targets);

   ostringstream query;

   query << "SELECT MEASUREMENT_HEADER.name, MEASUREMENT_HEADER.measurement_tag"
         << " FROM MEASUREMENT_HEADER WHERE " << conditional << ";";

   //cout << query.str() << endl;

   TSQLResult * result = db->Query(query.str().c_str());

   const int nrows = result->GetRowCount();

   TSQLRow * row = 0;

   ostringstream key_list;

   for (int i = 0; i < nrows; i++) {

      row = result->Next();

      key_list << row->GetField(0) << "," << row->GetField(1) << ";";

      delete row;
  }
  delete result;

  cout << key_list.str() << endl;

  return key_list.str();
 }
 //__________________________________________________________________________

};

} // nuvld namespace
} // genie namespace

#endif 

