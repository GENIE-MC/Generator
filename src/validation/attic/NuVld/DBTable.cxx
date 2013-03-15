//____________________________________________________________________________
/*
 Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory - Jan 12, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Dec 08, 2008 - CA
   NuValidator package refurbishment. Removed the neugen3 dependency. 
   Moved all sources to $GENIE/src/ValidationTools/NuVld.
   Some clases have been renamed.
 @ Aug 25, 2009 - CA
   Removed redundant versions of ParserUtils.h and ParserStatus.h in favor of
   the ones in $GENIE/Conventions and $GENIE/Utils. Updated code accordingly.

*/
//____________________________________________________________________________ 

#include <cassert>
#include <sstream>

#include <TMath.h>
#include <TObjArray.h>

#include "Messenger/Messenger.h"
#include "ValidationTools/NuVld/DBTable.h"
#include "ValidationTools/NuVld/DBNuXSecTableRow.h"
#include "ValidationTools/NuVld/DBNuXSecTableFields.h"
#include "ValidationTools/NuVld/DBElDiffXSecTableRow.h"
#include "ValidationTools/NuVld/DBElDiffXSecTableFields.h"
#include "ValidationTools/NuVld/DBSFTableRow.h"
#include "ValidationTools/NuVld/DBSFTableFields.h"

using std::cout;
using std::endl;
using std::ostringstream;

namespace genie {
namespace nuvld {

//____________________________________________________________________________
/*
template<class T>
           ostream & operator << (ostream & stream, const DBTable<T> & table)
{
  typename vector<T *>::const_iterator row_iter;

  for(row_iter = table.fTable.begin();
             row_iter != table.fTable.end(); ++row_iter) cout << *(*row_iter);

  return stream;
}*/
//____________________________________________________________________________
template<> DBTable<DBElDiffXSecTableRow>::DBTable()
{
  fFields  = new DBElDiffXSecTableFields();
  fQueryStr = 0;
  fIdList = 0;
}
//____________________________________________________________________________
template<> DBTable<DBNuXSecTableRow>::DBTable()
{
  fFields   = new DBNuXSecTableFields();
  fQueryStr = 0;
  fIdList   = 0;
}
//____________________________________________________________________________
template<> DBTable<DBSFTableRow>::DBTable()
{
  fFields   = new DBSFTableFields();
  fQueryStr = 0;
  fIdList   = 0;
}
//____________________________________________________________________________
template<class T> DBTable<T>::DBTable(const DBTable<T> * table)
{
  LOG("NuVld", pDEBUG) << "Copying DBTableFields...";

  fFields  = new DBTableFields(table->fFields);

  LOG("NuVld", pDEBUG) << "Copying DBTableRows...";

  typename vector<T *>::const_iterator iter;
  for(iter = table->fTable.begin();
                   iter != table->fTable.end(); ++iter) this->AddRow( *iter );

  LOG("NuVld", pDEBUG) << "Copying DBMeasurementIdList...";

  fIdList = new DBMeasurementIdList(table->fIdList);

  LOG("NuVld", pDEBUG) << "Copying DBQueryString...";

  if( table->QueryString() )
     fQueryStr = new DBQueryString( *(table->QueryString()) );
  else {
     LOG("NuVld", pWARN) << "No DBQueryString attached - check the DBI";
     fQueryStr = 0;
  }

  LOG("NuVld", pDEBUG) << "...successfull";
}
//____________________________________________________________________________
template<class T> DBTable<T>::~DBTable()
{
  if(fIdList) delete fIdList;
}
//____________________________________________________________________________
template<class T> int DBTable<T>::NRows(void) const
{
  return (int) fTable.size();
}
//____________________________________________________________________________
template<class T> void DBTable<T>::AddRow(DBTableRow * row)
{
  T * trow = dynamic_cast<T *> (row);
  fTable.push_back( trow );
}
//____________________________________________________________________________
template<class T> const DBTableRow * DBTable<T>::Row(int irow) const
{
   if(irow >=0 && irow < this->NRows())
                            return dynamic_cast <DBTableRow *> (fTable[irow]);
   else return 0;
}
//____________________________________________________________________________
template<class T> void DBTable<T>::MergeWithTable(DBTableBase * xsectb)
{
   DBTable<T> * xsect = dynamic_cast<DBTable<T> *> (xsectb);
   typename vector<T *>::iterator row_iter;

   for(row_iter = xsect->fTable.begin();
                             row_iter != xsect->fTable.end(); ++row_iter) {
       T * row = *row_iter;
       this->AddRow(row);
   }
}
//____________________________________________________________________________
template<class T> DBTable<T> * DBTable<T>::Subset(
                              string experiment, string measurement_tag) const
{
// Create a subset DBTable with the measurements corresponding to an
// (experiment name, measurement tag) pair.
//
  LOG("NuVld", pDEBUG)
              << "Slicing DBTable<T> for exp = "
                             << experiment << ", mtag = " << measurement_tag;

  DBTable * subset_table = new DBTable<T>();

  typename vector<T *>::const_iterator row_iter;
  for(row_iter = fTable.begin(); row_iter != fTable.end(); ++row_iter) {

    T * row = *row_iter;
    if(experiment.compare(row->Experiment()) == 0 &&
            measurement_tag.compare(row->XmlMeasurementTag() ) == 0 ) {
                subset_table->AddRow(row);
    }
  }
  SLOG("NuVld", pDEBUG)
             << "Sliced DBTable<T> has " << subset_table->NRows() << " rows";

  return subset_table;
}
//____________________________________________________________________________
template<class T> void DBTable<T>::SetQueryString (
                                           const DBQueryString & query_string)
{
  if(fQueryStr) delete fQueryStr;
  fQueryStr = new DBQueryString( query_string );
}
//____________________________________________________________________________
template<class T> void DBTable<T>::SetDBMeasurementIdList(
                                                  DBMeasurementIdList * id_list)
{
  if(fIdList) delete fIdList;
  fIdList = id_list;
}
//____________________________________________________________________________
template<class T> void DBTable<T>::SaveQueryStringToFile(
                                         TDirectory * dir, string name) const
{
// Saving the DBTable<T> to a ROOT file. There is no need to save the actual
// data - the RDBMS should remain the
// We only need to save the query that was used to fill this table.

  dir->cd();

  if( fQueryStr ) {
     LOG("NuVld", pDEBUG) << "\n Writing DBQueryString: " << *fQueryStr;
     fQueryStr->Write(name.c_str());
  } else {
     LOG("NuVld", pERROR)
              << "No DBQueryString attached - DBTable Can not be saved";
  }
}
//____________________________________________________________________________
template<> TGraphAsymmErrors * DBTable<DBNuXSecTableRow>::GetGraph(
                                           const char * opt, const char * var)
{
  string option(opt);

  const vector<DBNuXSecTableRow *> & rows = this->Rows();

  const int npoints = rows.size();

  if(npoints == 0) return 0;

  double * e          = new double[npoints];
  double * e_err_m    = new double[npoints];
  double * e_err_p    = new double[npoints];
  double * xsec       = new double[npoints];
  double * xsec_err_m = new double[npoints];
  double * xsec_err_p = new double[npoints];

  vector<DBNuXSecTableRow *>::const_iterator row_iter;

  int ipoint = 0;

  for(row_iter = rows.begin(); row_iter != rows.end(); ++row_iter) {

     double E     = (*row_iter)->E();
     double E_min = (*row_iter)->Emin();
     double E_max = (*row_iter)->Emax();

     e[ipoint]      = E;
     xsec[ipoint]   = (*row_iter)->XSec();

     //---- check if we need to compute the errors on E
     if( option.find("noE") != string::npos ) {
        e_err_m[ipoint] = 0.;
        e_err_p[ipoint] = 0.;
     } else {
        e_err_m[ipoint] = E - E_min;
        e_err_p[ipoint] = E_max - E;
     }

     //----- check option = type of errors shown on xsec
     //      - all    : sqrt( stat_err^2 + syst_err^2 )
     //      - stat   : statistical error only
     //      - syst   : systematic error only
     //      - noXsec : show no error
     if( option.find("all") != string::npos ) {
        xsec_err_m[ipoint] = (*row_iter)->ErrM();
        xsec_err_p[ipoint] = (*row_iter)->ErrP();

     } else if ( option.find("stat") != string::npos ) {
        xsec_err_m[ipoint] = (*row_iter)->StatErrM();
        xsec_err_p[ipoint] = (*row_iter)->StatErrP();

     } else if ( option.find("syst") != string::npos ) {
        xsec_err_m[ipoint] = (*row_iter)->SystErrM();
        xsec_err_p[ipoint] = (*row_iter)->SystErrP();

     } else {
        xsec_err_m[ipoint] = 0.;
        xsec_err_p[ipoint] = 0.;
     }

     //------ check if we plot xsec or xsec/E
     if( option.find("scale-with-E") != string::npos ) {
          xsec[ipoint]       /= e[ipoint];
          xsec_err_m[ipoint] /= e[ipoint];
          xsec_err_p[ipoint] /= e[ipoint];
     }
     ipoint++;
  }
  return new TGraphAsymmErrors(npoints, e, xsec,
                                    e_err_m, e_err_p, xsec_err_m, xsec_err_p);
}
//____________________________________________________________________________
template<> TGraphAsymmErrors * DBTable<DBElDiffXSecTableRow>::GetGraph(
    const char * opt, const char * var)
{
  const vector<DBElDiffXSecTableRow *> & rows = this->Rows();

  const int npoints = rows.size();

  if(npoints == 0) return 0;

  double * x      = new double[npoints];
  double * dx     = new double[npoints];
  double * Sigma  = new double[npoints];
  double * dSigma = new double[npoints];

  vector<DBElDiffXSecTableRow *>::const_iterator row_iter;

  int ipoint = 0;

  for(row_iter = rows.begin(); row_iter != rows.end(); ++row_iter) {

     Sigma[ipoint] = (*row_iter)->Sigma();
     dx[ipoint]    = 0;

     //y-err
     if(strcmp(opt,"err") == 0)
         dSigma[ipoint] = (*row_iter)->dSigma();
     else
         dSigma[ipoint] = 0;

     //x-var
     if      ( strcmp(var,"E")       == 0)  x[ipoint] = (*row_iter)->E();
     else if ( strcmp(var,"EP")      == 0)  x[ipoint] = (*row_iter)->EP();
     else if ( strcmp(var,"Theta")   == 0)  x[ipoint] = (*row_iter)->Theta();
     else if ( strcmp(var,"Q2")      == 0)  x[ipoint] = (*row_iter)->Q2();
     else if ( strcmp(var,"W2")      == 0)  x[ipoint] = (*row_iter)->W2();
     else if ( strcmp(var,"Nu")      == 0)  x[ipoint] = (*row_iter)->Nu();
     else if ( strcmp(var,"Epsilon") == 0)  x[ipoint] = (*row_iter)->Epsilon();
     else if ( strcmp(var,"Gamma")   == 0)  x[ipoint] = (*row_iter)->Gamma();
     else if ( strcmp(var,"x")       == 0)  x[ipoint] = (*row_iter)->x();
     else                                   x[ipoint] = 0;

     ipoint++;
  }
  return new TGraphAsymmErrors(npoints, x, Sigma, dx, dx, dSigma, dSigma);
}
//____________________________________________________________________________
template<> MultiGraph * DBTable<DBNuXSecTableRow>::GetMultiGraph(
                                           const char * opt, const char * var)
{
  LOG("NuVld", pDEBUG) << "Getting multi-graph";

  string option(opt);

  MultiGraph * mgraph = new MultiGraph();
  const DBMeasurementIdList * id_list = this->IdList();

  for(unsigned int i = 0; i < id_list->NIds(); i++) {

    LOG("NuVld", pDEBUG) << "N-IDS = " << id_list->NIds();

    string exp = id_list->GetId(i)->Experiment();
    string tag = id_list->GetId(i)->XmlMeasurementTag();

    LOG("NuVld", pDEBUG)
             << "Getting/Adding graph..." << i << ": " << exp << "/" << tag;

    DBTable<DBNuXSecTableRow> * subset = this->Subset(exp, tag);
    TGraphAsymmErrors *      graph  = subset->GetGraph(opt, var);

    mgraph->AddGraph(id_list->GetId(i)->Reference().c_str(), graph);
    delete subset;
  }

  LOG("NuVld", pDEBUG) << "ngraphs = " << mgraph->NGraphs();
  return mgraph;
}
//____________________________________________________________________________
template<> MultiGraph * DBTable<DBElDiffXSecTableRow>::GetMultiGraph(
   const char * opt, const char * var)
{
  MultiGraph * mgraph = new MultiGraph();
  const DBMeasurementIdList * id_list = this->IdList();

  for(unsigned int i = 0; i < id_list->NIds(); i++) {

    string exp = id_list->GetId(i)->Experiment();
    string tag = id_list->GetId(i)->XmlMeasurementTag();

    LOG("NuVld", pDEBUG)
             << "Getting/Adding graph..." << i << ": " << exp << "/" << tag;

    DBTable<DBElDiffXSecTableRow> * subset = this->Subset(exp, tag);
    TGraphAsymmErrors *  graph  = subset->GetGraph(opt, var);

    mgraph->AddGraph(id_list->GetId(i)->Reference().c_str(), graph);
    delete subset;
  }

  LOG("NuVld", pDEBUG) << "ngraphs = " << mgraph->NGraphs();
  return mgraph;
}
//____________________________________________________________________________
template<> TGraphAsymmErrors * DBTable<DBSFTableRow>::GetGraph(
    const char * opt, const char * var)
{
  const vector<DBSFTableRow  *> & rows = this->Rows();

  const int npoints = rows.size();

  if(npoints == 0) return 0;

  double * x    = new double[npoints];
  double * dx   = new double[npoints];
  double * y    = new double[npoints];
  double * dyp  = new double[npoints];
  double * dym  = new double[npoints];

  vector<DBSFTableRow *>::const_iterator row_iter;

  int ipoint = 0;

  for(row_iter = rows.begin(); row_iter != rows.end(); ++row_iter) {

     //x-var
     if      ( strcmp(var,"Q2") == 0)  x[ipoint] = (*row_iter)->Q2();
     else if ( strcmp(var,"x" ) == 0)  x[ipoint] = (*row_iter)->x();
     else                              x[ipoint] = 0;

     y  [ipoint] = (*row_iter)->SF();
     dx [ipoint] = 0;
     dyp[ipoint] = 0;
     dym[ipoint] = 0;

     //y-err
     if( strcmp(opt,"err") == 0) {
          dyp[ipoint] = (*row_iter)->ErrP();
          dym[ipoint] = (*row_iter)->ErrM();
     }
     ipoint++;
  }
  return new TGraphAsymmErrors(npoints, x, y, dx, dx, dym, dyp);
}
//____________________________________________________________________________
template<> MultiGraph * DBTable<DBSFTableRow>::GetMultiGraph(
                                           const char * opt, const char * var)
{
  MultiGraph * mgraph = new MultiGraph();
  const DBMeasurementIdList * id_list = this->IdList();

  for(unsigned int i = 0; i < id_list->NIds(); i++) {

    string exp = id_list->GetId(i)->Experiment();
    string tag = id_list->GetId(i)->XmlMeasurementTag();

    LOG("NuVld", pDEBUG)
       << "Getting/Adding graph..." << i << ": " << exp << "/" << tag;

    DBTable<DBSFTableRow> * subset = this->Subset(exp, tag);
    TGraphAsymmErrors * graph  = subset->GetGraph(opt, var);

    mgraph->AddGraph(id_list->GetId(i)->Reference().c_str(), graph);
    delete subset;
  }

  LOG("NuVld", pDEBUG) << "ngraphs = " << mgraph->NGraphs();
  return mgraph;
}
//____________________________________________________________________________

// template specializations:
//
template class DBTable<DBNuXSecTableRow>;
template class DBTable<DBElDiffXSecTableRow>;
template class DBTable<DBSFTableRow>;

} // nuvld namespace
} // genie namespace

