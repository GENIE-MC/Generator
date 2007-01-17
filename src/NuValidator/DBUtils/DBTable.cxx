//_____________________________________________________________________________
/*!

\class    genie::nuvld::DBTable

\brief    A NuVld data-base table.

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>

\created  January 2004
*/
//_____________________________________________________________________________

#include <cassert>
#include <sstream>

#include <TMath.h>
#include <TObjArray.h>

#include "DBUtils/DBTable.h"
#include "DBUtils/vXSecTableRow.h"
#include "DBUtils/vXSecTableFields.h"
#include "DBUtils/eDiffXSecTableRow.h"
#include "DBUtils/eDiffXSecTableFields.h"
#include "DBUtils/SFTableRow.h"
#include "DBUtils/SFTableFields.h"
#include "Messenger/Messenger.h"
#include "XmlParser/ParserUtils.h"

using std::cout;
using std::endl;
using std::ostringstream;

namespace genie {
namespace nuvld {

//template ostream & operator
//                << (ostream & stream, const DBTable<vXSecTableRow> & table);
//template ostream & operator
//            << (ostream & stream, const DBTable<eDiffXSecTableRow> & table);
//template ostream & operator
//                   << (ostream & stream, const DBTable<SFTableRow> & table);

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
template<> DBTable<eDiffXSecTableRow>::DBTable()
{
  fFields  = new eDiffXSecTableFields();
  fQueryStr = 0;
  fIdList = 0;
}
//____________________________________________________________________________
template<> DBTable<vXSecTableRow>::DBTable()
{
  fFields   = new vXSecTableFields();
  fQueryStr = 0;
  fIdList   = 0;
}
//____________________________________________________________________________
template<> DBTable<SFTableRow>::DBTable()
{
  fFields   = new SFTableFields();
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

  LOG("NuVld", pDEBUG) << "Copying MeasurementIdList...";

  fIdList = new MeasurementIdList(table->fIdList);

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
            measurement_tag.compare(row->MeasurementTag() ) == 0 ) {
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
template<class T> void DBTable<T>::SetMeasurementIdList(
                                                  MeasurementIdList * id_list)
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
template<> TGraphAsymmErrors * DBTable<vXSecTableRow>::GetGraph(
                                           const char * opt, const char * var)
{
  string option(opt);

  const vector<vXSecTableRow *> & rows = this->Rows();

  const int npoints = rows.size();

  if(npoints == 0) return 0;

  double * e          = new double[npoints];
  double * e_err_m    = new double[npoints];
  double * e_err_p    = new double[npoints];
  double * xsec       = new double[npoints];
  double * xsec_err_m = new double[npoints];
  double * xsec_err_p = new double[npoints];

  vector<vXSecTableRow *>::const_iterator row_iter;

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
template<> TGraphAsymmErrors * DBTable<eDiffXSecTableRow>::GetGraph(
                                           const char * opt, const char * var)
{
  const vector<eDiffXSecTableRow *> & rows = this->Rows();

  const int npoints = rows.size();

  if(npoints == 0) return 0;

  double * x      = new double[npoints];
  double * dx     = new double[npoints];
  double * Sigma  = new double[npoints];
  double * dSigma = new double[npoints];

  vector<eDiffXSecTableRow *>::const_iterator row_iter;

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
template<> MultiGraph * DBTable<vXSecTableRow>::GetMultiGraph(
                                           const char * opt, const char * var)
{
  LOG("NuVld", pDEBUG) << "Getting multi-graph";

  string option(opt);

  MultiGraph * mgraph = new MultiGraph();
  const MeasurementIdList * id_list = this->IdList();

  for(unsigned int i = 0; i < id_list->NIds(); i++) {

    LOG("NuVld", pDEBUG) << "N-IDS = " << id_list->NIds();

    string exp = id_list->GetId(i)->Experiment();
    string tag = id_list->GetId(i)->MeasurementTag();

    LOG("NuVld", pDEBUG)
             << "Getting/Adding graph..." << i << ": " << exp << "/" << tag;

    DBTable<vXSecTableRow> * subset = this->Subset(exp, tag);
    TGraphAsymmErrors *      graph  = subset->GetGraph(opt, var);

    mgraph->AddGraph(id_list->GetId(i)->Reference().c_str(), graph);
    delete subset;
  }

  LOG("NuVld", pDEBUG) << "ngraphs = " << mgraph->NGraphs();
  return mgraph;
}
//____________________________________________________________________________
template<> MultiGraph * DBTable<eDiffXSecTableRow>::GetMultiGraph(
                                           const char * opt, const char * var)
{
  MultiGraph * mgraph = new MultiGraph();
  const MeasurementIdList * id_list = this->IdList();

  for(unsigned int i = 0; i < id_list->NIds(); i++) {

    string exp = id_list->GetId(i)->Experiment();
    string tag = id_list->GetId(i)->MeasurementTag();

    LOG("NuVld", pDEBUG)
             << "Getting/Adding graph..." << i << ": " << exp << "/" << tag;

    DBTable<eDiffXSecTableRow> * subset = this->Subset(exp, tag);
    TGraphAsymmErrors *          graph  = subset->GetGraph(opt, var);

    mgraph->AddGraph(id_list->GetId(i)->Reference().c_str(), graph);
    delete subset;
  }

  LOG("NuVld", pDEBUG) << "ngraphs = " << mgraph->NGraphs();
  return mgraph;
}
//____________________________________________________________________________
template<> TGraphAsymmErrors * DBTable<SFTableRow>::GetGraph(
                                           const char * opt, const char * var)
{
  const vector<SFTableRow *> & rows = this->Rows();

  const int npoints = rows.size();

  if(npoints == 0) return 0;

  double * x    = new double[npoints];
  double * dx   = new double[npoints];
  double * y    = new double[npoints];
  double * dyp  = new double[npoints];
  double * dym  = new double[npoints];

  vector<SFTableRow *>::const_iterator row_iter;

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
template<> MultiGraph * DBTable<SFTableRow>::GetMultiGraph(
                                           const char * opt, const char * var)
{
  MultiGraph * mgraph = new MultiGraph();
  const MeasurementIdList * id_list = this->IdList();

  for(unsigned int i = 0; i < id_list->NIds(); i++) {

    string exp = id_list->GetId(i)->Experiment();
    string tag = id_list->GetId(i)->MeasurementTag();

    LOG("NuVld", pDEBUG)
       << "Getting/Adding graph..." << i << ": " << exp << "/" << tag;

    DBTable<SFTableRow> * subset = this->Subset(exp, tag);
    TGraphAsymmErrors *   graph  = subset->GetGraph(opt, var);

    mgraph->AddGraph(id_list->GetId(i)->Reference().c_str(), graph);
    delete subset;
  }

  LOG("NuVld", pDEBUG) << "ngraphs = " << mgraph->NGraphs();
  return mgraph;
}
//____________________________________________________________________________

// template specializations:
//
template class DBTable<vXSecTableRow>;
template class DBTable<eDiffXSecTableRow>;
template class DBTable<SFTableRow>;

} // nuvld namespace
} // genie namespace

