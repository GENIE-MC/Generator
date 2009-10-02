//____________________________________________________________________________
/*!

\program gtune_mafit

\brief   Cross section tuning utility: Fits MA QEL

\syntax  gtune_mafit -h host -u username -p password -d sets_to_fit

         where 
	  -h specifies the MySQL hostname and dbase
	  -u specifies the MySQL username
	  -p specifies the MySQL password
	  -d specifies which data sets to fit 
             (see list below, input as a comma separated list)

         Example:
         gtune_mafit -h mysql://localhost/NuScat -u costas -p mypass1 -d 0,3
                      
\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory

         Hugh Gallagher <gallag \at minos.phy.tufts.edu>
         Tufts University

\created June 06, 2008 

\cpright Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#include <cassert>
#include <cstdlib>
#include <sstream>
#include <algorithm>
#include <string>
#include <vector>

#include <TFile.h>
#include <TMinuit.h>
#include <TMath.h>
#include <TSQLServer.h>

#include "Algorithm/AlgFactory.h"
#include "Algorithm/AlgConfigPool.h"
#include "Base/XSecAlgorithmI.h"
#include "Conventions/Units.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodes.h"
#include "Registry/Registry.h"
#include "Utils/StringUtils.h"
#include "ValidationTools/NuVld/DBStatus.h"
#include "ValidationTools/NuVld/DBI.h"
#include "ValidationTools/NuVld/DBTable.h"
#include "ValidationTools/NuVld/DBQueryString.h"
#include "ValidationTools/NuVld/MultiGraph.h"

using std::vector;
using std::string;
using std::ostringstream;

using namespace genie;
using namespace genie::nuvld;
using namespace genie::utils;

// constants
//
// .............................................
// Fitted data sets :
// ---------------------------------------------
// ID   DESCRIPTION
// 0    neutrino      QEL [all]
// 1    neutrino      QEL [light targets]
// 2    neutrino      QEL [heavy targets]
// 3    anti-neutrino QEL [all]
// 4    anti-neutrino QEL [light targets]
// 5    anti-neutrino QEL [heavy targets]
// ..............................................

const int    kNDataSets = 6;
const double kEmin      =  0.1; // GeV
const double kEmax      = 30.0; // GeV

// keys (experiment id / measurement id) for extracting the
// cross section data from the NuValidator MySQL data-base
const char * kKeyList[kNDataSets] = {
/* 0 */ "ANL_12FT,1;ANL_12FT,3;BEBC,12;BNL_7FT,3;FNAL_15FT,3;Gargamelle,2;SERP_A1,0;SERP_A1,1;SKAT,8",
/* 1 */ "ANL_12FT,1;ANL_12FT,3;BEBC,12;BNL_7FT,3;FNAL_15FT,3",
/* 2 */ "Gargamelle,2;SERP_A1,0;SERP_A1,1;SKAT,8",
/* 3 */ "BNL_7FT,2;Gargamelle,3;Gargamelle,5;SERP_A1,2;SKAT,9",
/* 4 */ "BNL_7FT,2",
/* 5 */ "Gargamelle,3;Gargamelle,5;SERP_A1,2;SKAT,9"
};

// func prototypes
//
void          Init           (int argc, char ** argv);
string        GetArgument    (int argc, char ** argv, const char * option);
double        GetXSecGENIE   (double E, int imode);
void          GetXSecData    (void);
DBQueryString FormQuery      (const char * key_list);
void          DoTheFit       (void);
void          FitFunc        (Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t iflag);
void          Save           (string filename);

// globals
//
vector<int>                 gEnabledDataSets;
DBI *                       gDBI = 0;                   ///< dbase interface
DBTable<DBNuXSecTableRow> * gXSecData     [kNDataSets]; ///< fitted data
MultiGraph *                gXSecDataGrph [kNDataSets]; ///< fitted data as graphs
const XSecAlgorithmI *      gQELXSecModel;
Interaction *               gNuMuN;
Interaction *               gNuMuBarP;

//____________________________________________________________________________
int main(int argc, char ** argv)
{
  // Initialize
  Init(argc, argv);

  // Get data from the NuVld data-base
  GetXSecData();

  // Fit Ma QEL using MINUT
  DoTheFit();

  // Save plots etc
  Save("maqel_fit.root");

  return 0;
}
//____________________________________________________________________________
void Init(int argc, char ** argv)
{
  // get command line arguments
  string url      = GetArgument(argc, argv, "-h");
  string username = GetArgument(argc, argv, "-u");
  string passwd   = GetArgument(argc, argv, "-p");
  string datasets = GetArgument(argc, argv, "-d");

  vector<string> dsvec = str::Split(datasets,",");
  vector<string>::const_iterator it = dsvec.begin();
  for( ; it != dsvec.end(); ++it) {
    gEnabledDataSets.push_back( atoi(it->c_str()) );
  }

  // establish connection with the NuValidator data-base and create a
  // data-base interface
  TSQLServer * sql_server = TSQLServer::Connect(
       url.c_str(), username.c_str(), passwd.c_str());
  assert(sql_server && sql_server->IsConnected());
  gDBI = new DBI(sql_server);
  assert(gDBI);

  // Get the cross section models that we will be using
  AlgFactory * algf = AlgFactory::Instance();
  gQELXSecModel = 
     dynamic_cast<const XSecAlgorithmI *>(
             algf->GetAlgorithm("genie::QELPXSec","CC-Default"));

  LOG("MaQELFit", pNOTICE) 
     << "Got algorithm: " << gQELXSecModel->Id();

  // Create interaction objects to drive the cross section algorithm
  gNuMuN    = Interaction::QELCC(kPdgTgtFreeN, kPdgNeutron, kPdgNuMu);
  gNuMuBarP = Interaction::QELCC(kPdgTgtFreeP, kPdgProton,  kPdgAntiNuMu);
}
//____________________________________________________________________________
string GetArgument(int argc, char ** argv, const char * option)
{
// get the command line argument following the input 'option'

  for(int iarg = 0; iarg < argc-1; iarg++) {
    string argument(argv[iarg]);
    if (argument.compare(option) == 0 ) return string(argv[++iarg]);
  }
  return "";
}
//____________________________________________________________________________
void GetXSecData (void)
{
// download cross section data from the NuVld MySQL dbase

  for(int i=0; i<kNDataSets; i++) {
    gXSecData[i] = new DBTable<DBNuXSecTableRow>;

    DBQueryString query  = FormQuery(kKeyList[i]);
    DBStatus_t    status = gDBI->FillTable(gXSecData[i], query);

    assert(status == eDbu_OK);
  }

 // convert all filled-in DBTables to simple ROOT graphs (in fact multi-graphs: 
 // a collection of graphs, one per different expt/measurement id that makes up
 // the data-base table - can be used if different expt measurements are allowed 
 // to float in the fit).
 for(int i=0; i<kNDataSets; i++) {
    gXSecDataGrph[i] = gXSecData[i]->GetMultiGraph("all-noE");
 }
}
//____________________________________________________________________________
DBQueryString FormQuery(const char * key_list)
{
// forms a DBQueryString for extracting neutrino cross section data from the 
// input key-list and for the input energy range
  
  ostringstream query_string;
  
  query_string 
     << "KEY-LIST:" << key_list
     << "$CUTS:Emin=" << kEmin << ";Emax=" << kEmax
     << "$DRAW_OPT:none$DB-TYPE:vN-XSec";
  
  DBQueryString query(query_string.str());
  
  return query;
}
//____________________________________________________________________________
double GetXSecGENIE(double E, int imode)
{
// get the GENIE cross section, for the current values of the fitted params, 
// for the input mode at the input energy

  Interaction * interaction = 0;

  if      (imode == 0 || imode == 1 || imode == 2) interaction = gNuMuN;
  else if (imode == 3 || imode == 4 || imode == 5) interaction = gNuMuBarP;

  if(!interaction) return 0;

  interaction -> InitStatePtr() -> SetProbeE(E);

  double xsec_qel = gQELXSecModel->Integral(interaction);
  xsec_qel /= (1E-38 * units::cm2);

  return xsec_qel;
}
//____________________________________________________________________________
void DoTheFit(void)
{
  // Initialize MINUIT
  //
  const int np = 1;
  TMinuit * minuit = new TMinuit(np);

  double arglist[10];
  int ierrflg = 0;

  arglist[0] = 1;
  minuit->mnexcm("SET ERR",arglist,1,ierrflg);

  float        value [np] = {  1.00   };
  float        min   [np] = {  0.50   };
  float        max   [np] = {  2.00   };
  float        step  [np] = {  0.01   };
  const char * pname [np] = { "MaQEL" };

  for(int i=0; i<np; i++) {
    LOG("MaQELFit",pDEBUG)
        << "** Setting fit param " << i
          << "(" << pname[i] << ") value = " << value[i]
             << ", range = [" << min[i] << ", " << max[i] <<"]";

    minuit->mnparm(
       i, pname[i], value[i], step[i], min[i], max[i], ierrflg);
  }

  minuit->SetFCN(FitFunc);

  // MINUIT minimization step
  ierrflg    = 0;
  arglist[0] = 500;
  arglist[1] = 1.;
  minuit->mnexcm("MIGRAD",arglist,2,ierrflg);

  // Print results
  Double_t amin,edm,errdef;
  Int_t nvpar,nparx,icstat;
  minuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  minuit->mnprin(3,amin);
}
//____________________________________________________________________________
void FitFunc (
        Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t /*iflag*/)
{
// the MINUIT fit function

  //
  // update fit parameters
  //

  AlgConfigPool * confp = AlgConfigPool::Instance();
  Registry * user_conf = confp->GlobalParameterList();
  user_conf->UnLock();

  LOG("MaQELFit", pNOTICE) << "Setting: MaQEL = " << par[0];
  user_conf->Set("QEL-Ma", par[0]);

  AlgFactory * algf = AlgFactory::Instance();
  algf->ForceReconfiguration();

  //
  // calculate chisq for the current set of fitted parameters  
  //
 
  double chisq = 0;

  // loop over all data sets included in the fit
  for(int imode=0; imode<kNDataSets; imode++) {

    // include current data set?
    vector<int>::const_iterator it = 
        find(gEnabledDataSets.begin(), gEnabledDataSets.end(), imode);
    bool skip = (it==gEnabledDataSets.end());
    if(skip) continue;

    LOG("MaQELFit", pNOTICE) << " *** Data Set : " << imode;	

    MultiGraph * mgr = gXSecDataGrph[imode];

    // loop over graphs in current data-set (one graph per experiment/publication in this data set)
    int ngr = mgr->NGraphs();
    for(int igr = 0; igr < ngr; igr++) {

       LOG("MaQELFit", pNOTICE) 
           << " Subset : " << igr+1 << "/" << ngr 
           << " --> " << imode << " / " << mgr->GetLegendEntry(igr);

       // loop over data-points in current graph
       int np = mgr->GetGraph(igr)->GetN();
       for (int ip=0; ip < np; ip++) {

          double E             = mgr->GetGraph(igr)->GetX()[ip];
          double xsec_data     = mgr->GetGraph(igr)->GetY()[ip];    // data
          double xsec_data_err = mgr->GetGraph(igr)->GetErrorY(ip); // err(data)
          double xsec_model    = GetXSecGENIE(E,imode);

          double delta = (xsec_data>0) ? (xsec_data - xsec_model) / xsec_data_err : 0.;
          chisq += delta*delta;

          LOG("MaQELFit", pNOTICE)
             << " > pnt " << ip+1 << "/" << np << " @ E = " << E << " GeV : Data = " << xsec_data
             << " +/- " << xsec_data_err << " x1E-38 cm^2, Model = " 
             << xsec_model << " x1E-38 cm^2 "
             << " >> Running chisq = " << chisq;

       } // graph points
    } // graph
  } // data set

  f = chisq;

  LOG("MaQELFit", pINFO) << "**** chisq = " << chisq;
}
//____________________________________________________________________________
void Save(string filename)
{
// post-fit write-out

  TFile fout(filename.c_str(), "recreate");
  fout.cd();

  // save fitted data-sets
  for(int imode=0; imode<kNDataSets; imode++) {
    MultiGraph * mgr = gXSecDataGrph[imode];
    for(unsigned int igr=0; igr<mgr->NGraphs(); igr++) {
      ostringstream name;
      name << "data_set_" << imode << "_" << igr;
      mgr->GetGraph(igr)->SetName(name.str().c_str());
      mgr->GetGraph(igr)->SetTitle(mgr->GetLegendEntry(igr).c_str());
      mgr->GetGraph(igr)->Write();
    }
  }

  // save best fit predictions
  for(int i=0; i<kNDataSets; i++) {
  }

  // write-out fitted params / errors  

  fout.Close();
}
//____________________________________________________________________________
