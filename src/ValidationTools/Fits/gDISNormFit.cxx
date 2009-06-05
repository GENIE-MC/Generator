//____________________________________________________________________________
/*!

\program gtune_disnorm

\brief   Cross section tuning utility: Fits DIS scale factor to high energy
         neutrino and anti-neutrino cross section data.

\syntax  gtune_disnorm -h host -u username -p password -d data_sets_to_fit

         where
          -h specifies the MySQL hostname and dbase
          -u specifies the MySQL username
          -p specifies the MySQL password
          -d specifies which data sets to fit
             (see list below, input as a comma separated list)

\example shell% export GDISABLECACHING=YES
         shell% gtune_disnorm -h mysql://localhost/NuScat -u costas -p psw -d 0,1
                      
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
#include "Conventions/Units.h"
#include "EVGDrivers/GEVGDriver.h"
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
// 0    neutrino      TOT [E>10]             
// 1    anti-neutrino TOT [E>10]            
// ..............................................

const int    kNDataSets = 2;
const double kEmin      = 10.0; // GeV
const double kEmax      = 90.0; // GeV

// keys (experiment id / measurement id) for extracting the
// cross section data from the NuValidator MySQL data-base
const char * kKeyList[kNDataSets] = {
/* 0 */ "ANL_12FT,2;ANL_12FT,4;BEBC,0;BEBC,2;BEBC,5;BEBC,8;BNL_7FT,0;BNL_7FT,4;CCFR,2;CCFRR,0;CHARM,0;CHARM,4;FNAL_15FT,1;FNAL_15FT,2;Gargamelle,0;Gargamelle,10;Gargamelle,12;IHEP_ITEP,0;IHEP_ITEP,2;IHEP_JINR,0;SKAT,0",
/* 1 */ "BEBC,1;BEBC,3;BEBC,6;BEBC,7;BNL_7FT,1;CCFR,3;CHARM,1;CHARM,5;FNAL_15FT,4;FNAL_15FT,5;Gargamelle,1;Gargamelle,11;Gargamelle,13;IHEP_ITEP,1;IHEP_ITEP,3;IHEP_JINR,1"
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
GEVGDriver *                gNuMuNdrv;
GEVGDriver *                gNuMuPdrv;
GEVGDriver *                gNuMuBarNdrv;
GEVGDriver *                gNuMuBarPdrv;

//____________________________________________________________________________
int main(int argc, char ** argv)
{
  // Initialize
  Init(argc, argv);

  // Get data from the NuVld data-base
  GetXSecData();

  // Fit DIS scale using MINUT
  DoTheFit();

  // Save plots etc
  Save("disnorm_fit.root");

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

  // Get GEVGDriver objects for calculating the total cross section 
  // without having to access the individual cross section models here
  gNuMuPdrv    = new GEVGDriver;
  gNuMuNdrv    = new GEVGDriver;
  gNuMuBarPdrv = new GEVGDriver;
  gNuMuBarNdrv = new GEVGDriver;

  gNuMuPdrv    -> Configure(kPdgNuMu,     1,1);
  gNuMuNdrv    -> Configure(kPdgNuMu,     0,1);
  gNuMuBarPdrv -> Configure(kPdgAntiNuMu, 1,1);
  gNuMuBarNdrv -> Configure(kPdgAntiNuMu, 0,1);
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
// Get the GENIE cross section, for the current values of the fitted params, 
// for the input mode at the input energy
// The easiest way to get the 'total' summing up the ~100's of processes is
// via the GEVGDriver (although it is inefficient because not all process
// cross section integrated each time were affected by the fitted params)
//
  double xsec_tot_isosc = 0;

  TLorentzVector p4(0,0,E,E);

  if(imode==0) {
      xsec_tot_isosc = 0.5 * (gNuMuNdrv->XSecSum(p4) + 
                              gNuMuPdrv->XSecSum(p4));
  } else
  if(imode==1) {
      xsec_tot_isosc = 0.5 * (gNuMuBarNdrv->XSecSum(p4) + 
                              gNuMuBarPdrv->XSecSum(p4));
  }

  xsec_tot_isosc /= (1E-38 * units::cm2);
  return xsec_tot_isosc;
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
    LOG("DISNormFit",pDEBUG)
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

  LOG("DISNormFit", pNOTICE) << "Setting: DIS-XSecScale = " << par[0];
  user_conf->Set("DIS-XSecScale", par[0]);

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

    LOG("DISNormFit", pNOTICE) << " *** Data Set : " << imode;	

    MultiGraph * mgr = gXSecDataGrph[imode];

    // loop over graphs in current data-set (one graph per experiment/publication in this data set)
    int ngr = mgr->NGraphs();
    for(int igr = 0; igr < ngr; igr++) {

       LOG("DISNormFit", pNOTICE) 
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

          LOG("DISNormFit", pNOTICE)
             << " > pnt " << ip+1 << "/" << np << " @ E = " << E << " GeV : Data = " << xsec_data
             << " +/- " << xsec_data_err << " x1E-38 cm^2, Model = " 
             << xsec_model << " x1E-38 cm^2 "
             << " >> Running chisq = " << chisq;

       } // graph points
    } // graph
  } // data set

  f = chisq;

  LOG("DISNormFit", pINFO) << "**** chisq = " << chisq;
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
