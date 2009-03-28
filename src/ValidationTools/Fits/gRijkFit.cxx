//____________________________________________________________________________
/*!

\program gtune_rijk

\brief   Cross section tuning utility: Fits parameters controlling the cross
         section in the RES/DIS transition region.

         Not a user utility:
         Running and interpreting results is an expert operation.

\syntax  gtune_rijk -h host -u username -p password -d data_sets_to_fit

         where
          -h specifies the MySQL hostname and dbase
          -u specifies the MySQL username
          -p specifies the MySQL password
          -d specifies which data sets to fit
             (see list below, input as a comma separated list)

\example shell% export GDISABLECACHING=YES
         shell% gtune_rijk -h mysql://localhost/NuScat -u costas -p pswd -d 0,1,2,3
                      
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
#include <TH1.h>
#include <TFile.h>
#include <TLorentzVector.h>

#include "Algorithm/AlgFactory.h"
#include "Algorithm/AlgConfigPool.h"
#include "Base/XSecAlgorithmI.h"
#include "BaryonResonance/BaryonResList.h"
#include "Conventions/Units.h"
#include "EVGDrivers/GEVGDriver.h"
#include "Interaction/SppChannel.h"
#include "Messenger/Messenger.h"
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
// ...............................
// Fitted data sets :
// -------------------------------
// 0 -> total xsec
// 1 -> v p -> l- p pi+ xsec 
// 2 -> v n -> l- p pi0 xsec
// 3 -> v n -> l- n pi+ xsec
// 4 -> v n -> l- p pi+ pi- xsec
// 5 -> v p -> l- p pi+ pi0 xsec
// 6 -> v p -> l- n pi+ pi+ xsec
// ...............................

const int    kNDataSets = 7;
const double kEmin      =  0.1; // GeV
const double kEmax      = 20.0; // GeV

// keys (experiment id / measurement id) for extracting the
// cross section data from the NuValidator MySQL data-base
const char * kKeyList[kNDataSets] = {
/* 0 */ "ANL_12FT,2;ANL_12FT,4;BEBC,0;BEBC,2;BEBC,5;BEBC,8;BNL_7FT,0;BNL_7FT,4;CCFR,2;CCFRR,0;CHARM,0;CHARM,4;FNAL_15FT,1;FNAL_15FT,2;Gargamelle,0;Gargamelle,10;Gargamelle,12;IHEP_ITEP,0;IHEP_ITEP,2;IHEP_JINR,0;SKAT,0",
/* 1 */ "ANL_12FT,0;ANL_12FT,5;ANL_12FT,8;BEBC,4;BEBC,9;BEBC,13;BNL_7FT,5;FNAL_15FT,0;Gargamelle,4;SKAT,4;SKAT,5",
/* 2 */ "ANL_12FT,6;ANL_12FT,9;BNL_7FT,6;SKAT,6",
/* 3 */ "ANL_12FT,7;ANL_12FT,10;BNL_7FT,7;SKAT,7",
/* 4 */ "ANL_12FT,11;BNL_7FT,8",
/* 5 */ "ANL_12FT,12",
/* 6 */ "ANL_12FT,13"
};

// func prototypes
//
void          Init           (int argc, char ** argv);
string        GetArgument    (int argc, char ** argv, const char * option);
double        GetXSecGENIE   (double E, int imode);
void          GetXSecData    (void);
Interaction * GetExclusive   (int imode, bool res);
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
BaryonResList               gResList;
const XSecAlgorithmI *      gRESXSecModel;
const XSecAlgorithmI *      gDISXSecModel;
GEVGDriver *                gNuMuPdrv;
GEVGDriver *                gNuMuNdrv;

//____________________________________________________________________________
int main(int argc, char ** argv)
{
  // Initialize
  Init(argc, argv);

  // Get data from the NuVld data-base
  GetXSecData();

  // Fit Rijk using MINUT
  DoTheFit();

  // Save plots etc 
  Save("rijk_fit.root");

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

  gRESXSecModel = 
     dynamic_cast<const XSecAlgorithmI *>(
             algf->GetAlgorithm("genie::ReinSeghalRESPXSec","Default"));
  gDISXSecModel = 
     dynamic_cast<const XSecAlgorithmI *>(
             algf->GetAlgorithm("genie::DISPartonModelPXSec","CC-Default"));

  LOG("RijkFit", pNOTICE) << "Got algorithm: " << gRESXSecModel->Id();
  LOG("RijkFit", pNOTICE) << "Got algorithm: " << gDISXSecModel->Id();

  // Set all fitted GENIE Rijk parameters to 1 & reconfigure
  AlgConfigPool * confp = AlgConfigPool::Instance();
  Registry * user_conf = confp->GlobalParameterList();
  user_conf->UnLock();

  user_conf->Set("DIS-HMultWgt-vp-CC-m2", 1.);
  user_conf->Set("DIS-HMultWgt-vp-CC-m3", 1.);
  user_conf->Set("DIS-HMultWgt-vn-CC-m2", 1.);
  user_conf->Set("DIS-HMultWgt-vn-CC-m3", 1.);

  algf->ForceReconfiguration();

  // List of resonances to take into account
  string resonanes = user_conf->GetString("ResonanceNameList");
  gResList.DecodeFromNameList(resonanes);

  //
  gNuMuPdrv = new GEVGDriver;
  gNuMuNdrv = new GEVGDriver;
  gNuMuPdrv->Configure(kPdgNuMu,1,1);
  gNuMuNdrv->Configure(kPdgNuMu,0,1);
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
// get the GENIE cross section, for the current Rijk parameters, for the input
// mode (matching the fitted data ids; see header) at the input energy

  //
  // 0 -> total xsec
  //
  if(imode == 0) {
    TLorentzVector p4(0,0,E,E);
    double xsec_tot = 0.5 * (gNuMuPdrv->XSecSum(p4) + gNuMuNdrv->XSecSum(p4));
    xsec_tot /= (1E-38 * units::cm2);

    LOG("RijkFit", pDEBUG) 
       << "xsec(total; E = " << E << " GeV) = " << xsec_tot << " 1E-38 * cm2";

    return xsec_tot;
  }


  //
  // single pion channels
  //
  // 1 -> v p -> l- p pi+ xsec 
  // 2 -> v n -> l- p pi0 xsec
  // 3 -> v n -> l- n pi+ xsec
  //
  else if(imode == 1 || imode == 2 || imode == 3) 
  {
     double xsec_1pi     = 0;
     double xsec_1pi_res = 0;
     double xsec_1pi_dis = 0;
 
     Interaction * res = GetExclusive(imode, true );
     Interaction * dis = GetExclusive(imode, false);

     SppChannel_t spp_channel = kSppNull;
     if      (imode == 1)  spp_channel = kSpp_vp_cc_10100; 
     else if (imode == 2)  spp_channel = kSpp_vn_cc_10010; 
     else if (imode == 3)  spp_channel = kSpp_vn_cc_01100; 

     // calculate the resonance singe pi contribution
     //
     unsigned int nres = gResList.NResonances();
     for(unsigned int ires = 0; ires < nres; ires++) 
     {
        Resonance_t resonance_id = gResList.ResonanceId(ires); 

        res -> ExclTagPtr()   -> SetResonance(resonance_id); 
        res -> InitStatePtr() -> SetProbeE(E);

        // xsec for exciting this resonance
        double xsec_res = gRESXSecModel->Integral(res);

        // resonance contribution to the 1pi channel
        double br  = SppChannel::BranchingRatio(spp_channel, resonance_id);  // branching ratio  
        double igg = SppChannel::IsospinWeight (spp_channel, resonance_id);  // isospin Glebsch-Gordon coeff.

        xsec_1pi_res += (xsec_res * br * igg);
     }

     // calculate the DIS singe pi contribution
     //
     dis -> InitStatePtr() -> SetProbeE(E);
     xsec_1pi_dis = gDISXSecModel->Integral(dis);

     // sum-up
     xsec_1pi_res /= (1E-38 * units::cm2);
     xsec_1pi_dis /= (1E-38 * units::cm2);
     xsec_1pi = xsec_1pi_res + xsec_1pi_dis;

     LOG("RijkFit", pDEBUG) 
        << "xsec(1pi / mode = " << imode << ", E = " << E << " GeV) = " 
        << xsec_1pi << " 1E-38 * cm2 = (" << xsec_1pi_dis << " [dis] + " << xsec_1pi_res << " [res])";

     delete dis;
     delete res;

     return xsec_1pi;

  } // modes 1,2,3


  //
  // multi pion channels
  //
  // 4 -> v n -> l- p pi+ pi- xsec
  // 5 -> v p -> l- p pi+ pi0 xsec
  // 6 -> v p -> l- n pi+ pi+ xsec
  //

  else if(imode == 4 || imode == 5 || imode == 6) {

     double xsec_2pi     = 0;
     double xsec_2pi_res = 0;
     double xsec_2pi_dis = 0;

     Interaction * res = GetExclusive(imode, true );
     Interaction * dis = GetExclusive(imode, false);

     // calculate the RES 2pi contribution
     // ...

     // calculate the DIS 2pi contribution
     dis -> InitStatePtr() -> SetProbeE(E);
     xsec_2pi_dis = gDISXSecModel->Integral(dis);

     // sum-up
     xsec_2pi_res /= (1E-38 * units::cm2);
     xsec_2pi_dis /= (1E-38 * units::cm2);
     xsec_2pi = xsec_2pi_res + xsec_2pi_dis;

     LOG("RijkFit", pDEBUG) 
        << "xsec(2pi / mode = " << imode << ", E = " << E << " GeV) = " 
        << xsec_2pi << " 1E-38 * cm2 = (" << xsec_2pi_dis << " [dis] + " << xsec_2pi_res << " [res])";

     delete dis;
     delete res;

     return xsec_2pi;

  } // modes 4,5,6

  return 0;
}
//____________________________________________________________________________
Interaction * GetExclusive (int imode, bool res)
{
  Interaction * in = 0;

  // 1 -> v p -> l- p pi+ xsec 
  if(imode==1) {
    if(res) { in = Interaction::RESCC(kPdgTgtFreeP, kPdgProton,  kPdgNuMu); }
    else    { in = Interaction::DISCC(kPdgTgtFreeP, kPdgProton,  kPdgNuMu); }
    in->ExclTagPtr()->SetNPions(1,0,0);
    in->ExclTagPtr()->SetNNucleons(1,0);
  }

  // 2 -> v n -> l- p pi0 xsec
  else if(imode==2) {
    if(res) { in = Interaction::RESCC(kPdgTgtFreeN, kPdgNeutron,  kPdgNuMu); }
    else    { in = Interaction::DISCC(kPdgTgtFreeN, kPdgNeutron,  kPdgNuMu); }
    in->ExclTagPtr()->SetNPions(0,1,0);
    in->ExclTagPtr()->SetNNucleons(1,0);
  }

  // 3 -> v n -> l- n pi+ xsec
  else if(imode==3) {
    if(res) { in = Interaction::RESCC(kPdgTgtFreeN, kPdgNeutron,  kPdgNuMu); }
    else    { in = Interaction::DISCC(kPdgTgtFreeN, kPdgNeutron,  kPdgNuMu); }
    in->ExclTagPtr()->SetNPions(1,0,0);
    in->ExclTagPtr()->SetNNucleons(0,1);
  }

  // 4 -> v n -> l- p pi+ pi- xsec
  else if(imode==4) {
    if(res) { in = Interaction::RESCC(kPdgTgtFreeN, kPdgNeutron,  kPdgNuMu); }
    else    { in = Interaction::DISCC(kPdgTgtFreeN, kPdgNeutron,  kPdgNuMu); }
    in->ExclTagPtr()->SetNPions(1,0,1);
    in->ExclTagPtr()->SetNNucleons(1,0);
  }

  // 5 -> v p -> l- p pi+ pi0 xsec
  else if(imode==5) {
    if(res) { in = Interaction::RESCC(kPdgTgtFreeP, kPdgProton,  kPdgNuMu); }
    else    { in = Interaction::DISCC(kPdgTgtFreeP, kPdgProton,  kPdgNuMu); }
    in->ExclTagPtr()->SetNPions(1,1,0);
    in->ExclTagPtr()->SetNNucleons(1,0);
  }

  // 6 -> v p -> l- n pi+ pi+ xsec
  else if(imode==6) {
    if(res) { in = Interaction::RESCC(kPdgTgtFreeP, kPdgProton,  kPdgNuMu); }
    else    { in = Interaction::DISCC(kPdgTgtFreeP, kPdgProton,  kPdgNuMu); }
    in->ExclTagPtr()->SetNPions(2,0,0);
    in->ExclTagPtr()->SetNNucleons(0,1);
  }

  return in;  
}
//____________________________________________________________________________
void DoTheFit(void)
{
  // Initialize MINUIT
  //
  const int np = 4;
  TMinuit * minuit = new TMinuit(np);

  double arglist[10];
  int ierrflg = 0;

  arglist[0] = 1;
  minuit->mnexcm("SET ERR",arglist,1,ierrflg);

  float        value [np] = { 0.10,          1.00,           0.30,           1.00           };
  float        min   [np] = { 0.00,          0.00,           0.00,           0.00           };
  float        max   [np] = { 2.00,          2.00,           2.00,           2.00           };
  float        step  [np] = { 0.05,          0.05,           0.05,           0.05           };
  const char * pname [np] = {"R(vp/CC/1pi)", "R(vp/CC/2pi)", "R(vn/CC/1pi)", "R(vn/CC/2pi)" };

  for(int i=0; i<np; i++) {
    LOG("RijkFit",pDEBUG)
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

  LOG("RijkFit", pNOTICE) << "Setting: R(vp/CC/1pi) = " << par[0];
  LOG("RijkFit", pNOTICE) << "Setting: R(vp/CC/2pi) = " << par[1];
  LOG("RijkFit", pNOTICE) << "Setting: R(vn/CC/1pi) = " << par[2];
  LOG("RijkFit", pNOTICE) << "Setting: R(vn/CC/2pi) = " << par[3];

  user_conf->Set("DIS-HMultWgt-vp-CC-m2", par[0]);
  user_conf->Set("DIS-HMultWgt-vp-CC-m3", par[1]);
  user_conf->Set("DIS-HMultWgt-vn-CC-m2", par[2]);
  user_conf->Set("DIS-HMultWgt-vn-CC-m3", par[3]);

  AlgFactory * algf = AlgFactory::Instance();
  algf->ForceReconfiguration();

  //
  // calculate a global chisq for the current set of prediction
  //
 
  double chisq = 0;

  // loop over all data sets included in the fit
  for(int imode=0; imode<kNDataSets; imode++) {

    // include current data set?
    vector<int>::const_iterator it =
        find(gEnabledDataSets.begin(), gEnabledDataSets.end(), imode);
    bool skip = (it==gEnabledDataSets.end());
    if(skip) continue;

    LOG("RijkFit",pNOTICE) << " *** Data Set : " << imode;	

    MultiGraph * mgr = gXSecDataGrph[imode];

    // loop over graphs in current data-set (one graph per experiment/publication in this data set)
    int ngr = mgr->NGraphs();
    for(int igr = 0; igr < ngr; igr++) {

       LOG("RijkFit",pNOTICE) 
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

          LOG("RijkFit", pNOTICE)
             << " > pnt " << ip+1 << "/" << np << " @ E = " << E << " GeV : Data = " << xsec_data
             << " +/- " << xsec_data_err << " x1E-38 cm^2, Model = " 
             << xsec_model << " x1E-38 cm^2 "
             << " >> Running chisq = " << chisq;

       } // graph points
    } // graph
  } // data set

  f = chisq;

  LOG("NuVldFit", pINFO) << "**** chisq = " << chisq;
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
