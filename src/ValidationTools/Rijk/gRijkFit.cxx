//____________________________________________________________________________
/*!

\program gtune_rijk

\brief   Cross section tuning utility: Fits the Rijk parameters controlling
         the cross section in the RES/DIS transition region.
         Evolved from an older NuValidator app. that was working with neugen3.

\syntax  gtune_rijk -h host -u username -p password

\example gtune_rijk -h mysql://localhost/NuScat -u costas -p mypass1
                      
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
#include <sstream>
#include <string>

#include <TFile.h>
#include <TMinuit.h>
#include <TMath.h>
#include <TSQLServer.h>
#include <TH1.h>
#include <TFile.h>

#include "Algorithm/AlgFactory.h"
#include "Algorithm/AlgConfigPool.h"
#include "Base/XSecAlgorithmI.h"
#include "Conventions/Units.h"
#include "Messenger/Messenger.h"
#include "Registry/Registry.h"
#include "ValidationTools/NuVld/DBStatus.h"
#include "ValidationTools/NuVld/DBI.h"
#include "ValidationTools/NuVld/DBTable.h"
#include "ValidationTools/NuVld/DBQueryString.h"
#include "ValidationTools/NuVld/MultiGraph.h"

using std::string;
using std::ostringstream;

using namespace genie;
using namespace genie::nuvld;

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
void          Init         (int argc, char ** argv);
string        GetArgument  (int argc, char ** argv, const char * option);
double        GetXSecGENIE (double E, int imode);
void          GetXSecData  (void);
DBQueryString FormQuery    (const char * key_list);
void          DoTheFit     (void);
void          FitFunc      (Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t iflag);
void          Save         (string filename);

// globals
//
DBI * gDBI = 0; ///< NuValidator data base interface
DBTable<DBNuXSecTableRow> * gXSecData [kNDataSets]; ///< fitted data
MultiGraph *  gXSecDataGrph [kNDataSets]; ///< fitted data as graphs

//____________________________________________________________________________
int main(int argc, char ** argv)
{
  // Initialize
  Init(argc, argv);

  // Get data from the NuVld data-base
  GetXSecData();

  // Fit Rijk using MINUT
  DoTheFit();

  // Save plots etc for post-fic inspection
  Save("Rijk.root");

  return 0;
}
//____________________________________________________________________________
void Init(int argc, char ** argv)
{
  // get command line arguments
  string url      = GetArgument(argc, argv, "-h");
  string username = GetArgument(argc, argv, "-u");
  string passwd   = GetArgument(argc, argv, "-p");

  // establish connection with the NuValidator data-base and create a
  // data-base interface
  TSQLServer * sql_server = TSQLServer::Connect(
       url.c_str(), username.c_str(), passwd.c_str());
  assert(sql_server && sql_server->IsConnected());
  gDBI = new DBI(sql_server);
  assert(gDBI);

  // Get the cross section models that we will be using
  AlgFactory * algf = AlgFactory::Instance();

//  alg = dynamic_cast<const XSecAlgorithmI *>(algf->GetAlgorithm("QELPXSec","Default"));



  // Set all fitted GENIE Rijk parameters to 1 & reconfigure
  AlgConfigPool * confp = AlgConfigPool::Instance();
  Registry * user_conf = confp->GlobalParameterList();
  user_conf->UnLock();

  user_conf->Set("DIS-HMultWgt-vp-CC-m2", 1.);
  user_conf->Set("DIS-HMultWgt-vp-CC-m3", 1.);
  user_conf->Set("DIS-HMultWgt-vn-CC-m2", 1.);
  user_conf->Set("DIS-HMultWgt-vn-CC-m3", 1.);

  algf->ForceReconfiguration();
/*
  // Create all the interaction objects that should be used for driving the
  // cross section algorithms
  Interaction * res_cc_p = Interaction::RESCC();
  Interaction * res_cc_n = Interaction::RESCC();
  Interaction * dis_cc_p = Interaction::RESCC();
  Interaction * dis_cc_n = Interaction::RESCC();
*/
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

  double xsec = 0.;

  switch(imode) {

  // 0 -> total xsec
  case(0): 
  {
    break;
  }

  // 1 -> v p -> l- p pi+ xsec 
  case(1):
  {
    break;
  }

  // 2 -> v n -> l- p pi0 xsec
  case(2):
  {
    break;
  }

  // 3 -> v n -> l- n pi+ xsec
  case(3):
  {
    break;
  }

  // 4 -> v n -> l- p pi+ pi- xsec
  case(4):
  {
    break;
  }

  // 5 -> v p -> l- p pi+ pi0 xsec
  case(5):
  {
    break;
  }

  // 6 -> v p -> l- n pi+ pi+ xsec
  case(6):
  {
    break;
  }
  default:
    xsec = 0.;
  }

  xsec /= (1E-38 * units::cm2);

  LOG("RijkFit", pNOTICE) 
    << "xsec(E = " << E << " GeV, imode = " << imode << ") = " 
    << xsec << " 1E-38 * cm2";

  return xsec;
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

  float        value [np] = { 1.00,          1.00,           1.00,           1.00           };
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

    MultiGraph * mgr = gXSecDataGrph[imode];

    // loop over graphs in current data-set
    int ngr = mgr->NGraphs();
    for(int igr = 0; igr < ngr; igr++) {

       // loop over data-points in current graph
       int np = mgr->GetGraph(igr)->GetN();
       for (int ip=0; ip < np; ip++) {

          double E             = mgr->GetGraph(igr)->GetX()[ip];
          double xsec_data     = mgr->GetGraph(igr)->GetY()[ip];    // data
          double xsec_data_err = mgr->GetGraph(igr)->GetErrorY(ip); // err(data)
          double xsec_model    = GetXSecGENIE(E,imode);

          double delta = (xsec_data>0) ? (xsec_data - xsec_model) / xsec_data_err : 0.;
          chisq += delta*delta;

/*
          LOG("RijkFit",pINFO)
             << "    E = " << E << " -> Data = " << XSecD
                    << " +/- " << dXSecD << ", Model = " << XSecM
                                     << "/ Running chisq = " << chisq;
*/
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
  for(int i=0; i<kNDataSets; i++) {
  }

  // save best fit predictions
  for(int i=0; i<kNDataSets; i++) {
  }

  // write-out fitted params / errors  

  fout.Close();
}
//____________________________________________________________________________
