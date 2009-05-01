//____________________________________________________________________________
/*!

\program gvld_cohtest

\brief   

\syntax  gvld_cohtest -h host -u username -p password -d sets_to_fit

         where 
	  -h specifies the MySQL hostname and dbase
	  -u specifies the MySQL username
	  -p specifies the MySQL password
	  -d specifies which data sets to fit 
             (see list below, input as a comma separated list)

\example gtune_mafit -h mysql://localhost/NuScat -u costas -p mypass1 -d 0,3
                      
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
#include <TCanvas.h>
#include <TMinuit.h>
#include <TMath.h>
#include <TSQLServer.h>
#include <TGraph.h>
#include <TGraphErrors.h>

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
// 0    numu    NC COH pi / A = 20 
// 1    numu    CC COH pi / A = 20 
// 2    numubar CC COH pi / A = 20 
// 3    numu    NC COH pi / A = 27
// 4    numu    NC COH pi / A = 30 
// 5    numu    CC COH pi / A = 30 
// 6    numubar CC COH pi / A = 30 
// ..............................................

const int    kNDataSets = 7;
const double kEmin      =   0.1; // GeV
const double kEmax      = 200.0; // GeV

// keys (experiment id / measurement id) for extracting the
// cross section data from the NuValidator MySQL data-base
const char * kKeyList[kNDataSets] = {
/* 0 */ "CHARM,2",
/* 1 */ "BEBC,11;CHARM,6;FNAL_15FT,8",
/* 2 */ "BEBC,10;CHARM,7;FNAL_15FT,7",
/* 3 */ "AachenPadova,0",
/* 4 */ "Gargamelle,14;SKAT,3",
/* 5 */ "SKAT,1",
/* 6 */ "SKAT,2"
};
double kA[kNDataSets] = {
/* 0 */ 20,
/* 1 */ 20,
/* 2 */ 20,
/* 3 */ 27,
/* 4 */ 30,
/* 5 */ 30,
/* 6 */ 30
};
int kNeu[kNDataSets] = {
/* 0 */ kPdgNuMu,
/* 1 */ kPdgNuMu,
/* 2 */ kPdgAntiNuMu,
/* 3 */ kPdgNuMu,
/* 4 */ kPdgNuMu,
/* 5 */ kPdgNuMu,
/* 6 */ kPdgAntiNuMu
};
int kTgt[kNDataSets] = {
/* 0 */ 1000100200,
/* 1 */ 1000100200,
/* 2 */ 1000100200,
/* 3 */ 1000130270,
/* 4 */ 1000150300,
/* 5 */ 1000150300,
/* 6 */ 1000150300
};
InteractionType_t kCurr[kNDataSets] = {
/* 0 */ kIntWeakNC,
/* 1 */ kIntWeakCC,
/* 2 */ kIntWeakCC,
/* 3 */ kIntWeakNC,
/* 4 */ kIntWeakNC,
/* 5 */ kIntWeakCC,
/* 6 */ kIntWeakCC
};

// func prototypes
//
void          Init           (int argc, char ** argv);
string        GetArgument    (int argc, char ** argv, const char * option);
void          GetXSecGENIE   (void);
void          GetXSecData    (void);
DBQueryString FormQuery      (const char * key_list);
void          Save           (string filename);

// globals
//
vector<int>                 gEnabledDataSets;
DBI *                       gDBI = 0;                   ///< dbase interface
DBTable<DBNuXSecTableRow> * gXSecData     [kNDataSets]; ///< 
MultiGraph *                gXSecDataGrph [kNDataSets]; ///< 
const XSecAlgorithmI *      gCOHXSecModel;
TGraph *                    gXSecGENIE    [kNDataSets];

//____________________________________________________________________________
int main(int argc, char ** argv)
{
  // Initialize
  Init(argc, argv);

  // Get data from the NuVld data-base
  GetXSecData();

  // Get data from the NuVld data-base
  GetXSecGENIE();

  // Save plots etc
  Save("coh.root");

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
  gCOHXSecModel = 
     dynamic_cast<const XSecAlgorithmI *>(
        algf->GetAlgorithm("genie::ReinSeghalCOHPiPXSec","Default"));

  LOG("CohTest", pNOTICE) 
     << "Got algorithm: " << gCOHXSecModel->Id();
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
void GetXSecGENIE(void)
{
  const int np = 30;
  double dE = (kEmax-kEmin)/(np-1);

  for(int imode=0; imode<kNDataSets; imode++) {

     double E    [np];
     double xsec [np];

     Interaction * coh = 0;
     if(kCurr[imode] == kIntWeakCC) { coh = Interaction::COHCC(kTgt[imode], kNeu[imode]); }
     else                           { coh = Interaction::COHNC(kTgt[imode], kNeu[imode]); }

     for(int ie=0; ie<np; ie++) {
        E[ie] = kEmin + ie*dE;
        coh->InitStatePtr()->SetProbeE(E[ie]);
        xsec[ie] = gCOHXSecModel->Integral(coh);
        xsec[ie] /= (1E-40*units::cm2);
     }
     gXSecGENIE[imode] = new TGraph(np,E,xsec);
     delete coh;
  }
}
//____________________________________________________________________________
void Save(string filename)
{
  TFile fout(filename.c_str(), "recreate");
  fout.cd();

  // save data-sets
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

  // save corersponding predictions
  for(int imode=0; imode<kNDataSets; imode++) {
      ostringstream name;
      name << "genie_xsec_" << imode;
      gXSecGENIE[imode]->SetName(name.str().c_str());
      gXSecGENIE[imode]->Write();
  }

  fout.Close();
}
//____________________________________________________________________________
