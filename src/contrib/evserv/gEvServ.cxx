//____________________________________________________________________________
/*!

\program gevserv

\brief   GENIE v+A event generation server 

         Syntax :
           gevserv [-p port]

         Options :
           [] denotes an optional argument
           -p port number (default: 9090)

\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

\created September 18, 2007

\cpright Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#include <cassert>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <TSystem.h>
#include <TServerSocket.h>
#include <TSocket.h>
#include <TMessage.h>
#include <TBits.h>
#include <TMath.h>

#include "Conventions/Units.h"
#include "EVGCore/EventRecord.h"
#include "EVGDrivers/GEVGDriver.h"
#include "EVGDrivers/GEVGPool.h"
#include "EVGDrivers/GMCJMonitor.h"
#include "GHEP/GHepFlags.h"
#include "GHEP/GHepParticle.h"
#include "GHEP/GHepRecord.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"
#include "Numerical/Spline.h"
#include "PDG/PDGCodeList.h"
#include "Utils/XSecSplineList.h"
#include "Utils/StringUtils.h"
#include "Utils/CmdLnArgParser.h"

using std::string;
using std::vector;
using std::ostringstream;

using namespace genie;
using namespace genie::utils;

// ** Prototypes
//
void GetCommandLineArgs (int argc, char ** argv);
void PrintSyntax        (void);
void RunInitChecks      (void);
void HandleMesg         (string mesg);
void Handshake          (void);
void Configure          (string mesg);
void CalcTotalXSec      (string mesg);
void GenerateEvent      (string mesg);
void Shutdown           (void);

// ** Consts & Defaults
//
const int    kDefPortNum           = 9090;  // default port number
const string kHandshakeCmdRecv     = "RUB GENIE LAMP";
const string kHandshakeMesgSent    = "YOU HAVE 3 WISHES!";
const string kConfigCmdRecv        = "CONFIG";
const string kConfigOkMesgSent     = "CONFIG COMPLETED";
const string kConfigCmdLdSpl       = "load-splines";
const string kConfigCmdNeuList     = "neutrino-list";
const string kConfigCmdTgtList     = "target-list";
const string kXSecCmdRecv          = "XSEC";
const string kXSecCmdSent          = "XSECSPL";
const string kXSecOkMesgSent       = "XSEC SENT";
const string kEvgenCmdRecv         = "EVTVTX";
const string kEvgenHdrCmdSent      = "EVTREC";
const string kEvgenStdhepCmdSent   = "STDHEP";
const string kEvgenOkMesgSent      = "EVENT GENERATED";
const string kShutdownCmdRecv      = "SHUTDOWN";
const string kShutdownOkMesgSent   = "SHUTTING DOWN";
const string kErrNoConf            = "*** NOT CONFIGURED! ***";
const string kErrNoDriver          = "*** NO EVENT GENERATION DRIVER! ***";
const string kErrNoEvent           = "*** NULL OR UNPHYSICAL EVENT! ***";
const string kErr                  = "FAILED";

// ** User-specified options:
//
int gOptPortNum;   // port number

// ** Globals
//
TSocket * gSock       = 0;      // tcp/ip socket
bool      gShutDown   = false;  // 'shutting down?' flag
bool      gConfigured = false;  // 'am I configured?' flag
GEVGPool  gGPool;               // list of GENIE event generation drivers used in job

//____________________________________________________________________________
int main(int argc, char ** argv)
{
  // Parse command line arguments
  GetCommandLineArgs(argc,argv);

  // Run some checks
  RunInitChecks();

  // Open a server socket
  TServerSocket * serv_sock = new TServerSocket(gOptPortNum, kTRUE);

  // Accept a connection
  gSock = serv_sock->Accept();

  if(!gSock) exit(1);
  LOG("gevserv", pNOTICE) << "Listening on port: " << gOptPortNum;
  
  // Set no TCP/IP NODELAY
  int delay_ok = gSock->SetOption(kNoDelay,1);
  LOG("gevserv", pNOTICE) << "TCP_NODELAY > " << delay_ok;

  // Start listening for messages & take the corresponding actions

  while(1) {

    if(gShutDown) break;

    TMessage * mesg = 0;

    gSock->Recv(mesg);

    if(!mesg) continue;
    if(mesg->What() != kMESS_STRING) continue;

    char mesg_content[2048];
    mesg->ReadString(mesg_content, 2048);

    LOG("gevserv", pNOTICE) << "Processing mesg > " << mesg_content;

    HandleMesg(mesg_content);

  } // while(1)

  return 0;
}
//____________________________________________________________________________
void HandleMesg(string mesg)
{
  if(mesg.find(kHandshakeCmdRecv.c_str()) != string::npos) 
  { 
    Handshake();
  } 
  else
  if (mesg.find(kConfigCmdRecv.c_str()) != string::npos) 
  {
    Configure(mesg);
  } 
  else
  if (mesg.find(kXSecCmdRecv.c_str()) != string::npos) 
  {
    CalcTotalXSec(mesg);
  }
  else
  if (mesg.find(kEvgenCmdRecv.c_str()) != string::npos) 
  {
    GenerateEvent(mesg);
  } 
  else
  if (mesg.find(kShutdownCmdRecv.c_str()) != string::npos) 
  {
    Shutdown();
  }
}
//____________________________________________________________________________
void Handshake(void)
{
// Reply to client messages checking whether the event server is active
//
  LOG("gevserv", pNOTICE) 
         << "GENIE was pinged by a client (lamp was rubbed)! Responding...";
  
  gSock->Send(kHandshakeMesgSent.c_str());

  LOG("gevserv", pINFO) << "...done!";
}
//____________________________________________________________________________
void Configure(string mesg)
{
// Configure the GENIE event server
// ** Load splines 
//    - if the "load-splines" command is contained in the mesg 
//    - the splines are loaded from the the XML file specified in $GSPLOAD (server-side) 
//    the "load-splines" mesg is sent
// ** Specify the neutrino list
//    - adding a "neutrino-list=<comma separated list of pdg codes>" in the mesg
// ** Specify the target list
//    - adding a "target-list=<comma separated list of pdg codes>" in the mesg
//
// Note:
// - All nedded cross section splines for the specified neutrino/target lists must 
//   be available in the specified XML file. If not GENIE will attempt building the
//   the missing ones which may increase start-up overheads.
// - You can further control GENIE (suppress modes, set messenger verbosity, ...)
//   by setting all the std GENIE env vars at the server side. 
//   See the GENIE web site.

  LOG("gevserv", pNOTICE)  << "Configuring GENIE event server";

  mesg = str::FilterString(kConfigCmdRecv, mesg); 
  mesg = str::FilterString(":", mesg); 
  mesg = str::TrimSpaces(mesg);             

  LOG("gevserv", pNOTICE) << "Configure options: " << mesg;

  // Autoload splines from the XML file pointed at the $GSPLOAD env. var.
  // (if set at the server side)
  //
  if(mesg.find(kConfigCmdLdSpl) != string::npos) {
     XSecSplineList * xspl = XSecSplineList::Instance();
     xspl->AutoLoad();

     mesg.erase(mesg.find(kConfigCmdLdSpl),12);
     mesg = str::TrimSpaces(mesg);             
  }

  // Extract neutrino and target lists from the input mesg
  //
  bool allowdup = false;
  PDGCodeList neutrinos(allowdup);
  PDGCodeList targets(allowdup);

  vector<string> conf_opt_v = str::Split(mesg," ");
  vector<string>::iterator conf_opt_iter = conf_opt_v.begin();

  for( ; conf_opt_iter != conf_opt_v.end(); ++conf_opt_iter) {
    string conf_opt = *conf_opt_iter;
    LOG("gevserv", pNOTICE) 
          << "Processing config option: " << conf_opt;

    vector<string> sv = str::Split(conf_opt, "=");
    assert(sv.size()==2);
    string list_name     = sv[0];
    string particle_list = sv[1];

    vector<string> particles = str::Split(particle_list, ",");
    vector<string>::iterator particle_iter = particles.begin();

    for( ; particle_iter != particles.end(); ++particle_iter) {
      string particle_code_str = *particle_iter;
      int particle_code = atoi(particle_code_str.c_str());

      if(list_name.find(kConfigCmdNeuList) != string::npos) 
      {
         neutrinos.push_back(particle_code);
      } else 
      if(list_name.find(kConfigCmdTgtList) != string::npos) 
      {
	 targets.push_back(particle_code);
      }
    }
  }

  LOG("gevserv", pNOTICE) 
        << "Specified neutrino list: " << neutrinos;
  LOG("gevserv", pNOTICE) 
        << "Specified target list: "   << targets;


  // Loop over the specified neutrinos and targets and for each
  // possible pair create / configure a GENIE event generation driver.
  //
  PDGCodeList::const_iterator nuiter;
  PDGCodeList::const_iterator tgtiter;

  for(nuiter = neutrinos.begin(); nuiter != neutrinos.end(); ++nuiter) {
   for(tgtiter = targets.begin(); tgtiter != targets.end(); ++tgtiter) {

     int target_code   = *tgtiter;
     int neutrino_code = *nuiter;

     InitialState init_state(target_code, neutrino_code);

     LOG("gevserv", pNOTICE)
       << "\n\n ---- Creating a GEVGDriver object configured for init-state: "
       << init_state.AsString() << " ----\n\n";

     GEVGDriver * evgdriver = new GEVGDriver;
     evgdriver->Configure(init_state);
     evgdriver->UseSplines(); // will also check if all splines needed are loaded

     gGPool.insert( GEVGPool::value_type(init_state.AsString(), evgdriver) );

   } // targets
  } // neutrinos

  LOG("gevserv", pNOTICE)
       << "All necessary GEVGDriver object were pushed into GEVGPool\n";

  gConfigured = true;

  gSock->Send(kConfigOkMesgSent.c_str());

  LOG("gevserv", pINFO) << "...done!";
}
//____________________________________________________________________________
void CalcTotalXSec(string mesg)
{
  LOG("gevserv", pNOTICE) 
       << "Sending total xsec for enabled channels  - Input info : " << mesg;

  if(!gConfigured) {
      LOG("gevserv", pERROR) 
             << "Event server is not configured - Can not generate event";
      gSock->Send(kErrNoConf.c_str());
      gSock->Send(kErr.c_str());
      return;
  }

  // Extract info from the input mesg

  mesg = str::FilterString(kXSecCmdRecv, mesg); 
  mesg = str::FilterString(":", mesg); 
  mesg = str::TrimSpaces(mesg);             

  vector<string> sv = str::Split(mesg," "); 
  assert(sv.size()==2);

  int ipdgnu  = atoi(sv[0].c_str());  // neutrino code
  int ipdgtgt = atoi(sv[1].c_str());  // target code

  // Find the appropriate event generation driver for the given initial state

  InitialState init_state(ipdgtgt, ipdgnu);
  GEVGDriver * evg_driver = gGPool.FindDriver(init_state);
  if(!evg_driver) {
     LOG("gevserv", pERROR)
       << "No GEVGDriver object for init state: " << init_state.AsString();
     gSock->Send(kErrNoDriver.c_str());
     gSock->Send(kErr.c_str());
     return;
  }

  // Ask the event generation driver to sum up the splines for all enabled
  // channels

   LOG("gevserv", pNOTICE)
       << "Requesting total cross section for init state: " 
       << init_state.AsString();

  evg_driver->CreateXSecSumSpline (
     1000 /*nknots*/, 0.001 /*Emin*/, 300 /*Emax*/, true /*in-log*/);

  const Spline * total_xsec_spl = evg_driver->XSecSumSpline();
  assert(total_xsec_spl);

  // Send back the cross section data (use 200 MeV bins from 0.010 -> 200.010 GeV)

  double dE   =   0.200;
  double Emin =   0.010;
  double Emax = 200.010;
  int    np   = (int) TMath::Ceil((Emax-Emin)/dE);

  ostringstream xsec_hdr;
  xsec_hdr <<  kXSecCmdSent << ":" << np;
  gSock->Send(xsec_hdr.str().c_str());

  for(int ip=0; ip<np; ip++) {
     double E  = Emin + ip*dE;
     double xs = TMath::Max(0., total_xsec_spl->Evaluate(E) / (1E-38*units::cm2));
     
     ostringstream xsec_spl_knot;
     xsec_spl_knot << ip << " " << E << " " << Form("%15.8e",xs);
     gSock->Send(xsec_spl_knot.str().c_str());
  }

  gSock->Send(kXSecOkMesgSent.c_str());

  LOG("gevserv", pINFO) << "...done!";
}
//____________________________________________________________________________
void GenerateEvent(string mesg)
{
  LOG("gevserv", pNOTICE) << "Generating event - Input info : " << mesg;

  if(!gConfigured) {
      LOG("gevserv", pERROR) 
             << "Event server is not configured - Can not generate event";
      gSock->Send(kErrNoConf.c_str());
      gSock->Send(kErr.c_str());
      return;
  }

  // Extract info from the input mesg

  mesg = str::FilterString(kEvgenCmdRecv, mesg); 
  mesg = str::FilterString(":", mesg); 
  mesg = str::TrimSpaces(mesg);             

  vector<string> sv = str::Split(mesg," "); 

  assert(sv.size()==11);
  int    irun        = atoi(sv[0].c_str());  // just pass through
  int    ievt        = atoi(sv[1].c_str());  // ...
  int    ipdgnunoosc = atoi(sv[2].c_str());  // ...
  double vtxx        = atof(sv[3].c_str());  // ...
  double vtxy        = atof(sv[4].c_str());  // ...
  double vtxz        = atof(sv[5].c_str());  // ...
  int    ipdgnu      = atoi(sv[6].c_str());  // neutrino code
  int    ipdgtgt     = atoi(sv[7].c_str());  // target code
  double px          = atof(sv[8].c_str());  // neutrino px
  double py          = atof(sv[9].c_str());  // neutrino py
  double pz          = atof(sv[10].c_str()); // neutrino pz
  double E           = TMath::Sqrt(px*px + py*py + pz*pz);

  TLorentzVector p4(px,py,pz,E); 
     
  // Find the appropriate event generation driver for the given initial state

  InitialState init_state(ipdgtgt, ipdgnu);
  GEVGDriver * evg_driver = gGPool.FindDriver(init_state);
  if(!evg_driver) {
     LOG("gevserv", pERROR)
       << "No GEVGDriver object for init state: " << init_state.AsString();
     gSock->Send(kErrNoDriver.c_str());
     gSock->Send(kErr.c_str());
     return;
  }

  // Generate the requested event

  EventRecord * event = evg_driver->GenerateEvent(p4);

  // Check/print the generated event
  bool failed = (event==0) || event->IsUnphysical();
  if(failed) {
      LOG("gevserv", pWARN) 
              << "Failed to generate the requested event";
      gSock->Send(kErrNoEvent.c_str());
      gSock->Send(kErr.c_str());
      return;
  }
  LOG("gevserv", pINFO) << "Generated event: " << *event;

  // Extract some summary info & convert to what MINOS expects

  const Interaction * interaction = event->Summary();
  const ProcessInfo &  proc_info  = interaction->ProcInfo();
  const Kinematics &   kine       = interaction->Kine();

  int int_type  = -1;
  if      (proc_info.IsQuasiElastic())      int_type = 1;
  else if (proc_info.IsResonant())          int_type = 2;
  else if (proc_info.IsDeepInelastic())     int_type = 3;
  else if (proc_info.IsCoherent())          int_type = 4;
  else if (proc_info.IsInverseMuDecay())    int_type = 5;
  else if (proc_info.IsNuElectronElastic()) int_type = 6;

  int iaction   = -1;
  if      (proc_info.IsWeakNC()) iaction = 0;
  else if (proc_info.IsWeakCC()) iaction = 1;
  else                           iaction = 2; // cc+nc interference

  int nucleon = -1;
  GHepParticle * hitnucl = event->HitNucleon();
  if(hitnucl) {
    nucleon = hitnucl->Pdg();
  }

  int hitquark  = interaction->InitState().Tgt().HitQrkPdg();

  bool get_selected = true;
  double xbj_sel   = kine.x (get_selected);
  double y_sel     = kine.y (get_selected);
  double W2_sel    = TMath::Power(kine.W (get_selected), 2.);
  double q2_sel    = -1 * kine.Q2(get_selected);

  double tot_xsec  = event->XSec();
  double diff_xsec = event->DiffXSec();

  int ihadmode  = 0; // need to fill

  // Send back the event through the tcp/ip socket

  ostringstream hdr1, hdr2, hdr3, hdr4, stdhep_hdr;

  hdr1 
    << kEvgenHdrCmdSent << ": "
    << irun             << " " 
    << ievt             << " " 
    << ipdgnunoosc      << " "
    << vtxx             << " " 
    << vtxy             << " " 
    << vtxz;
  hdr2 
    << ipdgnu  << " " 
    << ipdgtgt << " " 
    << px      << " " 
    << py      << " " 
    << pz;
  hdr3 
    << int_type << " " 
    << iaction  << " " 
    << nucleon  << " " 
    << hitquark << " "
    << xbj_sel  << " " 
    << y_sel    << " "  
    << W2_sel   << " " 
    << q2_sel;
  hdr4 
    << tot_xsec  << " " 
    << diff_xsec << " " 
    << ihadmode;

  stdhep_hdr 
    << kEvgenStdhepCmdSent << ": " 
    << event->GetEntriesFast();

  gSock->Send(hdr1.str().c_str());
  gSock->Send(hdr2.str().c_str());
  gSock->Send(hdr3.str().c_str());
  gSock->Send(hdr4.str().c_str());
  gSock->Send(stdhep_hdr.str().c_str());

  unsigned int i=0;
  TIter event_iter(event);
  GHepParticle * p = 0;
  while ( (p = dynamic_cast<GHepParticle *>(event_iter.Next())) ) {

      ostringstream stdhep_entry;

      stdhep_entry 
  	 << i << " " << p->Status() << " " << p->Pdg() << " "
         << p->FirstMother()   << " " << p->LastMother()   << " "
         << p->FirstDaughter() << " " << p->LastDaughter() << " "
         << p->Px() << " " << p->Py() << " " << p->Pz() << " " << p->E() << " " 
         << p->Mass() << " "
         << p->Vx() << " " << p->Vy() << " " << p->Vz() << " " << p->Vt();

      gSock->Send(stdhep_entry.str().c_str());
      i++;
  }

  // Clean-up and report success

  delete event;

  gSock->Send(kEvgenOkMesgSent.c_str());

  LOG("gevserv", pINFO) << "...done!";
}
//____________________________________________________________________________
void Shutdown(void)
{
  LOG("gevserv", pNOTICE) << "Shutting GENIE event server down ...";

  gShutDown = true;

  gSock->Send(kShutdownOkMesgSent.c_str());

  LOG("gevserv", pINFO) << "...done!";
}
//____________________________________________________________________________
void RunInitChecks(void)
{
  if(gSystem->Getenv("GSPLOAD")) {
    string splines_filename = gSystem->Getenv("GSPLOAD");
    bool is_accessible = ! (gSystem->AccessPathName( splines_filename.c_str() ));
    if (!is_accessible) {
       LOG("gevserv", pWARN) 
          << "*** The file (" << splines_filename 
          << ") specified in $GSPLOAD doesn't seem to be available!";
       LOG("gevserv", pWARN) 
          << "*** Expect a significant start-up overhead!";
    }     
  } else {
     LOG("gevserv", pWARN) << "*** $GSPLOAD was not set!";
     LOG("gevserv", pWARN) << "*** Expect a significant start-up overhead!";
  }
}
//____________________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{
  LOG("gevserv", pNOTICE) << "Parsing command line arguments";

  CmdLnArgParser parser(argc,argv);

  // port number:
  if( parser.OptionExists('p') ) {
    LOG("gevserv", pINFO) << "Reading port number";
    gOptPortNum = parser.ArgAsInt('p');
  } else {
    LOG("gevserv", pINFO)
	<< "Unspecified port number - Using default (" << kDefPortNum << ")";
    gOptPortNum = kDefPortNum;
  }
}
//____________________________________________________________________________
void PrintSyntax(void)
{
  LOG("gevserv", pNOTICE)
    << "\n\n" << "Syntax:" << "\n"
    << "   gevserv [-p port] \n";
}
//____________________________________________________________________________

