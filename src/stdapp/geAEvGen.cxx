//____________________________________________________________________________
/*!

\program geAevgen **** in devel *****

\brief   Generate inclusive e+A scattering events.

         Syntax :
           geAevgen -n nevents -t theta -e energy -ep energy' ... 

         Options:
           [] Denotes an optional argument
            ...
            ...
 
\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory

\created June 14, 2007

\cpright Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#include <cassert>
#include <sstream>
#include <string>

#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>

#include "Conventions/Constants.h"
#include "EVGCore/EventRecord.h"
#include "GHEP/GHepParticle.h"
#include "GHEP/GHepStatus.h"
#include "Interaction/Interaction.h"
#include "Ntuple/NtpMCFormat.h"
#include "Ntuple/NtpMCTreeHeader.h"
#include "Ntuple/NtpMCEventRecord.h"
#include "Ntuple/NtpWriter.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGLibrary.h"
#include "PDG/PDGUtils.h"
#include "PDG/PDGCodes.h"
#include "Utils/CmdLineArgParserUtils.h"
#include "Utils/CmdLineArgParserException.h"

using std::ostringstream;
using std::string;

using namespace genie;
using namespace genie::constants;

// function prototypes
void   GetCommandLineArgs  (int argc, char ** argv);
void   PrintSyntax         (void);
void   CreateEventTemplate (void);

// command-line arguments
Long64_t gOptN;         // number of events to generate
int      gOptTarget;    //
double   gOptE;         // incoming e- energy E
double   gOptEp;        // outgoing e- energy E'
double   gOptTheta;     // scattering angle
int      gOptRunNu = 0; // run num.

EventRecord gEventTemplate;

//_________________________________________________________________________________
int main(int argc, char ** argv)
{
  //-- scan the command line arguments 
  GetCommandLineArgs(argc,argv);

  //-- initialize an Ntuple Writer
  NtpWriter ntpw(kNFGHEP, gOptRunNu);
  ntpw.Initialize("gntp-eA");

  CreateEventTemplate();
  LOG("geAevgen", pNOTICE) << "Output event template: "<< gEventTemplate;

  for(Long64_t iev = 0; iev < gOptN; iev++) {
    LOG("geAevgen", pNOTICE) << "Generating event: ....... " << iev;
 
    // init
    EventRecord event(gEventTemplate);


    // save eN event to the new event tree
    LOG("geAevgen", pNOTICE) << "Output eA event: "<< event;
    ntpw.AddEventRecord(iev++, &event);
  }

  //-- save event tree
  ntpw.Save();

  LOG("geAevgen", pINFO)  << "Done!";
  return 0;
}
//_________________________________________________________________________________
void CreateEventTemplate (void)
{
  PDGLibrary * pdglib = PDGLibrary::Instance();

  double m  = pdglib->Find(kPdgElectron)->Mass();
  double M  = pdglib->Find(gOptTarget)->Mass();
  double E  = gOptE;
  double p  = TMath::Sqrt( TMath::Max(0., E*E-m*m) );
  double Ep = gOptEp;
  double pp = TMath::Sqrt( TMath::Max(0., Ep*Ep-m*m) );
  double px = p-pp;
  double Ex = E+M-Ep;

  GHepParticle e  (kPdgElectron,     kIStInitialState,             -1,-1,-1,-1, 0,0,p, E,  0,0,0,0);
  GHepParticle A  (gOptTarget,       kIStInitialState,             -1,-1,-1,-1, 0,0,0, M,  0,0,0,0);
  GHepParticle ep (kPdgElectron,     kIStStableFinalState,          0,-1,-1,-1, 0,0,pp,Ep, 0,0,0,0);
  GHepParticle X  (kPdgHadronicBlob, kIStFinalStateNuclearRemnant,  1,-1,-1,-1, 0,0,px,Ex, 0,0,0,0);

  gEventTemplate.AddParticle(e);
  gEventTemplate.AddParticle(A);
  gEventTemplate.AddParticle(ep);
  gEventTemplate.AddParticle(X);
}
//_________________________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{
  LOG("geAevgen", pNOTICE) << "*** Parsing commad line arguments";

gOptN      = 10;         
gOptTarget = 1000260560;    
gOptE      = 1;        
gOptEp     = 0.8;        
gOptTheta  = 30;     
gOptRunNu  = 0; 

/*
  // number of events
  try {
    LOG("geAevgen", pINFO) << "Reading number of events to generate";
    gOptN = genie::utils::clap::CmdLineArgAsInt(argc,argv,'n');
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("geAevgen", pFATAL)
        << "Unspecified number of events to generate";
      PrintSyntax();
    }
  }
*/
}
//_________________________________________________________________________________
void PrintSyntax(void)
{

}
//_________________________________________________________________________________
