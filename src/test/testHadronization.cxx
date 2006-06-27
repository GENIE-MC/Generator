//____________________________________________________________________________
/*!

\program testHadronization

\brief   test program used for testing/debugging the KNO & PYTHIA hadronizers

        Syntax :
           testHadronization -n nevents -t test -a hadronizer -c config [-q]

         Options :
           -n  number of events
           -t  test type 
                 0: multiplicities
                 1: multiplicities + phase space decay
           -a  hadronizer (algorithm name, eg genie::KNOHadrinization)
           -c  hadronizer config set
           -q  set hit quark (needed for PYTHIA, not needed for KNO)

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created June 20, 2004

\cpright Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
         All rights reserved.
         For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#include <sstream>
#include <cassert>
#include <cstdlib>
#include <string>
#include <vector>

#include <TFile.h>
#include <TDirectory.h>
#include <TTree.h>
#include <TH1F.h>
#include <TLorentzVector.h>
#include <TMCParticle6.h>
#include <TMath.h>
#include <TClonesArray.h>
#include <TIterator.h>

#include "Algorithm/Algorithm.h"
#include "Algorithm/AlgFactory.h"
#include "Fragmentation/HadronizationModelI.h"
#include "Interaction/ProcessInfo.h"
#include "Interaction/InitialState.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "Utils/FragmRecUtils.h"
#include "Utils/CmdLineArgParserUtils.h"
#include "Utils/CmdLineArgParserException.h"

using std::string;
using std::vector;
using std::endl;
using std::ostringstream;

using namespace genie;

void PrintSyntax        (void);
void GetCommandLineArgs (int argc, char ** argv);

void testMultiplicities   (int n, const HadronizationModelI * m);
void testPhaseSpaceDecayer(int n, const HadronizationModelI * m);

void FillQrkArray(InteractionType_t it, int nu, 
             int * QrkCode, bool * SeaQrk, int nmax, int & nqrk);

TFile * gOutFile   = 0;
int     gNEvents   = -1;
int     gTestId    = -1;
bool    gSetHitQrk = false;
string  gHadAlg    = "";
string  gHadConfig = "";
//____________________________________________________________________________
int main(int argc, char ** argv)
{
  GetCommandLineArgs(argc, argv);

  AlgFactory * algf = AlgFactory::Instance();

  const HadronizationModelI * model = 
          dynamic_cast<const HadronizationModelI *> (
                                  algf->GetAlgorithm(gHadAlg, gHadConfig));
  assert(model);

  gOutFile = new TFile("./genie-hadronization.root","recreate");

  //-- test hadronization model model multiplicities
  if(gTestId==0) testMultiplicities(gNEvents, model);

  //-- test hadronization model (full events)
  if(gTestId==1) testPhaseSpaceDecayer(gNEvents, model);

  gOutFile->Close();

  return 0;
}
//____________________________________________________________________________
void testMultiplicities(int nevents, const HadronizationModelI * model)
{
  const int kNNu     = 2;
  const int kNNuc    = 2;
  const int kNQrkMax = 4;
  const int kNW      = 8;

  int    CcNc[2]        = { 1, 2 };
  int    NuCode[kNNu]   = { kPdgNuMu, kPdgNuMuBar };
  int    NucCode[kNNuc] = { kPdgProton, kPdgNeutron };
  double W[kNW] = { 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0 };

  int  QrkCode[kNQrkMax];
  bool SeaQrk [kNQrkMax];

  gOutFile->cd();

  TTree * hmult = new TTree("hmult","hadronizer multiplicities");

  int   br_iev;    // event number
  int   br_nuc;    // hit nucleon PDG code
  int   br_neut;   // neutrino PDG code
  int   br_qrk;    // hit quark PDG code
  int   br_sea;    // hit quark is from the sea=1 (valence=0)
  int   br_ccnc;   // CC=1, NC=2
  float br_W;      // hadronic invariant mass
  int   br_np;     // number of generated p
  int   br_nn;     // number of generated n
  int   br_npip;   // number of generated pi+
  int   br_npim;   // number of generated pi-
  int   br_npi0;   // number of generated pi0
  int   br_nKp;    // number of generated K+
  int   br_nKm;    // number of generated K-
  int   br_nK0;    // number of generated K0

  hmult->Branch("iev",   &br_iev,   "iev/I"); 
  hmult->Branch("nuc",   &br_nuc,   "nuc/I");
  hmult->Branch("neut",  &br_neut,  "neut/I");
  hmult->Branch("qrk",   &br_qrk,   "qrk/I");
  hmult->Branch("sea",   &br_sea,   "sea/I");
  hmult->Branch("ccnc",  &br_ccnc,  "ccnc/I");
  hmult->Branch("W",     &br_W,     "W/F");
  hmult->Branch("np",    &br_np,    "np/I");
  hmult->Branch("nn",    &br_nn,    "nn/I");
  hmult->Branch("npip",  &br_npip,  "npip/I");
  hmult->Branch("npim",  &br_npim,  "npim/I");
  hmult->Branch("npi0",  &br_npi0,  "npi0/I");
  hmult->Branch("nKp",   &br_nKp,   "nKp/I");
  hmult->Branch("nKm",   &br_nKm,   "nKm/I");
  hmult->Branch("nK0",   &br_nK0,   "nK0/I");

  // CC/NC loop
  for(int iccnc=0; iccnc<2; iccnc++) { 
    InteractionType_t it = (CcNc[iccnc]==1) ? kIntWeakCC : kIntWeakNC;

    // neutrino & hit nucleon loops
    for(int inu=0; inu<kNNu; inu++) {
      for(int inuc=0; inuc<kNNuc; inuc++) {

        InitialState init (26,56,NuCode[inu]);
        ProcessInfo  proc (kScDeepInelastic, it);
        Interaction  intr (init, proc);

        intr.GetInitialStatePtr()->GetTargetPtr()->SetStruckNucleonPDGCode(NucCode[inuc]);

        // hit quark loop (if requested)
        int nqrk=1;
        if(gSetHitQrk) {
	  FillQrkArray(it, NuCode[inu], QrkCode, SeaQrk, kNQrkMax, nqrk);
        }
	for(int iqrk=0; iqrk<nqrk; iqrk++) {
           if(gSetHitQrk) {
             intr.GetInitialStatePtr()->GetTargetPtr()->SetStruckQuarkPDGCode(QrkCode[iqrk]);
             intr.GetInitialStatePtr()->GetTargetPtr()->SetStruckSeaQuark(SeaQrk[iqrk]);
           }

           LOG("main",pNOTICE) << "hadronizing: " << intr.AsString();

           for(int iw=0; iw<kNW; iw++) {
             intr.GetKinematicsPtr()->SetW(W[iw]);

             for(int in=0; in<nevents; in++) {
                TClonesArray * plist = model->Hadronize(&intr);
                assert(plist);

                br_iev  = in;
                br_nuc  = NucCode[inuc];
                br_neut = NuCode[inuc];
                br_qrk  = (gSetHitQrk) ? QrkCode[iqrk] : 0;
                br_sea  = (gSetHitQrk) ? SeaQrk[iqrk]  : 0;
                br_ccnc = CcNc[iccnc];
                br_W    = W[iw];
                br_np   = utils::fragmrec::NParticles(kPdgProton,  plist);
                br_nn   = utils::fragmrec::NParticles(kPdgNeutron, plist);
                br_npip = utils::fragmrec::NParticles(kPdgPiPlus,  plist);
                br_npim = utils::fragmrec::NParticles(kPdgPiMinus, plist);
                br_npi0 = utils::fragmrec::NParticles(kPdgPi0,     plist);
                br_nKp  = utils::fragmrec::NParticles(kPdgKPlus,   plist);
                br_nKm  = utils::fragmrec::NParticles(kPdgKMinus,  plist);
                br_nK0  = utils::fragmrec::NParticles(kPdgK0,      plist);

                hmult->Fill();

                plist->Delete();
                delete plist;
             }//n
          }//w
        }//qrk
      }//inuc
    }//inu
  }//cc/nc

  hmult->Write("hmult");
}
//____________________________________________________________________________
void testPhaseSpaceDecayer(int nevents, const HadronizationModelI * model)
{
  const int npmax = 100; // max number of particles
  const int kNW   = 2;   // n W values

  double W[kNW] = { 1.5, 4.0 };

  gOutFile->cd();

  TTree * hps = new TTree("hps","hadronizer events");

  int   br_iev;         // event number
  int   br_nuc;         // hit nucleon PDG code
  int   br_neut;        // neutrino PDG code
  int   br_qrk;         // hit quark PDG code
  int   br_ccnc;        // CC=1, NC=2
  float br_W;           // hadronic invariant mass
  int   br_np;          // number of particles in the fragmentation record
  int   br_pdg[npmax];  // PDG code of each particle
  int   br_ist[npmax];  // Status code of each particle
  float br_px [npmax];  // px of each particle
  float br_py [npmax];  // py of each particle
  float br_pz [npmax];  // pz of each particle
  float br_KE [npmax];  // kinematic energy of each particle

  hps->Branch("iev",   &br_iev,   "iev/I");
  hps->Branch("nuc",   &br_nuc,   "nuc/I");
  hps->Branch("neut",  &br_neut,  "neut/I");
  hps->Branch("qrk",   &br_qrk,   "qrk/I");
  hps->Branch("ccnc",  &br_ccnc,  "ccnc/I");
  hps->Branch("W",     &br_W,     "W/F");
  hps->Branch("np",    &br_np,    "np/I");
  hps->Branch("pdg",    br_pdg,   "pdg[np]/I");
  hps->Branch("ist",    br_ist,   "ist[np]/I");
  hps->Branch("px",     br_px,    "px[np]/F");
  hps->Branch("py",     br_py,    "py[np]/F");
  hps->Branch("pz",     br_pz,    "pz[np]/F");
  hps->Branch("KE",     br_KE,    "KE[np]/F");

  InitialState init (1056026000, kPdgNuMu);
  ProcessInfo  proc (kScDeepInelastic, kIntWeakCC);
  Interaction  intr (init, proc);

  intr.GetInitialStatePtr()->GetTargetPtr()->SetStruckNucleonPDGCode(kPdgProton);

  for(int iw=0; iw<kNW; iw++) {

    intr.GetKinematicsPtr()->SetW(W[iw]);

    for(int in=0; in<nevents; in++) {

       //rest prev branch
       for(int k=0; k<npmax; k++) {
	 br_pdg[k]=0; br_ist[k]=0;
         br_px[k]=0; br_py[k]=0; br_pz[k]=0; br_KE[k]=0;  
       }

       // hadronize
       TClonesArray * plist = model->Hadronize(&intr);
       assert(plist);

       br_iev   = in;
       br_nuc   = kPdgProton;
       br_neut  = kPdgNuMu;
       br_qrk   = 0;
       br_ccnc  = 1;
       br_W     = W[iw];
       br_np    = plist->GetEntries();

       TMCParticle * particle = 0;
       TIter particle_iter(plist);

       unsigned int i=0;

       while( (particle = (TMCParticle *) particle_iter.Next()) ) {
           br_pdg[i] = particle->GetKF();
           br_ist[i] = particle->GetKS();
           br_px[i]  = particle->GetPx();
           br_py[i]  = particle->GetPy();
           br_pz[i]  = particle->GetPz();
           br_KE[i]  = particle->GetEnergy() - particle->GetMass();

           hps->Fill();
           i++;
       } // particle-iterator

       plist->Delete();
       delete plist;
    }
  }
  hps->Write("hps");
}
//____________________________________________________________________________
void FillQrkArray(InteractionType_t it, int nu, 
                           int * QrkCode, bool * SeaQrk, int nmax, int & nqrk)
{
// utility method: create/fill array with all possible hit quarks
//

  for(int i=0; i<nmax; i++) {
    QrkCode[i] = -1;
  }
  if(it==kIntWeakNC) {
      nqrk=4;
      assert(nqrk<=nmax);
      QrkCode[0] = kPdgUQuark;  SeaQrk[0] = false;
      QrkCode[1] = kPdgUQuark;  SeaQrk[1] = true;
      QrkCode[2] = kPdgDQuark;  SeaQrk[2] = false;
      QrkCode[3] = kPdgDQuark;  SeaQrk[3] = true;
  }
  else if (it==kIntWeakCC) {
      nqrk=2;
      assert(nqrk<=nmax);
      if(pdg::IsNeutrino(nu)) {
         QrkCode[0] = kPdgDQuark;  SeaQrk[0] = false;
         QrkCode[1] = kPdgDQuark;  SeaQrk[1] = true;
      } else if(pdg::IsAntiNeutrino(nu)) {
         QrkCode[0] = kPdgUQuark;  SeaQrk[0] = false;
         QrkCode[1] = kPdgUQuark;  SeaQrk[1] = true;
      } else {
  	 exit(1);
      }
  } else {
      exit(1);
  }
}
//____________________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{
// Parse the command line arguments

  //number of events:
  try {
    LOG("Main", pINFO) << "Reading number of events to generate";
    gNEvents = genie::utils::clap::CmdLineArgAsInt(argc,argv,'n');
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("Main", pFATAL) << "Number of events was not specified";
      PrintSyntax();
      exit(1);
    }
  }

  //test id:
  try {
    LOG("Main", pINFO) << "Reading hadronizer test id";
    gTestId = genie::utils::clap::CmdLineArgAsInt(argc,argv,'t');
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("Main", pFATAL) << "No test id was specified";
      PrintSyntax();
      exit(1);
    }
  }

  // hadronizer:
  try {
    LOG("Main", pINFO) << "Reading hadronization algorithm name";
    gHadAlg = genie::utils::clap::CmdLineArgAsString(argc,argv,'a');
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("Main", pINFO) << "No hadronization algorithm was specified";
      PrintSyntax();
      exit(1);
    }
  }

  // hadronizer config:
  try {
    LOG("Main", pINFO) << "Reading hadronization algorithm config name";
    gHadConfig = genie::utils::clap::CmdLineArgAsString(argc,argv,'c');
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("Main", pINFO) << "No hadronization algorithm config was specified";
      PrintSyntax();
      exit(1);
    }
  }

  // set struck quark?
  try {
    LOG("Main", pINFO) << "reading struck quark option";
    gSetHitQrk = genie::utils::clap::CmdLineArgAsBool(argc,argv,'q');
  } catch(exceptions::CmdLineArgParserException e) {
      LOG("Main", pINFO) << "Using default option for setting hit quark";
  }
}
//____________________________________________________________________________
void PrintSyntax(void)
{
  LOG("Main", pNOTICE)
    << "\n\n" << "Syntax:" << "\n"
    << "  testHadronization -n nevents -t test -a hadronizer -c config [-q]\n";
}
//____________________________________________________________________________



