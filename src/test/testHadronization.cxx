//____________________________________________________________________________
/*!

\program testHadronization

\brief   test program used for testing/debugging the KNO & PYTHIA hadronizers

        Syntax :
           testHadronization -n nevents -t test -a hadronizer -c config [-q]

         Options :
           -n  number of events
           -a  hadronizer (algorithm name, eg genie::KNOHadronization)
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

void FillQrkArray(InteractionType_t it, int nu, 
             int * QrkCode, bool * SeaQrk, int nmax, int & nqrk);

TFile * gOutFile   = 0;
int     gNEvents   = -1;
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

  const int npmax = 100; // max number of particles

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

  TTree * hadnt = new TTree("hadnt","hadronizer multiplicities");

  int   br_iev;         // event number
  int   br_nuc;         // hit nucleon PDG code
  int   br_neut;        // neutrino PDG code
  int   br_qrk;         // hit quark PDG code
  int   br_sea;         // hit quark is from the sea=1 (valence=0)
  int   br_ccnc;        // CC=1, NC=2
  float br_W;           // hadronic invariant mass
  int   br_np;          // number of generated p
  int   br_nn;          // number of generated n
  int   br_npip;        // number of generated pi+
  int   br_npim;        // number of generated pi-
  int   br_npi0;        // number of generated pi0
  int   br_nKp;         // number of generated K+
  int   br_nKm;         // number of generated K-
  int   br_nK0;         // number of generated K0
  int   br_n;           // total number of generated particles
  int   br_pdg[npmax];  // PDG code of each particle
  int   br_ist[npmax];  // Status code of each particle
  float br_px [npmax];  // px of each particle
  float br_py [npmax];  // py of each particle
  float br_pz [npmax];  // pz of each particle
  float br_KE [npmax];  // kinetic energy of each particle
  float br_E  [npmax];  // energy of each particle
  float br_M  [npmax];  // mass of each particle
  float br_pL [npmax];  // pL of each particle
  float br_pT2[npmax];  // pT2 of each particle
  float br_xF [npmax];  // xF of each particle
  float br_z  [npmax];  // xF of each particle

  hadnt->Branch("iev",   &br_iev,   "iev/I"); 
  hadnt->Branch("nuc",   &br_nuc,   "nuc/I");
  hadnt->Branch("neut",  &br_neut,  "neut/I");
  hadnt->Branch("qrk",   &br_qrk,   "qrk/I");
  hadnt->Branch("sea",   &br_sea,   "sea/I");
  hadnt->Branch("ccnc",  &br_ccnc,  "ccnc/I");
  hadnt->Branch("W",     &br_W,     "W/F");
  hadnt->Branch("np",    &br_np,    "np/I");
  hadnt->Branch("nn",    &br_nn,    "nn/I");
  hadnt->Branch("npip",  &br_npip,  "npip/I");
  hadnt->Branch("npim",  &br_npim,  "npim/I");
  hadnt->Branch("npi0",  &br_npi0,  "npi0/I");
  hadnt->Branch("nKp",   &br_nKp,   "nKp/I");
  hadnt->Branch("nKm",   &br_nKm,   "nKm/I");
  hadnt->Branch("nK0",   &br_nK0,   "nK0/I");
  hadnt->Branch("n",     &br_n,     "n/I");
  hadnt->Branch("pdg",    br_pdg,   "pdg[n]/I");
  hadnt->Branch("ist",    br_ist,   "ist[n]/I");
  hadnt->Branch("px",     br_px,    "px[n]/F");
  hadnt->Branch("py",     br_py,    "py[n]/F");
  hadnt->Branch("pz",     br_pz,    "pz[n]/F");
  hadnt->Branch("KE",     br_KE,    "KE[n]/F");
  hadnt->Branch("E",      br_E,     "E[n]/F");
  hadnt->Branch("M",      br_M,     "M[n]/F");
  hadnt->Branch("pL",     br_pL,    "pL[n]/F");
  hadnt->Branch("pT2",    br_pT2,   "pT2[n]/F");
  hadnt->Branch("xF",     br_xF,    "xF[n]/F");
  hadnt->Branch("z",      br_z,     "z[n]/F");

  const int nnull_max=100;
  int nnull=0;

  // CC/NC loop
  for(int iccnc=0; iccnc<2; iccnc++) { 
    InteractionType_t it = (CcNc[iccnc]==1) ? kIntWeakCC : kIntWeakNC;

    // neutrino & hit nucleon loops
    for(int inu=0; inu<kNNu; inu++) {
      for(int inuc=0; inuc<kNNuc; inuc++) {

        InitialState init (26,56,NuCode[inu]);
        ProcessInfo  proc (kScDeepInelastic, it);
        Interaction  intr (init, proc);

        intr.InitStatePtr()->TgtPtr()->SetHitNucPdg(NucCode[inuc]);

        // hit quark loop (if requested)
        int nqrk=1;
        if(gSetHitQrk) {
	  FillQrkArray(it, NuCode[inu], QrkCode, SeaQrk, kNQrkMax, nqrk);
        }
	for(int iqrk=0; iqrk<nqrk; iqrk++) {
           if(gSetHitQrk) {
             intr.InitStatePtr()->TgtPtr()->SetHitQrkPdg(QrkCode[iqrk]);
             intr.InitStatePtr()->TgtPtr()->SetHitSeaQrk(SeaQrk[iqrk]);
           }

           LOG("main",pNOTICE) << "hadronizing: " << intr.AsString();

           for(int iw=0; iw<kNW; iw++) {
             intr.KinePtr()->SetW(W[iw]);

             for(int in=0; in<gNEvents; in++) {
                TClonesArray * plist = model->Hadronize(&intr);

                if(!plist) {
                 // don't count the current event and repeat
                 in--;
                 nnull++;
                 if(nnull>nnull_max) exit(1);
                 continue;
                }
                br_iev  = in;
                br_nuc  = NucCode[inuc];
                br_neut = NuCode[inu];
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
                br_n    = plist->GetEntries();

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
                   br_E[i]   = particle->GetEnergy();
                   br_M[i]   = particle->GetMass();
                   br_pL[i]  = particle->GetPz();
                   br_pT2[i] = TMath::Power(particle->GetPx(),2) + 
      		               TMath::Power(particle->GetPy(),2);
                   br_xF[i]  = particle->GetPz() / (W[iw]/2); 
		   br_z[i]   = particle->GetEnergy() / W[iw];

                   i++;
                } // particle-iterator

                hadnt->Fill();

                plist->Delete();
                delete plist;
             }//n
          }//w
        }//qrk
      }//inuc
    }//inu
  }//cc/nc

  hadnt->Write("hadnt");

  gOutFile->Close();

  return 0;
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
          << "  testHadronization -n nevents -a hadronizer -c config [-q]\n";
}
//____________________________________________________________________________



