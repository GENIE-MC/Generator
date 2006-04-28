//____________________________________________________________________________
/*!

\program testHadronization

\brief   test program used for testing/debugging the KNO & PYTHIA hadronizers

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

using std::string;
using std::vector;
using std::endl;
using std::ostringstream;

using namespace genie;

void testMultiplicities   (int n, const HadronizationModelI * m, string dn, string dt);
void testPhaseSpaceDecayer(int n, const HadronizationModelI * m, string dn, string dt);

TFile * gOutFile = 0;

//____________________________________________________________________________
int main(int argc, char ** argv)
{
  AlgFactory * algf = AlgFactory::Instance();

  const HadronizationModelI * KNO = 
     dynamic_cast<const HadronizationModelI *> (algf->GetAlgorithm(
                                     "genie::KNOHadronization", "Default"));

  gOutFile = new TFile("./genie-hadronization.root","recreate");

  //-- test KNO model multiplicities
  testMultiplicities(1000, KNO, "kno_mult", "KNO model multiplicities");

  //-- test KNO model phase space decayer
  testPhaseSpaceDecayer(100, KNO, "kno_decay"", "KNO phase space decayer");

  gOutFile->Close();

  return 0;
}
//____________________________________________________________________________
void testMultiplicities(
  int nevents, const HadronizationModelI * model, string dname, string dtitle)
{
  const int kNNu   = 2;
  const int kNNuc  = 2;
  const int kNW    = 8;

  int    CcNc[2]        = { 1, 2       };
  int    NuCode[kNNu]   = { 14, -14    };
  int    NucCode[kNNuc] = { 2212, 2112 };
  double W[kNW]         = { 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0 };
  //  double W[kNW]         = { 1.2, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0 };

  gOutFile->cd();

  TDirectory * dmult = new TDirectory(dname.c_str(), dtitle.c_str());
  dmult->cd();

  TTree * hmult = new TTree("hmult","hadronizer multiplicities");

  int   br_iev, br_nuc, br_neut, br_qrk, br_ccnc;
  int   br_np, br_nn, br_npip, br_npim, br_npi0, br_nKp, br_nKm, br_nK0;
  float br_W;

  hmult->Branch("iev",   &br_iev,   "iev/I");
  hmult->Branch("nuc",   &br_nuc,   "nuc/I");
  hmult->Branch("neut",  &br_neut,  "neut/I");
  hmult->Branch("qrk",   &br_qrk,   "qrk/I");
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

  for(int iccnc=0; iccnc<2; iccnc++) { 
    InteractionType_t it = (CcNc[iccnc]==1) ? kIntWeakCC : kIntWeakNC;

    for(int inu=0; inu<kNNu; inu++) {
      for(int inuc=0; inuc<kNNuc; inuc++) {

        InitialState init (26,56,NuCode[inu]);
        ProcessInfo  proc (kScDeepInelastic, it);
        Interaction  intr (init, proc);

        intr.GetInitialStatePtr()->GetTargetPtr()->SetStruckNucleonPDGCode(NucCode[inuc]);

        LOG("main",pNOTICE) << "hadronizing: " << intr.AsString();

        for(int iw=0; iw<kNW; iw++) {
           intr.GetKinematicsPtr()->SetW(W[iw]);

           for(int in=0; in<nevents; in++) {
              TClonesArray * plist = model->Hadronize(&intr);
              assert(plist);

              br_iev  = in;
              br_nuc  = NucCode[inuc];
              br_neut = NuCode[inuc];
              br_qrk  = 0;
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
      }//inuc
    }//inu
  }//cc/nc

  hmult->Write("hmult");
}
//____________________________________________________________________________
void testPhaseSpaceDecayer(
   int nevents, const HadronizationModelI * model, string dname, string dtitle)
{
  const int kNW = 2;  // n W values
   double W[kNW] = { 1.5, 4.0 };

  gOutFile->cd();

  TDirectory * dps = new TDirectory(dname.c_str(), dtitle.c_str());
  dps->cd();

  TTree * hps = new TTree("hps","hadronizer phase space decay");

  int   br_nuc, br_neut, br_qrk, br_ccnc, br_hpdgc;
  float br_W, br_px, br_py, br_pz, br_KE;

  hps->Branch("nuc",   &br_nuc,   "nuc/I");
  hps->Branch("neut",  &br_neut,  "neut/I");
  hps->Branch("qrk",   &br_qrk,   "qrk/I");
  hps->Branch("ccnc",  &br_ccnc,  "ccnc/I");
  hps->Branch("W",     &br_W,     "W/F");
  hps->Branch("hpdgc", &br_hpdgc, "hpdgc/I");
  hps->Branch("px",    &br_px,    "px/F");
  hps->Branch("py",    &br_py,    "py/F");
  hps->Branch("pz",    &br_pz,    "pz/F");
  hps->Branch("KE",    &br_KE,    "KE/F");

  InitialState init (1056026000,kPdgNuMu);
  ProcessInfo  proc (kScDeepInelastic, kIntWeakCC);
  Interaction  intr (init, proc);

  intr.GetInitialStatePtr()->GetTargetPtr()->SetStruckNucleonPDGCode(kPdgProton);

  for(int iw=0; iw<kNW; iw++) {

    intr.GetKinematicsPtr()->SetW(W[iw]);

    for(int in=0; in<nevents; in++) {

       TClonesArray * plist = model->Hadronize(&intr);
       assert(plist);

       TMCParticle * particle = 0;
       TIter particle_iter(plist);

       while( (particle = (TMCParticle *) particle_iter.Next()) ) {
           double E    = particle->GetEnergy();
           double Px   = particle->GetPx();
           double Py   = particle->GetPy();
           double Pz   = particle->GetPz();
           double m    = particle->GetMass();
           int    pdgc = particle->GetKF();

           br_nuc   = kPdgProton;
           br_neut  = kPdgNuMu;
           br_qrk   = 0;
           br_ccnc  = 1;
           br_hpdgc = pdgc;
           br_W     = W[iw];
           br_px    = Px;
           br_py    = Py;  
           br_pz    = Pz;
           br_KE    = E-m;

           hps->Fill();
       } // particle-iterator

       plist->Delete();
       delete plist;
    }
  }
  hps->Write("hps");
}
//____________________________________________________________________________

