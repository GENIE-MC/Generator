//____________________________________________________________________________
/*!

\program gtestCOHNCGamma

\brief   Crude v0 program used for testing/debugging the COH NCGamma model

\author  Jon Sensenig inspired by gtestDISSF test file

\created July 23, 2020

*/
//____________________________________________________________________________

#include <string>
#include <vector>

#include <TFile.h>
#include <TTree.h>

#include "Framework/Algorithm/Algorithm.h"
#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/Utils/CmdLnArgParser.h"
#include "Framework/Utils/KineUtils.h"
#include "Physics/Coherent/XSection/COHFormFactorI.h"
#include "Physics/Coherent/XSection/COHFormFactorMap.h"
#include "Physics/Coherent/XSection/DeVriesFormFactor.h"
#include "Physics/Coherent/XSection/FourierBesselFFCalculator.h"

using namespace genie;
using namespace std;


void gtestCOHNCGamma  ();
void FormFactorTest(const TFile & f);
void FFInterpolateTest(const TFile & f);
void DeltaFFTest(const TFile & f);
void RTest();
void XSecTest(const TFile & f);


//__________________________________________________________________________
void gtestCOHNCGamma ()
{
  TFile file("./coh_ncgamma.root","recreate");

  // Modules under test
  FormFactorTest(file);    // DeVries form factors
  FFInterpolateTest(file); // Interpolation and Fourier Bessel coeffs.
  DeltaFFTest(file);       // Delta transition form factors
  RTest();                 // Delta current matrix traces
  XSecTest(file);          // Top level cross section

  file.Close();
}

//__________________________________________________________________________
void FormFactorTest(const TFile & f)
{

  cout << "\n ============= Beginning DeVries FF Test ============== \n" << endl;

  // -- define the out TTree
  TTree * ncgamma = new TTree("FF_DeVries","COH NC Gamma");
  int    brTgt;   // neutrino PDG code
  double brpFF;   // proton FF
  double brnFF;   // neutron FF
  double brQ;     // Q momentum
  bool   brNuc;    // Nucleus check

  ncgamma->Branch("tgt",  &brTgt,  "tgt/I");
  ncgamma->Branch("pFF",  &brpFF,  "pFF/D");
  ncgamma->Branch("nFF",  &brnFF,  "nFF/D");
  ncgamma->Branch("qmom",  &brQ,  "qmom/D");
  ncgamma->Branch("has_nuc",  &brNuc,  "has_nuc/B");

  // -- request the COH NC Gamma model
  AlgFactory * algf = AlgFactory::Instance();
  AlgId id("COHFormFactorMap","Default");
  const Algorithm * alg = algf->GetAlgorithm(id);
  const COHFormFactorMap * ffmap_alg = dynamic_cast<const COHFormFactorMap *>(alg);

  const int kNTgt  = 5;
  const int kNQ  = 20;
  int target   [kNTgt] = {1000010030, 1000060120, 1000020030, 1000180400, 1000260560}; // = { 3H, 12C, 3He, 40Ar, 56Fe };
  double Qmom    [kNQ] = {0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1};
  double pFF_arr [kNQ];
  double nFF_arr [kNQ];

  for(int itgt=0; itgt<kNTgt; itgt++) {

    brTgt = target[itgt];
    brNuc = ffmap_alg->HasNucleus(brTgt);
    if( !brNuc ) {
      cout << "Target " << brTgt << " not found " << endl;
      continue;
    }

    for(int iQ=0; iQ<kNQ; iQ++) {
      brQ = Qmom[iQ];
      brpFF = ffmap_alg->ProtonFF( brQ, brTgt );
      brnFF = ffmap_alg->NeutronFF( brQ, brTgt );

      pFF_arr[iQ] = brpFF;
      nFF_arr[iQ] = brnFF;

      if ( brpFF != brnFF ) { cout << "Warning Proton and Neutron FF differ!" <<endl; }
      cout << " Target " << brTgt << " Q = " << brQ << " Proton FF = " << brpFF << "  n_ff = " << brnFF << endl;

      ncgamma->Fill();
    }
  }

  TGraph * pFF_Q = new TGraph(kNQ, Qmom, pFF_arr);
  TGraph * nFF_Q = new TGraph(kNQ, Qmom, nFF_arr);

  pFF_Q->SetTitle("Proton FF vs Q");
  nFF_Q->SetTitle("neutron FF vs Q");
  pFF_Q->GetXaxis()->SetTitle("Q (GeV)");
  nFF_Q->GetXaxis()->SetTitle("Q (GeV)");
  pFF_Q->GetYaxis()->SetTitle("FF");
  nFF_Q->GetYaxis()->SetTitle("FF");

  ncgamma->Write();
  pFF_Q->Write();
  nFF_Q->Write();

  cout << "End of DeVries FF test!" << endl;
}

//__________________________________________________________________________
void FFInterpolateTest(const TFile & f)
{

  cout << "\n ============= Beginning FF Interpolation Test ============== \n" << endl;

  // -- define the out TTree
  TTree * ff_int = new TTree("FF_Int","COH NC Gamma");
  int    brTgt;   // neutrino PDG code
  double brpFF;   // proton FF
  double brnFF;   // neutron FF
  double brQ;     // Q momentum
  bool   brNuc;    // Nucleus check

  ff_int->Branch("tgt",  &brTgt,  "tgt/I");
  ff_int->Branch("pFF",  &brpFF,  "pFF/D");
  ff_int->Branch("nFF",  &brnFF,  "nFF/D");
  ff_int->Branch("qmom",  &brQ,  "qmom/D");
  ff_int->Branch("has_nuc",  &brNuc,  "has_nuc/B");

  // -- request the COH NC Gamma model
  AlgFactory * algf = AlgFactory::Instance();
  AlgId id("COHFormFactorInterpolation","Default");
  const Algorithm * alg = algf->GetAlgorithm(id);
  const COHFormFactorInterpolation * ffint_alg = dynamic_cast<const COHFormFactorInterpolation *>(alg);

  const int kNTgt  = 5;
  const int kNQ  = 20;
  int target   [kNTgt] = {1000010030, 1000060120, 1000020030, 1000180400, 1000260560}; // = { 3H, 12C, 3He, 40Ar, 56Fe };
  double Qmom    [kNQ]   = {0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1};
  double pFF_arr [kNQ];
  double nFF_arr [kNQ];

  for(int itgt=0; itgt<kNTgt; itgt++) {

    brTgt = target[itgt];
    brNuc = ffint_alg->HasNucleus(brTgt);
    if( !brNuc ) {
      cout << "Target " << brTgt << " not found " << endl;
      continue;
    }

    for(int iQ=0; iQ<kNQ; iQ++) {
      brQ = Qmom[iQ];
      brpFF = ffint_alg->ProtonFF( brQ, brTgt );
      brnFF = ffint_alg->NeutronFF( brQ, brTgt );

      pFF_arr[iQ] = brpFF;
      nFF_arr[iQ] = brnFF;

      if ( brpFF != brnFF ) { cout << "Warning Proton and Neutron FF differ!" <<endl; }
      cout << "Q " << brQ << " Proton FF = " << brpFF << "  n_ff = " << brnFF << endl;

      ff_int->Fill();
    }
  }

  TGraph * pFF_Q = new TGraph(kNQ, Qmom, pFF_arr);
  TGraph * nFF_Q = new TGraph(kNQ, Qmom, nFF_arr);

  pFF_Q->SetTitle("Interpolated Proton FF vs Q");
  nFF_Q->SetTitle("neutron FF vs Q");
  pFF_Q->GetXaxis()->SetTitle("Q (GeV)");
  nFF_Q->GetXaxis()->SetTitle("Q (GeV)");
  pFF_Q->GetYaxis()->SetTitle("FF");
  nFF_Q->GetYaxis()->SetTitle("FF");

  ff_int->Write();
  pFF_Q->Write();
  nFF_Q->Write();

  cout << "End of FF Interpolation test!" << endl;
}

//__________________________________________________________________________
void DeltaFFTest(const TFile & f)
{

  cout << "\n ============= Beginning Delta Transition FF Test ============== \n" << endl;

  // -- define the out TTree
  TTree * dt_ff = new TTree("DeltaTransFF","COH NC Gamma");
  double brC3V;   //  FF
  double brC3VNC; //  FF
  double brC5ANC; //  FF
  double brQ;     // Q momentum

  dt_ff->Branch("C3V",   &brC3V,   "C3V/D");
  dt_ff->Branch("C3VNC", &brC3VNC, "C3VNC/D");
  dt_ff->Branch("C5ANC", &brC5ANC, "C5ANC/D");
  dt_ff->Branch("qmom",  &brQ,     "qmom/D");

  // -- request the COH NC Gamma model
  AlgFactory * algf = AlgFactory::Instance();
  AlgId id("DeltaTransitionFormFactor","Default");
  const Algorithm * alg = algf->GetAlgorithm(id);
  const DeltaTransitionFormFactor * dt_ff_alg = dynamic_cast<const DeltaTransitionFormFactor *>(alg);

  const int kNQ  = 20;
  double Qmom      [kNQ]   = {0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1};
  double C3V_arr   [kNQ];
  double C3VNC_arr [kNQ];
  double C5ANC_arr [kNQ];

  for(int iQ=0; iQ<kNQ; iQ++) {

    brQ     = Qmom[iQ];
    brC3V   = dt_ff_alg->C3V(brQ);
    brC3VNC = dt_ff_alg->C3VNC(brQ);
    brC5ANC = dt_ff_alg->C5ANC(brQ);
                                                                                                              
    C3V_arr[iQ]   = brC3V;
    C3VNC_arr[iQ] = brC3VNC;
    C5ANC_arr[iQ] = brC5ANC;
                                                                                                              
    cout << "Q = " << brQ << " C3V = " << brC3V << " C3VNC = " << brC3VNC << " C5ANC = " << brC5ANC <<  endl;
    dt_ff->Fill();
  }

  TGraph * C3V   = new TGraph(kNQ, Qmom, C3V_arr);
  TGraph * C3VNC = new TGraph(kNQ, Qmom, C3VNC_arr);
  TGraph * C5ANC = new TGraph(kNQ, Qmom, C5ANC_arr);

  C3V->SetTitle("C3V vs Q");
  C3VNC->SetTitle("C3VNC vs Q");
  C5ANC->SetTitle("C5ANC vs Q");
  C3V->GetXaxis()->SetTitle("Q (GeV)");
  C3VNC->GetXaxis()->SetTitle("Q (GeV)");
  C5ANC->GetXaxis()->SetTitle("Q (GeV)");
  C3V->GetYaxis()->SetTitle("C3V");
  C3VNC->GetYaxis()->SetTitle("C3VNC");
  C5ANC->GetYaxis()->SetTitle("C5ANC");

  dt_ff->Write();
  C3V->Write();
  C3VNC->Write();
  C5ANC->Write();

  cout << "End of Delta Transistion FF test!" << endl;
}

//__________________________________________________________________________
void RTest()
{

  cout << "\n ============= Beginning Delta Current Trace Test ============== \n" << endl;

  int target = 1000180400; // 40Ar
  int probe  = 14;         // nu mu
  int prod   = 22;         // gamma
  double E   = 1.1;

  // -- request the COH NC Gamma model
  AlgFactory * algf = AlgFactory::Instance();
  AlgId id("COHDeltaCurrent","Default"); // combine
  const Algorithm * algR = algf->GetAlgorithm(id);
  const COHDeltaCurrent * Dcurr_alg = dynamic_cast<const COHDeltaCurrent *>(algR);

  // Get the DeVries form factors
  AlgId idFF("COHFormFactorMap","Default");
  const Algorithm * algFF = algf->GetAlgorithm(idFF);
  const COHFormFactorMap * ffmap_alg = dynamic_cast<const COHFormFactorMap *>(algFF);

  // Set up the interaction and kinematics
  Interaction * i = genie::Interaction::COHNC(target, probe, prod, E);
  Kinematics * kine = i->KinePtr();

  // Naively set kinematic values
  kine->SetHadSystP4(0.1,0.1,0.3,0.7);
  kine->SetFSLeptonP4(0.1,0.1,0.5,0.9);

  const int kNQ     = 8;
  double Qmom [kNQ] = {0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9}; 

  cout << "Resonance Delta 1232 code " << Dcurr_alg->Resonance() << endl;

  for(int iQ=0; iQ<kNQ; iQ++) {
    kine->SetQ2(pow( Qmom[iQ], 2 ));
    utils::math::GTrace R = Dcurr_alg->R(i, ffmap_alg);

    for ( unsigned int i=0; i<R.size(); i++ ) {
      for ( unsigned int j=0; j<R.size(); j++ ) {
        cout << "Q " <<  Qmom[iQ] << " Trace  [" << i << j << "] Re " << R[i][j].real() << " Im " << R[i][j].imag() << endl;
      }
    }
  }

  cout << "End of Delta Current Trace test!" << endl;
}

//__________________________________________________________________________
void XSecTest(const TFile & f)
{

  cout << "\n ============= Beginning X-Section Test ============== \n" << endl;

    // -- define the out TTree
  TTree * ncgamma_xsec = new TTree("XSec","COH NC Gamma");
  double brXsec; // cross section
  double brInt;  // Integral
  double brQ;    // Q
  double brE;    // energy

  ncgamma_xsec->Branch("xsec",  &brXsec,  "xsec/D");
  ncgamma_xsec->Branch("integral",  &brInt,  "integral/D");
  ncgamma_xsec->Branch("q",  &brQ,  "q/D");
  ncgamma_xsec->Branch("e",  &brE,  "e/D");

  int target = 1000180400; // 40Ar
  int probe  = 14;         // nu mu
  int prod   = 22;         // gamma

  // -- request the COH NC Gamma model
  AlgFactory * algf = AlgFactory::Instance();
  AlgId id("AlvarezRusoSalaCOHGammaPXSec","Default");
  const Algorithm * algXsec = algf->GetAlgorithm(id);
  const AlvarezRusoSalaCOHGammaPXSec * xsec_alg = dynamic_cast<const AlvarezRusoSalaCOHGammaPXSec *>(algXsec);

  // Not used for now so let it initialize to whatever
  KinePhaseSpace_t t;

  const int kNQ     = 20;
  double Qmom [kNQ] = {0.01,0.015,0.02,0.025,0.03,0.035,0.04,0.045,0.05,0.055,0.06,0.65,0.07,0.075,0.08,0.085,0.09,0.1,0.15,0.2};
  double Evec [kNQ] = {0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0};
  double Xsec_arr [kNQ];
  double Int_arr  [kNQ];

  for(int iQ=0; iQ<kNQ; iQ++) {
    for(int iE=0; iE<kNQ; iE++) {
      
      brQ = Qmom[iQ];
      brE = Evec[iE] + 0.5;

      // Set up the interaction and kinematics
      Interaction * i = Interaction::COHNC(target, probe, prod, brE);
      Kinematics * kine = i->KinePtr();
                                                                    
      // Naively set kinematic values
      kine->SetHadSystP4(0.01,0.01,0.03,brE*0.5);
      kine->SetFSLeptonP4(0.01,0.01,0.05,brE*0.2);
      kine->Setx(0.1);
      kine->Sety(0.15);

      bool vp = xsec_alg->ValidProcess(i);
      bool vk = xsec_alg->ValidKinematics(i);
      if ( (!vp) || (!vk) ) { 
        cout << "Valid Process " << vp << " Valid Kinematics " << vk << endl; 
        continue;
      }

      kine->SetQ2(pow( brQ, 2 ));

      brXsec = xsec_alg->XSec(i, t);
      brInt = xsec_alg->Integral(i);
      cout << "Q " << brQ << " E " << brE << " XSec " << brXsec << " XSec Integral " << brInt << endl;
     
      Xsec_arr[iE] = brXsec;
      Int_arr[iE] = brInt;
      ncgamma_xsec->Fill();
    } 
  }
 
  TGraph * gxsec = new TGraph(kNQ, Evec, Xsec_arr);
  TGraph * gint  = new TGraph(kNQ, Evec, Int_arr);

  gxsec->SetTitle("XSec vs E");
  gint->SetTitle("Integral vs E");
  gxsec->GetXaxis()->SetTitle("E (GeV)");
  gint->GetXaxis()->SetTitle("E (GeV)");
  gxsec->GetYaxis()->SetTitle("XSec");
  gint->GetYaxis()->SetTitle("Integral");

  ncgamma_xsec->Write();  
  gxsec->Write();
  gint->Write();

  cout << "End of NC Gamma XSec test!" << endl;
}
