// *********************************************************************
//  Validation of single kaon kinematics - Martti Nirkko (28th Nov 2014)
//  Compile and run in terminal:     root -l -b -q validation.C+g
// *********************************************************************

#include "code/singlekaon_xsec.cxx"
#include <TH1D.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TMath.h>
#include <TStyle.h>

void validation() {
  // Validate kinematical distributions of single kaon production
  // Have: E_nu, T_l, T_K, cos(theta_l), phi_Kq (generated randomly)
  // Want: Generate differential xsec histograms for these variables
  
  // ***************************
  // *   STEERING PARAMETERS   *
  // ***************************
  const int COMP = 10;                  // time complexity (runs with O(n^4)...)
  const int pLeptonIsUsed = 0;          // use momentum instead of kinetic energy
  const int thetaIsUsed = 0;            // use angle instead of cosine
  const double pi = TMath::Pi();
  
  // NOTE: To compare with the output of d4sigma_hist.C, use the same input variables above. Here, 
  //       for plotting the differential 1D cross-section, I want to use Tlep and costheta.
  
  const int type = 2;                   // lepton:   1=electron, 2=muon, 3=tau
  const int reac = 3;                   // reaction: 1=NN, 2=NP, 3=PP
  
  const double NEvents = 1e6;           // Number of events to generate
  
  // Initialise random number seed
  srand (time(NULL));
  
  double Enu;                           // neutrino energy [GeV]
  printf("Please enter neutrino energy: ");
  scanf("%lf", &Enu);
  printf("Trying to find input file for E_nu = %3.1lf GeV...\n", Enu);
  singlekaon_xsec *inst = new singlekaon_xsec();
  inst->init(Enu, type, reac);
  
  // Get output from d4sigma_plot.C to compare histograms
  std::string fname = Form("data/d4sigma_plot_%3.1lfGeV.root", Enu);
  TFile *infile = new TFile(fname.c_str());
  if (!infile) {
    printf("ERROR: File not found!");
    return;
  }
  TH1D *dTl = (TH1D*)infile->Get("dsigma_dTlepton");
  TH1D *dTk = (TH1D*)infile->Get("dsigma_dTkaon");
  TH1D *dct = (TH1D*)infile->Get("dsigma_dcostheta");
  TH1D *dph = (TH1D*)infile->Get("dsigma_dphikq");
  
  std::string strl;
  double ml = 0.;                     // lepton mass [GeV]
  if      (type==1) { ml = 0.510998928e-3;  strl = "e"; }
  else if (type==2) { ml = 0.1056583715;    strl = "#mu"; }
  else if (type==3) { ml = 1.77682;         strl = "#tau"; }
  else {std::cout<<"ERROR: Invalid lepton type!"<<std::endl; return;}
  
  std::string strK;
  double mK = 0.;                    // kaon mass [GeV]
  if      (reac==1) { mK = 0.493677; strK = "K^{+}"; }
  else if (reac==2) { mK = 0.497614; strK = "K^{0}"; }
  else if (reac==3) { mK = 0.493677; strK = "K^{+}"; }
  else {std::cout<<"ERROR: Invalid reaction!"<<std::endl; return;}
  
  std::string strN0, strN1;
  double mN = 0.;                    // nucleon mass [GeV]
  if      (reac==1) { mN = 0.939565379; strN0 = "n"; strN1 = "n"; }
  else if (reac==2) { mN = 0.939565379; strN0 = "n"; strN1 = "p"; }
  else if (reac==3) { mN = 0.938272046; strN0 = "p"; strN1 = "p"; }
  else {std::cout<<"ERROR: Invalid reaction!"<<std::endl; return;}
  
  std::string reaction = Form( "%.2lf GeV #nu_{%s} %s #rightarrow %s^{-} %s %s", Enu, strl.c_str(), 
                                strN0.c_str(), strl.c_str(), strN1.c_str(), strK.c_str() );
  
  // Number of bins
  const int NBINS = 10*COMP;
    
  // Initialize histograms
  TH1D *hTK = new TH1D("hTK", "Kaon kinetic energy",     NBINS,    0., Enu-mK-ml);
  TH1D *hTl = new TH1D("hTl", "Lepton kinetic energy",   NBINS,    0., Enu-mK-ml);
  TH1D *hct = new TH1D("hct", "Lepton polar angle",      NBINS,   -1., 1.);
  TH1D *hph = new TH1D("hph", "Kaon azimuthal angle",    NBINS/5, -pi, pi);
  
  TH1D *hQ2 = new TH1D("hQ2", "Momentum transfer",       NBINS,    0., Enu-mK-ml);
  TH1D *hTN = new TH1D("hTN", "Nucleon kinetic energy",  NBINS,    0., Enu-mK-ml);
  TH1D *hcl = new TH1D("hcl", "Lepton polar angle",      NBINS,   -1., 1.);
  TH1D *hcK = new TH1D("hcK", "Kaon polar angle",        NBINS,   -1., 1.);
  TH1D *hcN = new TH1D("hcN", "Nucleon polar angle",     NBINS,   -1., 1.);
  TH1D *hpl = new TH1D("hpl", "Lepton azimuthal angle",  NBINS,   -pi, pi);
  TH1D *hpK = new TH1D("hpK", "Kaon azimuthal angle",    NBINS,   -pi, pi);
  TH1D *hpN = new TH1D("hpN", "Nucleon azimuthal angle", NBINS,   -pi, pi);
  
  // Histograms for momenta (measured by detector)
  double maxmom = sqrt(pow(Enu-mK-ml + mN, 2) - mN*mN);
  TH1D *hmoml = new TH1D("hmoml", "Lepton momentum",  NBINS,  0., maxmom);
  TH1D *hmomK = new TH1D("hmomK", "Kaon momentum",    NBINS,  0., maxmom);
  TH1D *hmomN = new TH1D("hmomN", "Nucleon momentum", NBINS,  0., maxmom);
  TH1D *hmomlx = new TH1D("hmomlx", "Lepton momentum",  NBINS,  0., maxmom);
  TH1D *hmomKx = new TH1D("hmomKx", "Kaon momentum",    NBINS,  0., maxmom);
  TH1D *hmomNx = new TH1D("hmomNx", "Nucleon momentum", NBINS,  0., maxmom);
  TH1D *hmomly = new TH1D("hmomly", "Lepton momentum",  NBINS,  0., maxmom);
  TH1D *hmomKy = new TH1D("hmomKy", "Kaon momentum",    NBINS,  0., maxmom);
  TH1D *hmomNy = new TH1D("hmomNy", "Nucleon momentum", NBINS,  0., maxmom);
  TH1D *hmomlz = new TH1D("hmomlz", "Lepton momentum",  NBINS,  0., maxmom);
  TH1D *hmomKz = new TH1D("hmomKz", "Kaon momentum",    NBINS,  0., maxmom);
  TH1D *hmomNz = new TH1D("hmomNz", "Nucleon momentum", NBINS,  0., maxmom);
  
  // Initialise all needed parameters
  double TK_max, TK, Tl_max, Tl, pl_max, pl, costh_l, phi_l, phi_Kq, xsec;
  double El, theta, Q2;
  TLorentzVector P4_nu, P4_lep, q;
  TVector3 lep3, p_l, p_Q, p_K, p_N, p_Kq, p_Nq;
  double TN, EK, pK, pN; 
  double pQ, costh_Kq, theta_Kq, rot_Z, rot_X;
  double costh_K, costh_N, phi_K, phi_N;
  int sign;
  
  // Unit 3-vectors and neutrino momentum
  const TVector3 e1(1,0,0);
  const TVector3 e2(0,1,0);
  const TVector3 e3(0,0,1);
  const TVector3 p_nu = Enu*e3;
  
  // Generate large number of events
  for (int i=0; i<(int)NEvents; i++) {
    xsec = -999.;
    while(xsec <= 0) {                                        // retry until valid kinematics are found
      TK_max  = Enu - mK - ml;                                // maximal allowed kaon kinetic energy
      TK      = TK_max*(rand()/(double)RAND_MAX);             // kaon kinetic energy [0, TK_max]
      Tl_max  = Enu - mK - ml;                                // maximal allowed lepton kinetic energy (uncorrelated!)
      pl_max  = sqrt((Tl_max+ml)*(Tl_max+ml)-ml*ml);          // maximal allowed lepton momentum
      if (pLeptonIsUsed) {
        pl      = pl_max*(rand()/(double)RAND_MAX);           // lepton momentum [0, pl_max]
        Tl      = sqrt(pl*pl+ml*ml)-ml;
      } else {
        Tl      = Tl_max*(rand()/(double)RAND_MAX);           // lepton kinetic energy [0, Tl_max]
      }
      if (thetaIsUsed) {
        theta   = pi*(rand()/(double)RAND_MAX);               // lepton polar angle [0, pi]
        costh_l = cos(theta);
      } else {
        costh_l = 2.*(rand()/(double)RAND_MAX)-1.;            // cosine of lepton polar angle [-1, 1]
        theta = acos(costh_l);
      }
      phi_l   = pi*(2.*(rand()/(double)RAND_MAX)-1.);         // lepton azimuthal angle [-pi, pi]
      phi_Kq  = pi*(2.*(rand()/(double)RAND_MAX)-1.);         // kaon azimuthal angle [-pi, pi]
      xsec = inst->diffxsec(Tl,TK,theta,phi_Kq);              // 4D-differential cross-section
    }
    
    // INFO: In dsig_dQ2.f, the value for cos(theta) is changed to the following:
    //   cos(theta)=(2.0*Enu*El-aml*aml+aq2)/(2.0*Enu*alepvec);
    // aq2 is the (negative) integration variable for Q^2, going from 0 to -Enu
    
    // Calculate total energy and momentum of particles
    El = Tl + ml; pl = sqrt(El*El - ml*ml);
    EK = TK + mK; pK = sqrt(EK*EK - mK*mK);
    
    // Calculate Q^2 from lepton kinematics
    P4_nu = TLorentzVector(0.,0.,Enu,Enu);
    lep3 = TVector3(0.,0.,0.);
    lep3.SetMagThetaPhi(pl,theta,0.);
    P4_lep = TLorentzVector(lep3,El);
    q = P4_nu - P4_lep;
    Q2 = -q.Mag2();                                 // --> Q^2 histogram
    
    // Calculate 3-momentum of lepton
    p_l[0] = pl*sin(theta)*cos(phi_l);
    p_l[1] = pl*sin(theta)*sin(phi_l);
    p_l[2] = pl*costh_l;
    
    // Calculate 3-momentum of momentum transfer to hadronic system
    p_Q = -p_l;
    p_Q[2] += Enu;
    pQ = p_Q.Mag();
    TVector3 dirQ = p_Q.Unit();
    
    // In the momentum transfer plane, new angles... (see eq.17 of notes)
    costh_Kq = (pQ*pQ + pK*pK + mN*mN - (Enu-El-EK+mN)*(Enu-El-EK+mN))/(2.*pQ*pK);
    theta_Kq = acos(costh_Kq);
    
    // Get rotation angle from z-axis to Q^2-direction
    rot_Z = acos(e3*(p_Q.Unit()));         // rotate counter-clockwise, towards x-axis by [0,pi]
    
    // Get rotation angle from x-axis to Q^2-direction (projected to transverse plane)
    sign = 0;
    if (p_Q[1] > 0) sign = +1;
    else            sign = -1;
    TVector3 p_Qt(p_Q[0],p_Q[1],0);
    rot_X = sign*acos(e1*(p_Qt.Unit()));   // rotate around z-axis by [-pi,pi]
    
    // Calculate 3-momentum of kaon (in Q-frame)
    p_Kq[0] = pK*sin(theta_Kq)*cos(phi_Kq);
    p_Kq[1] = pK*sin(theta_Kq)*sin(phi_Kq);
    p_Kq[2] = pK*costh_Kq;
    
    // Calculate 3-momentum of nucleon (in Q-frame)
    p_Nq = -p_Kq;
    p_Nq[2] += pQ;
    
    // Rotate particles into correct frame
    p_K = p_Kq;
    p_K.RotateY(rot_Z);
    p_K.RotateZ(rot_X);
    p_N = p_Nq;
    p_N.RotateY(rot_Z);
    p_N.RotateZ(rot_X);
    
    // Kinetic energy of nucleon
    pN = p_N.Mag();
    TN = sqrt(pN*pN+mN*mN) - mN;                    // --> TN histogram
    
    // Momenta of particles
    hmoml->Fill(pl, xsec);
    hmomK->Fill(pK, xsec);
    hmomN->Fill(pN, xsec);
    hmomlx->Fill(p_l[0], xsec);
    hmomKx->Fill(p_K[0], xsec);
    hmomNx->Fill(p_N[0], xsec);
    hmomly->Fill(p_l[1], xsec);
    hmomKy->Fill(p_K[1], xsec);
    hmomNy->Fill(p_N[1], xsec);
    hmomlz->Fill(p_l[2], xsec);
    hmomKz->Fill(p_K[2], xsec);
    hmomNz->Fill(p_N[2], xsec);
    
    // Calculate polar angles w.r.t. neutrino
    costh_l = p_l.CosTheta();
    costh_K = p_K.CosTheta();
    costh_N = p_N.CosTheta();
    
    // CAUTION: TESTING!!!
    // Rotate to lepton frame
    p_l.RotateUz(p_l.Unit());
    p_K.RotateUz(p_l.Unit());
    p_N.RotateUz(p_l.Unit());
    
    // Calculate azimuthal angles w.r.t. neutrino (z-axis)
    phi_l = p_l.Phi() - p_l.Phi();
    phi_K = p_K.Phi() - p_l.Phi();
    phi_N = p_N.Phi() - p_l.Phi();
        
    // Multiply Jacobians, if needed
    //if (thetaIsUsed)    xsec *= sin(theta);     // costh -> theta
    if (!pLeptonIsUsed) xsec *= El/pl;      // pl -> Tlep
    
    // Weight entries by xsec
    hTK->Fill(TK,      xsec);
    hTl->Fill(Tl,      xsec);
    hct->Fill(costh_l, xsec);
    hph->Fill(phi_Kq,  xsec);
    
    hQ2->Fill(Q2, xsec);
    hTN->Fill(TN, xsec);
    hcl->Fill(costh_l, xsec);
    hcK->Fill(costh_K, xsec);
    hcN->Fill(costh_N, xsec);
    hpl->Fill(phi_l, xsec);   // should always be zero
    hpK->Fill(phi_K, xsec);
    hpN->Fill(phi_N, xsec);
    
    if (i%((int)NEvents/20)==0) cout << "Processing... (" << i/((int)NEvents/100) << "%)" << endl;
  }
  
  // Scale histograms to cross-section
  double totalxsec = dTk->Integral();
  hTK->Scale(totalxsec/hTK->Integral());
  hTl->Scale(totalxsec/hTl->Integral());
  hct->Scale(totalxsec/hct->Integral());
  hph->Scale(totalxsec/hph->Integral());
  hQ2->Scale(totalxsec/hQ2->Integral());
  hTN->Scale(totalxsec/hTN->Integral());
  hcl->Scale(totalxsec/hcl->Integral());
  hcK->Scale(totalxsec/hcK->Integral());
  hcN->Scale(totalxsec/hcN->Integral());
  hpl->Scale(totalxsec/hpl->Integral());
  hpK->Scale(totalxsec/hpK->Integral());
  hpN->Scale(totalxsec/hpN->Integral());
  
  // Make lines thicker & colourful
  hTK->SetLineWidth(2); hTK->SetLineColor(kRed);
  hTl->SetLineWidth(2); hTl->SetLineColor(kRed);
  hct->SetLineWidth(2); hct->SetLineColor(kRed);
  hph->SetLineWidth(2); hph->SetLineColor(kRed);
  hQ2->SetLineWidth(2); //hQ2->SetLineColor(kBlue);
  hTN->SetLineWidth(2); hTN->SetLineColor(kOrange);
  hcl->SetLineWidth(2); hcl->SetLineColor(kRed);
  hcK->SetLineWidth(2); hcK->SetLineColor(kBlack);
  hcN->SetLineWidth(2); hcN->SetLineColor(kOrange);
  hpl->SetLineWidth(2); hpl->SetLineColor(kRed);
  hpK->SetLineWidth(2); hpK->SetLineColor(kBlack);
  hpN->SetLineWidth(2); hpN->SetLineColor(kOrange);
  hmoml->SetLineWidth(2); hmoml->SetLineColor(kRed);
  hmomK->SetLineWidth(2); hmomK->SetLineColor(kBlack);
  hmomN->SetLineWidth(2); hmomN->SetLineColor(kOrange);
  
  hmomlx->SetLineWidth(2); hmomlx->SetLineColor(kRed);
  hmomKx->SetLineWidth(2); hmomKx->SetLineColor(kBlack);
  hmomNx->SetLineWidth(2); hmomNx->SetLineColor(kOrange);
  hmomly->SetLineWidth(2); hmomly->SetLineColor(kRed);
  hmomKy->SetLineWidth(2); hmomKy->SetLineColor(kBlack);
  hmomNy->SetLineWidth(2); hmomNy->SetLineColor(kOrange);
  hmomlz->SetLineWidth(2); hmomlz->SetLineColor(kRed);
  hmomKz->SetLineWidth(2); hmomKz->SetLineColor(kBlack);
  hmomNz->SetLineWidth(2); hmomNz->SetLineColor(kOrange);
  
  // Play with stat box
  gStyle->SetOptStat(0);
  //gStyle->SetStatX(0.5);
  //gStyle->SetStatY(0.88);
  
  // Set logarithmic scale on/off
  int logscale = 0;
  if(!logscale) {
    hTK->SetMinimum(0);
    hTl->SetMinimum(0);
    hct->SetMinimum(0);
    hph->SetMinimum(0);
    
    hQ2->SetMinimum(0);
    hTN->SetMinimum(0);
    hcl->SetMinimum(0); hcK->SetMinimum(0); hcN->SetMinimum(0);
    hpl->SetMinimum(0); hpK->SetMinimum(0); hpN->SetMinimum(0);
    
    hmoml->SetMinimum(0); hmomK->SetMinimum(0); hmomN->SetMinimum(0);
  }
  
  hTK->SetMaximum(1.2*max(hTK->GetMaximum(),dTk->GetMaximum()));
  hTl->SetMaximum(1.2*max(hTl->GetMaximum(),dTl->GetMaximum()));
  hct->SetMaximum(1.2*max(hct->GetMaximum(),dct->GetMaximum()));
  hph->SetMaximum(1.2*max(hph->GetMaximum(),dph->GetMaximum()));
  
  std::string outname = Form("%3.1lfGeV_", Enu);
  
  // Plot distributions of primary variables
  TCanvas *c1 = new TCanvas("c1", "Kinematics of single kaon production", 1200, 900);
  c1->SetTitle(reaction.c_str());
  c1->Divide(2,2);
  TPad* p[4];
  for (int k=0; k<4; k++) {
    p[k] = (TPad*)c1->GetPad(k+1);
    if (logscale) p[k]->SetLogy();
  }
  c1->cd(1);
  hTl->SetXTitle("T_{#mu} [GeV]");
  hTl->Draw(); dTl->Draw("same");
  c1->cd(2);
  hTK->SetXTitle("T_{K} [GeV]");
  hTK->Draw(); dTk->Draw("same");
  c1->cd(3);
  hct->SetXTitle("cos(#theta_{#mu}) [ ]");
  hct->Draw(); dct->Draw("same");
  c1->cd(4);
  hph->SetXTitle("#phi_{Kq} [rad]");
  hph->Draw(); dph->Draw("same");
  c1->Print(("images/kinematics_"+outname+"01_input.png").c_str()); c1->Close();
  
  // Plot distributions of secondary variables (e.g. Q^2, cos(theta_kaon), ...)
  TCanvas *c2 = new TCanvas("c2", "Q^2", 1200, 900);
  c2->SetTitle(reaction.c_str());
  c2->Divide(2,2);
  //TPad* p[4];
  for (int k=0; k<4; k++) {
    p[k] = (TPad*)c2->GetPad(k+1);
    if (logscale) p[k]->SetLogy();
  }
  c2->cd(1);
  hQ2->SetTitle("Momentum transfer to nucleon"); hQ2->SetXTitle("Q^{2} [GeV^{2}]");
  hQ2->SetMaximum(1.2*hQ2->GetMaximum());
  hQ2->Draw();
  c2->cd(2);
  hTl->SetTitle("Kinetic energy of particle"); hTl->SetXTitle("T_{x} [GeV]");
  hTl->SetLineColor(kRed); hTK->SetLineColor(kBlack);
  hTl->SetMaximum(max(1.2*hTN->GetMaximum(),max(hTl->GetMaximum(),hTK->GetMaximum())));
  hTl->Draw(); hTK->Draw("same"); hTN->Draw("same");
  c2->cd(3);
  hcl->SetTitle("Polar angle w.r.t. neutrino"); hcl->SetXTitle("cos(#theta_{#nux}) [ ]");
  hcl->SetMaximum(1.2*max(hcN->GetMaximum(),max(hcl->GetMaximum(),hcK->GetMaximum())));
  hcl->Draw(); hcK->Draw("same"); hcN->Draw("same");
  c2->cd(4);
  hpl->SetTitle("Azimuthal angle w.r.t. lepton"); hpl->SetXTitle("#phi_{lx} [rad]");
  hpl->SetMaximum(1.2*max(hpK->GetMaximum(),hpN->GetMaximum()));
  hpl->Draw(); hpK->Draw("same"); hpN->Draw("same");
  c2->Print(("images/kinematics_"+outname+"02_output.png").c_str()); c2->Close();
  
  // Plot momenta of final state particles (measured by the detectors)
  TCanvas *c3 = new TCanvas("c3", "Momenta in lab frame", 1200, 900);
  c3->SetTitle(reaction.c_str());
  c3->Divide(2,2);
  for (int k=0; k<4; k++) {
    p[k] = (TPad*)c3->GetPad(k+1);
    if (logscale) p[k]->SetLogy();
  }
  c3->cd(1);
  hmomlx->SetTitle("X-momentum"); hmomlx->SetXTitle("p_{x} [GeV/c]");
  hmomlx->SetMaximum(1.2*max(hmomKx->GetMaximum(),hmomlx->GetMaximum()));
  hmomlx->Draw(); hmomKx->Draw("same"); hmomNx->Draw("same");
  c3->cd(2);
  hmomly->SetTitle("Y-momentum"); hmomly->SetXTitle("p_{y} [GeV/c]");
  hmomly->SetMaximum(1.2*max(hmomKy->GetMaximum(),hmomly->GetMaximum()));
  hmomly->Draw(); hmomKy->Draw("same"); hmomNy->Draw("same");
  c3->cd(3);
  hmomlz->SetTitle("Z-momentum"); hmomlz->SetXTitle("p_{z} [GeV/c]");
  hmomlz->SetMaximum(1.2*max(hmomKz->GetMaximum(),hmomlz->GetMaximum()));
  hmomlz->Draw(); hmomKz->Draw("same"); hmomNz->Draw("same");
  c3->cd(4);
  hmoml->SetTitle("Total momentum"); hmoml->SetXTitle("|p| [GeV/c]");
  hmoml->SetMaximum(1.2*max(hmomK->GetMaximum(),hmoml->GetMaximum()));
  hmoml->Draw(); hmomK->Draw("same"); hmomN->Draw("same");
  c3->Print(("images/kinematics_"+outname+"03_momenta.png").c_str()); c3->Close();
  
  // Generate rootfile for GENIE comparison
  TFile* outfile = new TFile(("data/kinematics_"+outname+"validation.root").c_str(), "RECREATE");
  // Kinetic energies of final state particles
  hTl->Write("T_lepton");
  hTK->Write("T_kaon");
  hTN->Write("T_nucleon");
  // Azimuthal angles w.r.t. neutrino
  hct->Write("cos(theta_lepton)");
  hcK->Write("cos(theta_kaon)");
  hcN->Write("cos(theta_nucleon)");
  // Radial angles w.r.t. neutrino (z-axis) and lepton (x-axis)
  hpl->Write("phi_l");
  hpK->Write("phi_kaon");
  hpN->Write("phi_nucleon");
  // Other variables related to momentum transfer to nucleus
  hph->Write("phi_Kq");
  hQ2->Write("Q_squared");
  outfile->Close();
  
  if (hTl) delete hTl;
  if (hTK) delete hTK;
  if (hTN) delete hTN;
  if (hct) delete hct;
  if (hcK) delete hcK;
  if (hcN) delete hcN;
  if (hph) delete hph;
  if (hpl) delete hpl;
  if (hpK) delete hpK;
  if (hpN) delete hpN;
  if (hQ2) delete hQ2;
  if (hmoml) delete hmoml;
  if (hmomK) delete hmomK;
  if (hmomN) delete hmomN;
  if (hmomlx) delete hmomlx;
  if (hmomKx) delete hmomKx;
  if (hmomNx) delete hmomNx;
  if (hmomly) delete hmomly;
  if (hmomKy) delete hmomKy;
  if (hmomNy) delete hmomNy;
  if (hmomlz) delete hmomlz;
  if (hmomKz) delete hmomKz;
  if (hmomNz) delete hmomNz;
  if (inst) delete inst;
}

