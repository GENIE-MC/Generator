//____________________________________________________________________________
/*
 Copyright (c) 2003-2016, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Tingjun Yang (Stanford Univ.)

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Mar 05, 2009 - TY
   Was first added in v2.5.1.

*/
//____________________________________________________________________________

#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TLorentzVector.h>
#include <TProfile.h>
#include <TTree.h>

#include "EVGCore/EventRecord.h"
#include "GHEP/GHepParticle.h"
#include "Ntuple/NtpMCEventRecord.h"
#include "Ntuple/NtpMCRecHeader.h"
#include "Ntuple/NtpMCTreeHeader.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodes.h"
#include "validation/Hadronization/HadPlots.h"

using namespace genie;
using namespace genie::mc_vs_data;

//____________________________________________________________________________
HadPlots::HadPlots()
{

}
//____________________________________________________________________________
HadPlots::HadPlots(string name)
{
  modelName = name;
}
//____________________________________________________________________________
HadPlots::~HadPlots()
{

}
//____________________________________________________________________________
void HadPlots::SetName(string name)
{
  modelName = name;
}
//____________________________________________________________________________
void HadPlots::LoadData(string mcfile)
{
  mcFiles.push_back(mcfile);
}
//____________________________________________________________________________
void HadPlots::BookHists()
{
  LOG("gvldtest", pDEBUG) << "Booking histograms";

  const int nw = 14;
  const double W2bins[nw+1] = {1,2,3,4,6,8,12,16,23,32,45,63,90,125,225};
  const double Wbins [nw+1] = {1,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7,8,9,10};

  for (int i = 0; i<kMaxFiles; i++){

    for (int j = 0; j<5; j++){
      kno[j][i] = new TH1D();
      kno[j][i]->SetBins(100,0,4);
    }
    for (int j = 0; j<5; j++){
      npi0_nneg[j][i] = new TProfile();
      npi0_nneg[j][i]->SetBins(9,-1.5,7.5);
    } 
    for (int j = 0; j<3; j++){
      npi0_nch[j][i] = new TProfile();
      npi0_nch[j][i]->SetBins(19,-1.5,17.5);
    }

    npi0_nm[i] = new TProfile();
    npi0_nm[i]->SetBins(7,-1.5,5.5);
    npi0_nm_lo[i] = new TProfile();
    npi0_nm_lo[i]->SetBins(7,-1.5,5.5);
    npi0_nm_hi[i] = new TProfile();
    npi0_nm_hi[i]->SetBins(7,-1.5,5.5);
    
    xf_pip[i] = new TH1D();
    xf_pim[i] = new TH1D();
    xf_pip[i]->SetBins(40,-1,1);
    xf_pim[i]->SetBins(40,-1,1);
    fxf_pip[i] = new TH1D();
    fxf_pim[i] = new TH1D();
    fxf_pip[i]->SetBins(40,-1,1);
    fxf_pim[i]->SetBins(40,-1,1);

    Fxf_pos1[i] = new TH1D();
    Fxf_pro1[i] = new TH1D();
    Fxf_pos_kno1[i] = new TH1D();
    Fxf_pos2[i] = new TH1D();
    Fxf_neg1[i] = new TH1D();
    Fxf_neg_kno1[i] = new TH1D();
    Fxf_neg2[i] = new TH1D();
    Fxf_pos1[i]->SetBins(40,-1,1);
    Fxf_pro1[i]->SetBins(40,-1,1);
    Fxf_pos_kno1[i]->SetBins(40,-1,1);
    Fxf_pos2[i]->SetBins(40,-1,1);
    Fxf_neg1[i]->SetBins(40,-1,1);
    Fxf_neg_kno1[i]->SetBins(40,-1,1);
    Fxf_neg2[i]->SetBins(40,-1,1);

    Fxf_pos1_hi[i] = new TH1D();
    Fxf_pro1_hi[i] = new TH1D();
    Fxf_pos2_hi[i] = new TH1D();
    Fxf_neg1_hi[i] = new TH1D();
    Fxf_neg2_hi[i] = new TH1D();
    Fxf_pos1_hi[i]->SetBins(40,-1,1);
    Fxf_pro1_hi[i]->SetBins(40,-1,1);
    Fxf_pos2_hi[i]->SetBins(40,-1,1);
    Fxf_neg1_hi[i]->SetBins(40,-1,1);
    Fxf_neg2_hi[i]->SetBins(40,-1,1);

    z_pos[i] = new TH1D();
    z_neg[i] = new TH1D();
    z_pos[i]->SetBins(20,0,1);
    z_neg[i]->SetBins(20,0,1);

    z_E1[i] = new TH1D();
    z_E2[i] = new TH1D();
    z_E3[i] = new TH1D();
    z_E1[i]->SetBins(20,0,1);
    z_E2[i]->SetBins(20,0,1);
    z_E3[i]->SetBins(20,0,1);
    
    z1_pos[i] = new TH1D();
    z1_pro[i] = new TH1D();
    z1_neg[i] = new TH1D();
    z1_pos[i]->SetBins(20,0,1);
    z1_pro[i]->SetBins(20,0,1);
    z1_neg[i]->SetBins(20,0,1);
    z2_pos[i] = new TH1D();
    z2_pro[i] = new TH1D();
    z2_neg[i] = new TH1D();
    z2_pos[i]->SetBins(20,0,1);
    z2_pro[i]->SetBins(20,0,1);
    z2_neg[i]->SetBins(20,0,1);
  
    pt2_xf1[i] = new TH1D();
    pt2_xf2[i] = new TH1D();
    pt2_xf3[i] = new TH1D();
    pt2_xf4[i] = new TH1D();
    pt2_xf1[i]->SetBins(24,0,1.2);
    pt2_xf2[i]->SetBins(24,0,1.2);
    pt2_xf3[i]->SetBins(24,0,1.2);
    pt2_xf4[i]->SetBins(24,0,1.2);

    pt2_pip[i] = new TH1D();
    pt2_pim[i] = new TH1D();
    pt2_pip[i]->SetBins(24,0,1.2);
    pt2_pim[i]->SetBins(24,0,1.2);

    pt2_W2_F[i] = new TProfile();
    pt2_W2_B[i] = new TProfile();
    pt2_W2_F[i]->SetBins(nw-1,W2bins);
    pt2_W2_B[i]->SetBins(nw-1,W2bins);
  
    pt_W_F[i] = new TProfile();
    pt_W_B[i] = new TProfile();
    pt_W[i] = new TProfile();
    pt_W_F[i]->SetBins(nw,Wbins);
    pt_W_B[i]->SetBins(nw,Wbins);
    pt_W[i]->SetBins(nw,Wbins);
    
    pt2_xf_loW[i] = new TProfile();
    pt2_xf_hiW[i] = new TProfile();
    pt2_xf_loW[i]->SetBins(20,-1,1);
    pt2_xf_hiW[i]->SetBins(20,-1,1);

    neta_W[i] = new TProfile();
    neta_W[i] -> SetBins(40,0,10);
    neta_W_F[i] = new TProfile();
    neta_W_F[i] -> SetBins(40,0,10);
  }
}
//____________________________________________________________________________
void HadPlots::Analyze()
{
  LOG("gvldtest", pNOTICE) 
    << "Analyzing input files for model: " << modelName;

  if (!mcFiles.size()){
    LOG("gvldtest", pERROR) << "No MC files available";
    exit(0);
  }

  BookHists();

  // Loop over mc files
  for (unsigned imc = 0; imc < mcFiles.size(); imc++){
    LOG("gvldtest", pDEBUG) << "*** Trying File Number " << imc ;
    TFile fin(mcFiles[imc].c_str(),"READ");
    TTree * er_tree = 0;
    er_tree = dynamic_cast <TTree *>( fin.Get("gtree")  );
    if (!er_tree) {
      LOG("gvldtest", pERROR) 
         << "Null input GHEP event tree found in file: " << mcFiles[imc];
      return;
    }
    NtpMCEventRecord * mcrec = 0;
    er_tree->SetBranchAddress("gmcrec", &mcrec);   
    if (!mcrec) {
      LOG("gvldtest", pERROR) << "Null MC record";
      return;
    }
    Long64_t nmax = er_tree->GetEntries();
    if (nmax<0) {
      LOG("gvldtest", pERROR) << "Number of events = 0";
      return;
    }

    LOG("gvldtest", pNOTICE) 
     << "*** Analyzing: " << nmax << " events found in file " << mcFiles[imc];

    int im = -1; //interaction type, 0:vp, 1:vn, 2:vbarp, 3:vbarn

    const int nw = 14;
    double gerrx[100] = {0.0};
    double W2lo[nw] = { 1, 2, 3, 4, 6,  8, 12, 16, 23, 32, 45, 63,  90, 125 };
    double W2hi[nw] = { 2, 3, 4, 6, 8, 12, 16, 23, 32, 45, 63, 90, 125, 225 };
    
    double aW2[nw] = {0.0};
    double aW[nw]  = {0.0};
    
    double nch[nw]      = {0.0};   // no. of charge particles
    double errnch[nw]   = {0.0};    
    double nchpi[nw]    = {0.0};   // no. of charge pions    
    double nv[nw]       = {0};     // no. of neutrino interactions
    double nv2[nw]      = {0};    
    double pnch[13][nw] = {{0.}};
    double pnchnv[nw]   = {0.};    
    double nneg[nw]     = {0.0};   // mean no. of negative particles
    double errnneg[nw]  = {0.0};    
    double Dneg[nw]     = {0.0};   // dispersion D_=sqrt(<n_**2>-<n_>**2)         
    double D[nw]        = {0.0};   // dispersion of charged hadrons        
    double D_nch[nw]    = {0.0};   // D/<n_ch>         
    double nchkno[5]    = {0.0};   // calculate <n>
    double nkno[5]      = {0.0};        
    //pi0 distributions
    double npizero[nw]  = {0.0};   // no. of pi0's
    double errnpi0[nw]  = {0.0}; 
    double Dnpi0[nw]    = {0.0};    
    //F/B multiplicity
    double nchf[nw]     = {0.0};   // no. of forward charged particles
    double nchb[nw]     = {0.0};   // no. of backward charged particles    
    double nposf[nw]    = {0.0};
    double nnegf[nw]    = {0.0};
    double nposb[nw]    = {0.0};
    double nnegb[nw]    = {0.0};
    
    double cut1_nv = 0;
    double cut2_nv = 0;
    double cut4_nv = 0;
    
    double nv_E[3] = {0};
    
    double cut5_nv[nw] = {0.};
    double cut5a_nv = 0;
    double cut5b_nv = 0;
    double cut5d_nv = 0;
    double cut5e_nv = 0;
    double cut5f_nv = 0;
    double cut5g_nv = 0;
        
    //reweight to nubar spectrum
    double nubarw_data[10] = {44.5,303.7,458,463.1,446.6,316.3,272.9,169.3,325.3,60.3};
    double nubarw_mc[10] = {0.};

    for(Long64_t iev = 0; iev < nmax; iev++) {
      
      LOG("gvldtest", pDEBUG) << "*** Analyzing event: " << iev;

      er_tree->GetEntry(iev);

      NtpMCRecHeader rec_header = mcrec->hdr;
      EventRecord &  event      = *(mcrec->event);
      
      LOG("gvldtest", pDEBUG) << event;
      
      // go further only if the event is physical
      bool is_unphysical = event.IsUnphysical();
      if(is_unphysical) {
	LOG("gvldtest", pINFO) << "Skipping unphysical event";
	mcrec->Clear();
	continue;
      }

      // input particles and primary lepton
      GHepParticle * neutrino = event.Probe();
      assert(neutrino);
      GHepParticle * target = event.Particle(1); 
      assert(target);
      GHepParticle * fsl = event.FinalStatePrimaryLepton();
      assert(fsl);
      //GHepParticle * fsh = event.FinalStateHadronicSystem();
      GHepParticle * fsh = event.Particle(3);
      assert(fsh);
      GHepParticle * hitnucl = event.HitNucleon();
      if (!hitnucl) continue; 

      // find the interaction type
      if (neutrino->Pdg() == kPdgNuMu)
      {
	if      (target->Pdg() == kPdgProton ) { im = 0; }
	else if (target->Pdg() == kPdgNeutron) { im = 1; }
      }
      else 
      if (neutrino->Pdg() == kPdgAntiNuMu)
      {
	if      (target->Pdg() == kPdgProton ) { im = 2; }
	else if (target->Pdg() == kPdgNeutron) { im = 3; }
      }
      if (im<0)
      {
	LOG("gvldtest", pERROR) 
           << "Unexpected interaction: neutrino = " << neutrino->Pdg()
           << " target = " << target->Pdg();
	return;
      }
      
      double M = hitnucl->Mass(); //mass of struck nucleon

      TLorentzVector k1 = *(neutrino->P4());
      TLorentzVector k2 = *(fsl->P4());
      TLorentzVector p1 = *(hitnucl->P4());
      TLorentzVector p2 = *(fsh->P4());
      TLorentzVector q  = k1-k2; 
      double v  = (q*p1)/M;         // v (E transfer in hit nucleon rest frame)
      double Q2 = -1 * q.M2();      // momemtum transfer
      double y  = v*M/(k1*p1);      // Inelasticity, y = q*P1/k1*P1
      double W2 = M*M + 2*M*v - Q2; // Hadronic Invariant mass ^ 2
      double W  = TMath::Sqrt(W2);
      //double Phs = sqrt(pow(p2.Px(),2)+pow(p2.Py(),2)+pow(p2.Pz(),2));
      v+=M;  //measured v

      LOG("gvldtest", pDEBUG) 
        << "Q2 = " << Q2 << ", W = " << W << ", y = " << y << ", v = " << v;

      int np        = 0;
      int nn        = 0;
      int npip      = 0;
      int npim      = 0;
      int npi0      = 0;
      int nkp       = 0;
      int nkm       = 0;
      int nk0       = 0;
      int npbar     = 0;
      int ncharged  = 0;
      int npositive = 0;
      int nnegative = 0;
      int neta      = 0;

      GHepParticle * p = 0;
      TIter event_iter(&event);
      
      while((p=dynamic_cast<GHepParticle *>(event_iter.Next()))){
	if (p->Pdg() == kPdgEta) neta++;
	if (p->Status() == kIStStableFinalState && p->Pdg()!=kPdgMuon && p->Pdg()!=kPdgAntiMuon) {
	  if (p->Pdg() == kPdgProton)     np++;
	  if (p->Pdg() == kPdgNeutron)    nn++;
	  if (p->Pdg() == kPdgPiP)        npip++;
	  if (p->Pdg() == kPdgPiM)        npim++;
	  if (p->Pdg() == kPdgPi0)        npi0++;
	  if (p->Pdg() == kPdgKP)         nkp++;
	  if (p->Pdg() == kPdgKM)         nkm++;
	  if (p->Pdg() == kPdgK0)         nk0++;
	  if (p->Pdg() == kPdgAntiProton) npbar++;
	  if (p->Charge()){
	    ncharged++;
	    if (p->Charge()>0){
	      npositive++;
	    }
	    if (p->Charge()<0){
	      nnegative++;
	    }
	  }
	}
      }
      
      LOG("gvldtest",pDEBUG) 
        << "np = " << np << ", nn = " << nn
        << ", npip = " << npip << ", npim = " << npim << ", npi0 = " << npi0;

      double weight = 1.;

      //cuts in Phys.Rev.D25,624 (1982)
      bool cut6 = y>=0.1&&y<0.8;
      int ip_nubar = -1;
      if (W>1.0&&W<=1.5) ip_nubar = 0;
      else if (W>1.5&&W<=2.0) ip_nubar = 1;
      else if (W>2.0&&W<=2.5) ip_nubar = 2;
      else if (W>2.5&&W<=3.0) ip_nubar = 3;
      else if (W>3.0&&W<=3.5) ip_nubar = 4;
      else if (W>3.5&&W<=4.0) ip_nubar = 5;
      else if (W>4.0&&W<=4.5) ip_nubar = 6;
      else if (W>4.5&&W<=5.0) ip_nubar = 7;
      else if (W>5.0&&W<=7.5) ip_nubar = 8;
      else if (W>7.5&&W<=10.0) ip_nubar = 9;
      if (ip_nubar!=-1&&cut6) nubarw_mc[ip_nubar]+= weight;

      int ipos = -1; //W2 bin
      for (int idx = 0; idx<nw; idx++){
	if (W2>=W2lo[idx]&&W2<W2hi[idx]){
	  ipos = idx;
	}
      }
      if (ipos==-1) continue; //out of the energy range

      aW2[ipos] += weight*W2;

      nch[ipos]     += weight*(ncharged);
      errnch[ipos]  += weight*pow((double)ncharged,2);
      nchpi[ipos]   += weight*(npip+npim);
      nneg[ipos]    += weight*(nnegative);
      errnneg[ipos] += weight*pow((double)nnegative,2);
      nv[ipos]      += weight;
      nv2[ipos]     += weight*weight;
      npizero[ipos] += weight*npi0;
      errnpi0[ipos] += weight*pow((double)npi0,2);
      //prepare for KNO
      int nchtot = ncharged;
      //if (nchtot%2!=0) cout<<"Warning: nch = "<<nchtot<<endl;
      if (nchtot<=12&&nchtot>=0){
	pnchnv[ipos]+=weight;
	pnch[nchtot][ipos]+=weight;
      }

      if (W<2){
	nchkno[0] += weight*(ncharged);
	nkno[0]   += weight;
      }
      else if (W<3){
	nchkno[1] += weight*(ncharged);
	nkno[1]   += weight;
      }
      else if (W<4){
	nchkno[2] += weight*(ncharged);
	nkno[2]   += weight;
      }
      else if (W<5){
	nchkno[3] += weight*(ncharged);
	nkno[3]   += weight;
      }
      else {
	nchkno[4] += weight*(ncharged);
	nkno[4]   += weight;
      }

      //<pi0> <n-> correlation
      if (W>3&&W<=4){
	npi0_nneg[0][imc]->Fill(nnegative,npi0,weight);
      }
      else if (W>4&&W<=5){
	npi0_nneg[1][imc]->Fill(nnegative,npi0,weight);
      }
      else if (W>5&&W<=7){
	npi0_nneg[2][imc]->Fill(nnegative,npi0,weight);
      }
      else if (W>7&&W<=10){
	npi0_nneg[3][imc]->Fill(nnegative,npi0,weight);
      }
      else if (W>10&&W<=14){
	npi0_nneg[4][imc]->Fill(nnegative,npi0,weight);
      }	

      if (W<4){
	npi0_nm_lo[imc]->Fill(nnegative,npi0,weight);
      }
      else if (W>4){
	npi0_nm_hi[imc]->Fill(nnegative,npi0,weight);
      }
      
      if (W<4){
	npi0_nch[0][imc]->Fill(ncharged,npi0,weight);
      }
      else if (W<8){
	npi0_nch[1][imc]->Fill(ncharged,npi0,weight);
      }
      else if (W<15){
	npi0_nch[2][imc]->Fill(ncharged,npi0,weight);
      }
      
    }//end of first loop
    
    for (int i = 0; i<nw; i++){
      if (nv[i]){
	nch[i]   /= nv[i];
	D[i]      = sqrt(errnch[i]/nv[i]-pow(nch[i],2));
	errnch[i] = D[i]/sqrt(nv2[i]);
	nchpi[i] /= nv[i]*2; //1/2*(<pi+>+<pi->)
	nneg[i]  /= nv[i];
	Dneg[i]   = sqrt(errnneg[i]/nv[i]-pow(nneg[i],2));
	errnneg[i]= Dneg[i]/sqrt(nv2[i]);
	npizero[i]  /= nv[i];
	Dnpi0[i]  = sqrt(errnpi0[i]/nv[i]-pow(npizero[i],2));
	errnpi0[i]= Dnpi0[i]/sqrt(nv2[i]);
	D_nch[i]  = D[i]/nch[i];
	aW2[i]   /=nv[i];
	for (int j = 0; j<13; j++){
	  pnch[j][i] /= pnchnv[i];
	}
      }
    }

    //<nch> for KNO
    for (int i = 0; i<5; i++){
      nchkno[i]/=nkno[i];
    }

    //second loop
    for(Long64_t iev = 0; iev < nmax; iev++){
      
      er_tree->GetEntry(iev);
      
      NtpMCRecHeader rec_header = mcrec->hdr;
      EventRecord &  event      = *(mcrec->event);

      // go further only if the event is physical
      bool is_unphysical = event.IsUnphysical();
      if(is_unphysical) {
	LOG("gvldtest", pINFO) << "Skipping unphysical event";
	mcrec->Clear();
	continue;
      }

      //input particles and primary lepton
      GHepParticle * neutrino = event.Probe();
      assert(neutrino);
      GHepParticle * target = event.Particle(1); 
      assert(target);
      GHepParticle * fsl = event.FinalStatePrimaryLepton();
      assert(fsl);
      //GHepParticle * fsh = event.FinalStateHadronicSystem();
      GHepParticle * fsh = event.Particle(3);
      assert(fsh);
      GHepParticle * hitnucl = event.HitNucleon();
      if (!hitnucl) continue; //
      const Interaction * interaction = event.Summary();
      const ProcessInfo &  proc_info  = interaction->ProcInfo();

      double M = hitnucl->Mass(); //mass of struck nucleon
      double Ev = neutrino->E();  //neutrino energy
      double El = fsl->E();       //lepton energy

      TLorentzVector k1 = *(neutrino->P4());
      TLorentzVector k2 = *(fsl->P4());
      TLorentzVector p1 = *(hitnucl->P4());
      TLorentzVector p2 = *(fsh->P4());
      TLorentzVector q  = k1-k2; 
      double v  = (q*p1)/M;         // v (E transfer in hit nucleon rest frame)
      double Q2 = -1 * q.M2();      // momemtum transfer
      double x  = 0.5*Q2/(M*v);     // Bjorken x
      double y  = v*M/(k1*p1);      // Inelasticity, y = q*P1/k1*P1
      double W2 = M*M + 2*M*v - Q2; // Hadronic Invariant mass ^ 2
      double W  = TMath::Sqrt(W2);
      double Phs = sqrt(pow(p2.Px(),2)+pow(p2.Py(),2)+pow(p2.Pz(),2));
      double beta = Phs/p2.E();
      double gamma = p2.E()/sqrt(pow(p2.E(),2)-pow(Phs,2));
      v+=M;  //measured v

      int np = 0;
      int nn = 0;
      int npip = 0;
      int npim = 0;
      int npi0 = 0;
      int nkp = 0;
      int nkm = 0;
      int nk0 = 0;
      int npbar = 0;
      int ncharged = 0;
      int npositive = 0;
      int nnegative = 0;
      int neta = 0;

      vector<int> pid;
      vector<double> px;
      vector<double> py;
      vector<double> pz;
      vector<double> eng;
      vector<double> mass;
      vector<double> charge;
      GHepParticle * p = 0;
      TIter event_iter(&event);
    
      while((p=dynamic_cast<GHepParticle *>(event_iter.Next()))){

        if (p->Pdg() == kPdgEta)        neta++;
	if (p->Status() == kIStStableFinalState && p->Pdg()!=kPdgMuon && p->Pdg()!=kPdgAntiMuon) {
	  if (p->Pdg() == kPdgProton)     np++;
	  if (p->Pdg() == kPdgNeutron)    nn++;
	  if (p->Pdg() == kPdgPiP)        npip++;
	  if (p->Pdg() == kPdgPiM)        npim++;
	  if (p->Pdg() == kPdgPi0)        npi0++;
	  if (p->Pdg() == kPdgKP)         nkp++;
	  if (p->Pdg() == kPdgKM)         nkm++;
	  if (p->Pdg() == kPdgK0)         nk0++;
	  if (p->Pdg() == kPdgAntiProton) npbar++;

	  if (p->Charge()){
	    ncharged++;
	    if (p->Charge()>0){
	      npositive++;
	    }
	    if (p->Charge()<0){
	      nnegative++;
	    }
	  }
	  pid.push_back(p->Pdg());
	  px.push_back(p->Px());
	  py.push_back(p->Py());
	  pz.push_back(p->Pz());
	  eng.push_back(p->E());
	  mass.push_back(p->Mass());
	  charge.push_back(p->Charge());
	}
      }
      
      double weight = 1.;
      
      if (W<2){
	kno[0][imc]->Fill((ncharged)/nchkno[0],weight);
      }
      else if (W<3){
	kno[1][imc]->Fill((ncharged)/nchkno[1],weight);
      }
      else if (W<4){
	kno[2][imc]->Fill((ncharged)/nchkno[2],weight);
      }
      else if (W<5){
	kno[3][imc]->Fill((ncharged)/nchkno[3],weight);
      }
      else {
	kno[4][imc]->Fill((ncharged)/nchkno[4],weight);
      }

      int ipos = -1;
      for (int idx = 0; idx<nw; idx++){
	if (W2>=W2lo[idx]&&W2<W2hi[idx]){
	  ipos = idx;
	}
      }
      if (ipos==-1) continue; //out of the energy range
      
      //cuts in Phys.Rev.D27,1 (1983)
      bool cut0 = true;
      
      //cuts in Nucl.Phys.B214,369 (1983)
      bool cut1 = W>3&&Ev>5&&El>3&&x>0.1;
      
      //cuts in Z.Phys.C24,119 (1984)
      //note: experiment did correct for proton mis-identification, excellent!
      bool cut2 = W2>5&&Q2>1;
      
      //cuts in Z.Phys.C27,239 (1985)
      bool cut3  = Q2>1&&x<0.95&&El>4;
      //bool cut3a = W2<50&&Q2>1&&x<0.95&&El>4;
      //bool cut3b = W2>50&&Q2>1&&x<0.95&&El>4;
      bool cut3c = W2>9&&W2<25&&Q2>1&&x<0.95&&El>4;
      bool cut3d = W2>40&&Q2>1&&x<0.95&&El>4;
      
      //cuts in Phys.Rev.D19,1 (1979)
      bool cut4 = false;
      if (im==0){
	//cut4 = fullcut1(k1,k2,pid.size(),&pid[0],&px[0],&py[0],&pz[0])&&im==0&&Ev>10&&x>0.05&&y<0.9&&W2>16;
	cut4 = im==0&&Ev>10&&x>0.05&&y<0.9&&W2>16;
      }
      else {
	cut4 = false;
      }
      
      //cuts in Phys.Rev.D24,1071 (1981)
      //note: experiment did not correct for proton mis-identification
      bool cut5  = y>0.1&&y<0.8;
      //cut5 = true;
      bool cut5a = y>0.1&&y<0.8&&W<4&&Q2>=1&&Q2<45;
      bool cut5b = y>0.1&&y<0.8&&W>=4&&Q2>=1&&Q2<45;
      bool cut5c = y>0.1&&y<0.8&&Q2>=1&&Q2<45;
      bool cut5d = y>0.1&&y<0.8&&W>2&&W<4&&Q2<45;
      //bool cut5d = true;
      bool cut5e = y>0.1&&y<0.8&&W>2&&W<4&&Q2>1&&Q2<45;
      bool cut5f = y>0.1&&y<0.8&&W>4&&W<10&&Q2<45;
      bool cut5g = y>0.1&&y<0.8&&W>4&&W<10&&Q2>1&&Q2<45;
      
      
      //cuts in Phys.Rev.D25,624 (1982)
      bool cut6 = y>=0.1&&y<0.8;
      //note: experiment did not correct for proton mis-identification
      double wei_nubarp = 1.;
      int ip_nubar = -1;
      if (W>1.0&&W<=1.5) ip_nubar = 0;
      else if (W>1.5&&W<=2.0) ip_nubar = 1;
      else if (W>2.0&&W<=2.5) ip_nubar = 2;
      else if (W>2.5&&W<=3.0) ip_nubar = 3;
      else if (W>3.0&&W<=3.5) ip_nubar = 4;
      else if (W>3.5&&W<=4.0) ip_nubar = 5;
      else if (W>4.0&&W<=4.5) ip_nubar = 6;
      else if (W>4.5&&W<=5.0) ip_nubar = 7;
      else if (W>5.0&&W<=7.5) ip_nubar = 8;
      else if (W>7.5&&W<=10.0) ip_nubar = 9;
      if (ip_nubar!=-1&&cut6){
	if (nubarw_mc[ip_nubar]>0) {
	  if (neutrino->Pdg()==kPdgAntiNuMu) wei_nubarp = weight*nubarw_data[ip_nubar]/nubarw_mc[ip_nubar];
	  npi0_nm[imc]->Fill(nnegative,npi0,wei_nubarp);
	}
      }
      
      //count no. of (weighted) events
      if (cut1) {
	cut1_nv += weight;
      }
      if (cut2) {
	cut2_nv += weight;
      }
      
      if (proc_info.IsDeepInelastic()){
	if (Ev<15){
	  nv_E[0]+=weight;
	}
	else if (Ev<30){
	  nv_E[1]+=weight;
	}
	else {
	  nv_E[2]+=weight;
	}
      }
      
      if (cut4){
	cut4_nv+=weight;
      }
      
      if (cut5){
	cut5_nv[ipos] += weight;
	aW[ipos] += weight*sqrt(W2);
      }

      if (cut5a){
	cut5a_nv+=weight;
      }
      if (cut5b){
	cut5b_nv+=weight;
      }
      if (cut5d){
	cut5d_nv+=weight;
      }
      if (cut5e){
	cut5e_nv+=weight;
      }
      if (cut5f){
	cut5f_nv+=weight;
      }
      if (cut5g){
	cut5g_nv+=weight;
      }

      double meanpt = 0.;
      double meanpt_f = 0.;
      double meanpt_b = 0.;
      double meanpt2_f = 0.;
      double meanpt2_b = 0.;
      double nhch = 0.;
      double nhch_f = 0.;
      double nhch_b = 0.;
      double nhch_f2 = 0.;
      double nhch_b2 = 0.;
      
      for (unsigned ipar = 0; ipar<pid.size(); ipar++){//loop through hadrons
	int hadcharge = 0;
	if (charge[ipar]>0) hadcharge = 1;
	if (charge[ipar]<0) hadcharge = -1;
	double ptot = sqrt(pow(px[ipar],2)+pow(py[ipar],2)+pow(pz[ipar],2));
	double Pz = (px[ipar]*p2.Px()+py[ipar]*p2.Py()+pz[ipar]*p2.Pz())/Phs;
	double pt = sqrt(pow(ptot,2)-pow(Pz,2)); 
	double Ecm = gamma*(eng[ipar]-beta*Pz);
	Pz = gamma*(Pz-beta*eng[ipar]); //Lorentz boost
	double xf = 2*Pz/W;
	double z = eng[ipar]/v;
	double xf_uncor = xf;  //assign pion mass to proton
	double z_uncor = z;
	if (abs(pid[ipar])==2212&&ptot>1){ //proton mis-identified as pion for momenta above 1GeV
	  double Pz_uncor = (px[ipar]*p2.Px()+py[ipar]*p2.Py()+pz[ipar]*p2.Pz())/Phs;
	  double eng_uncor = sqrt(pow(ptot,2)+pow(0.1396,2));
	  Pz_uncor = gamma*(Pz_uncor-beta*eng_uncor); //Lorentz boost
	  xf_uncor = 2*Pz_uncor/W;
	  z_uncor  = eng_uncor/v;
	}
	  
	if (hadcharge){
	  meanpt += pt;
	  nhch++;
	  if (xf>0){
	    meanpt_f += pt;
	    nhch_f++;
	  }
	  if (xf<0){
	    meanpt_b += pt;
	    nhch_b++;
	  }
	}
	
	if (cut0){
	  if (xf>0){//forward
	    if (hadcharge){//charged hadrons
	      nchf[ipos]+=weight;
	    }
	  }	  
	  else if (xf<0){//backward
	    if (hadcharge){//charged hadrons
	      nchb[ipos]+=weight;
	    }
	  }
	}
	  
	if (cut1){
	  if (pid[ipar]==211){//pi+
	    xf_pip[imc]->Fill(xf,weight);
	    fxf_pip[imc]->Fill(xf,weight*Ecm*2/3.1416/W);
	  }
	  else if (pid[ipar]==-211){//pi-
	    xf_pim[imc]->Fill(xf,weight);
	    fxf_pim[imc]->Fill(xf,weight*Ecm*2/3.1416/W);
	  }
	}
	  
	if (cut2){
	  if (hadcharge>0){//postive
	    z_pos[imc]->Fill(z,weight);
	  }
	  else if (hadcharge<0){//negtive
	    z_neg[imc]->Fill(z,weight);
	  }
	}
	  
	if (hadcharge&&proc_info.IsDeepInelastic()){
	  if (Ev<15){
	    z_E1[imc]->Fill(z_uncor,weight);
	  }
	  else if (Ev<30){
	    z_E2[imc]->Fill(z_uncor,weight);
	  }
	  else {
	    z_E3[imc]->Fill(z_uncor,weight);
	  }
//	  if (pid[ipar]==211||pid[ipar]==2212&&pow(eng[ipar],2)-pow(mass[ipar],2)>1){
//	    pt2_pip[imc]->Fill(pt*pt,weight);
//	  }
//	  if (pid[ipar]==-211){
//	    pt2_pim[imc]->Fill(pt*pt,weight);
//	  }
	}
	  
	
	if (cut3&&hadcharge){
	  if (xf>0.3) {
	    meanpt2_f += pt*pt;
	    nhch_f2 ++;
	  }
	  if (xf<-0.3){
	    meanpt2_b += pt*pt;
	    nhch_b2 ++;
	  }
	}

	if (cut3c&&hadcharge){
	  pt2_xf_loW[imc]->Fill(xf,pt*pt,weight);
	}
	if (cut3d&&hadcharge){
	  pt2_xf_hiW[imc]->Fill(xf,pt*pt,weight);
	}
	  
	if (hadcharge&&cut4){
	  if (xf_uncor>0.0&&xf_uncor<.1) pt2_xf1[imc]->Fill(pt*pt,weight);
	  if (xf_uncor>0.1&&xf_uncor<.3) pt2_xf2[imc]->Fill(pt*pt,weight);
	  if (xf_uncor>0.3&&xf_uncor<.6) pt2_xf3[imc]->Fill(pt*pt,weight);
	  if (xf_uncor>0.6&&xf_uncor<1.) pt2_xf4[imc]->Fill(pt*pt,weight);
	}

	if (cut5){
	  if (xf_uncor>0){
	    if (hadcharge>0){
	      nposf[ipos]+=weight;
	    }
	    else if (hadcharge<0){
	      nnegf[ipos]+=weight;
	    }
	  }
	  if (xf_uncor<0){
	    if (hadcharge>0){
	      nposb[ipos]+=weight;
	    }
	    else if (hadcharge<0){
	      nnegb[ipos]+=weight;
	    }
	  }
	}
	if (cut5a&&xf_uncor>0){
	  if (hadcharge>0) {
	    z1_pos[imc]->Fill(z_uncor,weight);
	    if (abs(pid[ipar])==2212) z1_pro[imc]->Fill(z_uncor,weight);
	  }
	  else if (hadcharge<0){
	    z1_neg[imc]->Fill(z_uncor,weight);
	  }
	}
	if (cut5b&&xf_uncor>=0){
	  if (hadcharge>0) {
	    z2_pos[imc]->Fill(z_uncor,weight);
	    if (abs(pid[ipar])==2212) z2_pro[imc]->Fill(z_uncor,weight);
	  }
	  else if (hadcharge<0){
	    z2_neg[imc]->Fill(z_uncor,weight);
	  }
	}
	if (cut5d){
	  if (hadcharge>0){
	    Fxf_pos1[imc]->Fill(xf_uncor,weight*Ecm*2/3.1416/W);
	    if (abs(pid[ipar])==2212) Fxf_pro1[imc]->Fill(xf_uncor,weight*Ecm*2/3.1416/W);
	    //if (ihadmod==2) Fxf_pos_kno1[imc]->Fill(xf_uncor,weight*Ecm*2/3.1416/W);
	    //Fxf_pos1[imc]->Fill(xf,weight);
	  }
	  if (hadcharge<0){
	    //if (pid[ipar]==-211){
	    Fxf_neg1[imc]->Fill(xf_uncor,weight*Ecm*2/3.1416/W);
	    //if (ihadmod==2) Fxf_neg_kno1[imc]->Fill(xf_uncor,weight*Ecm*2/3.1416/W);
	    
	  }
	}
	if (cut5e){
	  if (hadcharge>0){
	    Fxf_pos2[imc]->Fill(xf_uncor,weight*Ecm*2/3.1416/W);
	  }
	  if (hadcharge<0){
	    Fxf_neg2[imc]->Fill(xf_uncor,weight*Ecm*2/3.1416/W);
	  }
	}
	if (cut5f){
	  if (hadcharge>0){
	    Fxf_pos1_hi[imc]->Fill(xf_uncor,weight*Ecm*2/3.1416/W);
	    if (abs(pid[ipar])==2212) Fxf_pro1_hi[imc]->Fill(xf_uncor,weight*Ecm*2/3.1416/W);
	  }
	  if (hadcharge<0){
	    Fxf_neg1_hi[imc]->Fill(xf_uncor,weight*Ecm*2/3.1416/W);
	  }
	}
	if (cut5g){
	  if (hadcharge>0){
	    Fxf_pos2_hi[imc]->Fill(xf_uncor,weight*Ecm*2/3.1416/W);
	  }
	  if (hadcharge<0){
	    Fxf_neg2_hi[imc]->Fill(xf_uncor,weight*Ecm*2/3.1416/W);
	  }
	}
      }//loop through hadrons
      
      //test lower E 
      
      neta_W[imc]->Fill(W,neta);
      if(nhch_f){
	neta_W_F[imc]->Fill(W,neta);
      }
      /* else{
	LOG("gvldtest", pERROR) 
	  << "nhch_f =  " << nhch_f;
	  }*/
      if (cut3){
	if(nhch_f2){
	  meanpt2_f/=nhch_f2;
	  pt2_W2_F[imc]->Fill(W2,meanpt2_f,weight);
	}
	if(nhch_b2){
	  meanpt2_b/=nhch_b2;
	  pt2_W2_B[imc]->Fill(W2,meanpt2_b,weight);
	}
      }
      if (cut5c){
	if (nhch&&W<10){
	  meanpt/=nhch;
	  //meanpt2/=nhch;
	  pt_W[imc]->Fill(W,meanpt,weight);
	  //	    if (meanpt>2||meanpt<0||weight!=1)
	  //	      cout<<W<<" "<<meanpt<<" "<<weight<<endl;
	}
	if (nhch_f){
	  meanpt_f/=nhch_f;
	  pt_W_F[imc]->Fill(W,meanpt_f,weight);
	}
	if (nhch_b){
	  meanpt_b/=nhch_b;
	  pt_W_B[imc]->Fill(W,meanpt_b,weight);
	}
      }	    
    }//loop over all events

    for (int i = 0; i<nw; i++){
      if (nv[i]){
	nchf[i]/=nv[i];
	nchb[i]/=nv[i];
      }
      if (cut5_nv[i]){
	aW[i]/=cut5_nv[i];
	nposf[i]/=cut5_nv[i];
	nnegf[i]/=cut5_nv[i];
	nposb[i]/=cut5_nv[i];
	nnegb[i]/=cut5_nv[i];
      }
    }

    xf_pip[imc]->Scale(xf_pip[imc]->GetNbinsX()/
		  (xf_pip[imc]->GetXaxis()->GetXmax()-xf_pip[imc]->GetXaxis()->GetXmin())/cut1_nv);
    
    xf_pim[imc]->Scale(xf_pim[imc]->GetNbinsX()/
		  (xf_pim[imc]->GetXaxis()->GetXmax()-xf_pim[imc]->GetXaxis()->GetXmin())/cut1_nv);
    fxf_pip[imc]->Scale(fxf_pip[imc]->GetNbinsX()/
		   (fxf_pip[imc]->GetXaxis()->GetXmax()-fxf_pip[imc]->GetXaxis()->GetXmin())/cut1_nv);
    
    fxf_pim[imc]->Scale(fxf_pim[imc]->GetNbinsX()/
		   (fxf_pim[imc]->GetXaxis()->GetXmax()-fxf_pim[imc]->GetXaxis()->GetXmin())/cut1_nv);
    z_pos[imc]->Scale(z_pos[imc]->GetNbinsX()/
		 (z_pos[imc]->GetXaxis()->GetXmax()-z_pos[imc]->GetXaxis()->GetXmin())/cut2_nv);
    z_neg[imc]->Scale(z_neg[imc]->GetNbinsX()/
		 (z_neg[imc]->GetXaxis()->GetXmax()-z_neg[imc]->GetXaxis()->GetXmin())/cut2_nv);
    z1_pos[imc]->Scale(z1_pos[imc]->GetNbinsX()/
		  (z1_pos[imc]->GetXaxis()->GetXmax()-z1_pos[imc]->GetXaxis()->GetXmin())/cut5a_nv);
    z1_pro[imc]->Scale(z1_pro[imc]->GetNbinsX()/
		  (z1_pro[imc]->GetXaxis()->GetXmax()-z1_pro[imc]->GetXaxis()->GetXmin())/cut5a_nv);
    z1_neg[imc]->Scale(z1_neg[imc]->GetNbinsX()/
		  (z1_neg[imc]->GetXaxis()->GetXmax()-z1_neg[imc]->GetXaxis()->GetXmin())/cut5a_nv);
    z2_pos[imc]->Scale(z2_pos[imc]->GetNbinsX()/
		  (z2_pos[imc]->GetXaxis()->GetXmax()-z2_pos[imc]->GetXaxis()->GetXmin())/cut5b_nv);
    z2_pro[imc]->Scale(z2_pro[imc]->GetNbinsX()/
		  (z2_pro[imc]->GetXaxis()->GetXmax()-z2_pro[imc]->GetXaxis()->GetXmin())/cut5b_nv);
    z2_neg[imc]->Scale(z2_neg[imc]->GetNbinsX()/
		  (z2_neg[imc]->GetXaxis()->GetXmax()-z2_neg[imc]->GetXaxis()->GetXmin())/cut5b_nv);
    
    z_E1[imc]->Scale(20./nv_E[0]);
    z_E2[imc]->Scale(20./nv_E[1]);
    z_E3[imc]->Scale(20./nv_E[2]);
    
    pt2_pip[imc]->Scale(20./(nv_E[0]+nv_E[1]+nv_E[2]));
    pt2_pim[imc]->Scale(20./(nv_E[0]+nv_E[1]+nv_E[2]));
    
    pt2_xf1[imc]->Scale(pt2_xf1[imc]->GetNbinsX()/
		   (pt2_xf1[imc]->GetXaxis()->GetXmax()-pt2_xf1[imc]->GetXaxis()->GetXmin())/cut4_nv);
    pt2_xf2[imc]->Scale(pt2_xf2[imc]->GetNbinsX()/
		   (pt2_xf2[imc]->GetXaxis()->GetXmax()-pt2_xf2[imc]->GetXaxis()->GetXmin())/cut4_nv);
    pt2_xf3[imc]->Scale(pt2_xf3[imc]->GetNbinsX()/
		   (pt2_xf3[imc]->GetXaxis()->GetXmax()-pt2_xf3[imc]->GetXaxis()->GetXmin())/cut4_nv);
    pt2_xf4[imc]->Scale(pt2_xf4[imc]->GetNbinsX()/
		   (pt2_xf4[imc]->GetXaxis()->GetXmax()-pt2_xf4[imc]->GetXaxis()->GetXmin())/cut4_nv);
    Fxf_pos1[imc]->Scale(Fxf_pos1[imc]->GetNbinsX()/
		    (Fxf_pos1[imc]->GetXaxis()->GetXmax()-Fxf_pos1[imc]->GetXaxis()->GetXmin())/cut5d_nv);
    Fxf_pos_kno1[imc]->Scale(Fxf_pos_kno1[imc]->GetNbinsX()/
			(Fxf_pos_kno1[imc]->GetXaxis()->GetXmax()-Fxf_pos_kno1[imc]->GetXaxis()->GetXmin())/cut5d_nv);
    Fxf_pro1[imc]->Scale(Fxf_pro1[imc]->GetNbinsX()/
		    (Fxf_pro1[imc]->GetXaxis()->GetXmax()-Fxf_pro1[imc]->GetXaxis()->GetXmin())/cut5d_nv);
    Fxf_neg1[imc]->Scale(Fxf_neg1[imc]->GetNbinsX()/
		    (Fxf_neg1[imc]->GetXaxis()->GetXmax()-Fxf_neg1[imc]->GetXaxis()->GetXmin())/cut5d_nv);
    Fxf_neg_kno1[imc]->Scale(Fxf_neg_kno1[imc]->GetNbinsX()/
			(Fxf_neg_kno1[imc]->GetXaxis()->GetXmax()-Fxf_neg_kno1[imc]->GetXaxis()->GetXmin())/cut5d_nv);
    Fxf_pos2[imc]->Scale(Fxf_pos2[imc]->GetNbinsX()/
		    (Fxf_pos2[imc]->GetXaxis()->GetXmax()-Fxf_pos2[imc]->GetXaxis()->GetXmin())/cut5e_nv);
    Fxf_neg2[imc]->Scale(Fxf_neg2[imc]->GetNbinsX()/
		    (Fxf_neg2[imc]->GetXaxis()->GetXmax()-Fxf_neg2[imc]->GetXaxis()->GetXmin())/cut5e_nv);
    Fxf_pos1_hi[imc]->Scale(Fxf_pos1_hi[imc]->GetNbinsX()/
		       (Fxf_pos1_hi[imc]->GetXaxis()->GetXmax()-Fxf_pos1_hi[imc]->GetXaxis()->GetXmin())/cut5f_nv);
    Fxf_pro1_hi[imc]->Scale(Fxf_pro1_hi[imc]->GetNbinsX()/
		       (Fxf_pro1_hi[imc]->GetXaxis()->GetXmax()-Fxf_pro1_hi[imc]->GetXaxis()->GetXmin())/cut5f_nv);
    Fxf_neg1_hi[imc]->Scale(Fxf_neg1_hi[imc]->GetNbinsX()/
		       (Fxf_neg1_hi[imc]->GetXaxis()->GetXmax()-Fxf_neg1_hi[imc]->GetXaxis()->GetXmin())/cut5f_nv);
    Fxf_pos2_hi[imc]->Scale(Fxf_pos2_hi[imc]->GetNbinsX()/
		       (Fxf_pos2_hi[imc]->GetXaxis()->GetXmax()-Fxf_pos2_hi[imc]->GetXaxis()->GetXmin())/cut5g_nv);
    Fxf_neg2_hi[imc]->Scale(Fxf_neg2_hi[imc]->GetNbinsX()/
		       (Fxf_neg2_hi[imc]->GetXaxis()->GetXmax()-Fxf_neg2_hi[imc]->GetXaxis()->GetXmin())/cut5g_nv);
    
    nch_w[imc] = new TGraphErrors(nw-1,aW2,nch,gerrx,gerrx);
    nchpi_w[imc] = new TGraph(nw-1,aW2,nchpi);
    
    D_nneg[imc] = new TGraph(nw-1,nneg,Dneg);
    
    D_W2[imc] = new TGraph(nw-1,aW2,D_nch);
    
    npi0_w[imc] = new TGraphErrors(nw-1,aW2,npizero,gerrx,errnpi0);
    
    D_npi0[imc] = new TGraph(nw-1,npizero,Dnpi0);
    
    for (int i = 0; i<13; i++){
      pnch_w2[i][imc] = new TGraph(nw-1,aW2,pnch[i]);
    }

    // KNO distributions
    for (int i = 0; i<5; i++){
      kno[i][imc]->Sumw2();
      kno[i][imc]->Scale(1/kno[i][imc]->Integral());
      kno[i][imc]->Scale(nchkno[i]);
      kno[i][imc]->SetMarkerColor(2);
      kno[i][imc]->SetMarkerSize(0.8);
      kno[i][imc]->SetLineColor(2);
      kno[i][imc]->SetStats(0);
    }    
    kno[0][imc]->SetMarkerStyle(20);
    kno[1][imc]->SetMarkerStyle(21);
    kno[2][imc]->SetMarkerStyle(22);
    kno[3][imc]->SetMarkerStyle(23);
    kno[4][imc]->SetMarkerStyle(29);        

    nch_w_f[imc]   = new TGraph(nw-1,aW2,nchf);
    nch_w_b[imc]   = new TGraph(nw-1,aW2,nchb);    
    npos_w_f[imc]  = new TGraph(nw-2,aW,nposf);
    nneg_w_f[imc]  = new TGraph(nw-2,aW,nnegf);
    npos_w_b[imc]  = new TGraph(nw-2,aW,nposb);
    nneg_w_b[imc]  = new TGraph(nw-2,aW,nnegb);    
    npos_w2_f[imc] = new TGraph(nw-2,aW2,nposf);
    nneg_w2_f[imc] = new TGraph(nw-2,aW2,nnegf);
    npos_w2_b[imc] = new TGraph(nw-2,aW2,nposb);
    nneg_w2_b[imc] = new TGraph(nw-2,aW2,nnegb);

  }//loop over all mc files
  
}
//____________________________________________________________________________
