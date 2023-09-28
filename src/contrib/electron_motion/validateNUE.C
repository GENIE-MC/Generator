#include "TString.h" 
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"

#include "Framework/Ntuple/NtpMCTreeHeader.h" 
#include "Framework/Ntuple/NtpMCEventRecord.h" 
#include "Framework/EventGen/EventRecord.h" 
#include "Framework/ParticleData/BaryonResUtils.h"

#include "Framework/GHEP/GHepParticle.h"       
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/GHEP/GHepStatus.h"

#include "Framework/ParticleData/PDGCodes.h"
#include <string>
#include <iostream>
#include <cmath>

using namespace genie ;
using std::to_string;

double calc_velocity(double p, double m){
  //Calculate velocity from momentum
  return p/sqrt(pow(p,2)+pow(m,2));
}

void make_vdist(TString filename = "gntp.0.ghep.root", TString dir=""){
  TFile infile(dir+filename);
  //TFile infile('gntp.0.ghep.root');

  // Get the GENIE GHEP tree and set its branch address
  TTree * tree = dynamic_cast<TTree*> ( infile.Get("gtree") );

  std::string title_car = "Velocity Distributions";
  const char *char_title_car = title_car.c_str();

  std::string title_ang = "Angular Distributions";
  const char *char_title_ang = title_ang.c_str();


  //SEt histograms
  TH1F *vx_dist = new TH1F("vx",char_title_car,20,-0.15,0.15);
  TH1F *vy_dist = new TH1F("vy",char_title_car,20,-0.15,0.15);
  TH1F *vz_dist = new TH1F("vz",char_title_car,20,-0.15,0.15);

  //2D histogram with cos theta and phi of velocities
  TH2F *ang_dist = new TH2F("vxy",char_title_ang,20,-1,1,20,-3.1415,3.1415);

  NtpMCEventRecord * mcrec = 0;
  tree->SetBranchAddress( "gmcrec", &mcrec);

  // Event loop
  for(Long64_t i=0; i < tree->GetEntries(); i++){
    tree->GetEntry(i);

    EventRecord & event = *( mcrec->event );
    
    const Interaction & inter = *( event.Summary() ) ;

    const ProcessInfo & proc_info = inter.ProcInfo() ;

    if (proc_info.IsElectronScattering()){
      GHepParticle * electron = event.HitElectron(); //Find initial electron
      const TLorentzVector & pe = *( electron -> P4() ) ; //Get its 4 momentum
      //Get velocity from four momentum
      //double vx = calc_velocity(pe.Px(),electron->Mass());
      double vx = pe.Px()/electron->Energy();
      double vy = pe.Py()/electron->Energy();
      double vz = pe.Pz()/electron->Energy();

      //Get angles
      double costheta = vz/(sqrt(pow(vx,2)+pow(vy,2)+pow(vz,2)));
      double phi = vy/abs(vy)*acos(vx/(sqrt(pow(vx,2)+pow(vy,2))));

      //Fill histograms
      vx_dist->Fill(pe.Px());
      vy_dist->Fill(pe.Py());
      vz_dist->Fill(pe.Pz());

      //Fill 2d histogram with angles
      ang_dist->Fill(costheta,phi);

      //Print out values
      std::cout<<"Mass = "<<electron->Mass()<<std::endl;
      std::cout<<"px = "<<pe.Px()<< " py = "<<pe.Py()<<" pz = "<<pe.Pz()<<std::endl;
      std::cout<<"vx = "<<vx<< " vy = "<<vy<<" vz = "<<vz<<std::endl;
      std::cout<<"costheta = "<<costheta<<" phi = "<<phi<<std::endl;
      

    }//If is scattering

  }//end event loop

  //TFile *file1 = new TFile(dir+"/Plots/histos.root","recreate");
  TCanvas canvas("canvas");
  vx_dist->Draw("HIST");
  vx_dist->SetLineColor(kBlack);
  //vx_dist->SetLineWidth
  vx_dist->SetFillColorAlpha(kBlue,0.35);
  vx_dist->GetXaxis()->SetTitle("v [c]");
  vx_dist->SetStats(kFALSE);
  //vx_dist->Write();
  vy_dist->Draw("SAME");
  vy_dist->SetLineColor(kBlack);
  vy_dist->SetFillColorAlpha(kRed,0.35);
  //vy_dist->Write();
  vz_dist->Draw("SAME");
  vz_dist->SetLineColor(kBlack);
  vz_dist->SetFillColorAlpha(kGreen,0.35);
  //vz_dist->Write();
  //file1->Close();
  auto legend = new TLegend(0.1,0.7,0.3,0.9);
  //legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
  legend->AddEntry(vx_dist,"vx");
  legend->AddEntry(vy_dist,"vy");
  legend->AddEntry(vz_dist,"vz");
  legend->Draw();

  std::string file_out = "Plots/vdist.pdf";
  const char *fo = file_out.c_str();

  canvas.Print(fo);

  TCanvas canvas2("canvas2");
  ang_dist->Draw("COLZ");
  ang_dist->GetXaxis()->SetTitle("cos#theta");
  ang_dist->GetYaxis()->SetTitle("#phi");
  ang_dist->SetStats(kFALSE);
  std::string file_out2 = "Plots/angdist.pdf";
  const char *fo2 = file_out2.c_str();
  canvas2.Print(fo2);
}

void make_y_dist(TString filename = "gntp.0.ghep.root", TString dir=""){
  TFile infile(dir+filename);
  //TFile infile('gntp.0.ghep.root');

  // Get the GENIE GHEP tree and set its branch address
  TTree * tree = dynamic_cast<TTree*> ( infile.Get("gtree") );

  //SEt histograms
  TH1F *y_dist = new TH1F("y","y",20,0,1);

  NtpMCEventRecord * mcrec = 0;
  tree->SetBranchAddress( "gmcrec", &mcrec);

  // Event loop
  for(Long64_t i=0; i < tree->GetEntries(); i++){
    tree->GetEntry(i);

    EventRecord & event = *( mcrec->event );
    
    const Interaction & inter = *( event.Summary() ) ;

    const ProcessInfo & proc_info = inter.ProcInfo() ;

    if (proc_info.IsElectronScattering()){
      GHepParticle * electron = event.HitElectron(); //Find initial electron
      const TLorentzVector & pe = *( electron -> P4() ) ; //Get its 4 momentum

      //Fill histograms
      y_dist->Fill(inter.Kine().y(true));

    }//If is scattering

  }//end event loop

  TCanvas canvas("canvas");
  y_dist->Draw("HIST");
  y_dist->SetLineColor(kBlack);
  y_dist->SetFillColorAlpha(kBlue,0.85);
  y_dist->GetXaxis()->SetTitle("y");
  y_dist->SetStats(kFALSE);
  std::string file_out = "Plots/ydist.pdf";
  const char *fo = file_out.c_str();
  canvas.Print(fo);
}

void validateNUE(TString filename = "gntp.0.ghep.root",TString dir=""){
  //Make all distributions
  gSystem->Exec("mkdir -p Plots");
  make_vdist(filename,dir);
  make_y_dist(filename,dir);
}

void make_finale_tree(TString filename = "gntp.0.ghep.root", TString dir=""){
  std::vector<int> nu_pdgs = {12,14,-12,-14};


  //Initial electron
  std::vector<Float_t> initiale_E;
  std::vector<Float_t> initiale_px;
  std::vector<Float_t> initiale_py;
  std::vector<Float_t> initiale_pz;

  //Final electron
  std::vector<Float_t> finale_E;
  std::vector<Float_t> finale_px;
  std::vector<Float_t> finale_py;
  std::vector<Float_t> finale_pz;

  //Initial nu
  std::vector<Int_t> initialnu_pdg;
  std::vector<Float_t> initialnu_E;
  std::vector<Float_t> initialnu_px;
  std::vector<Float_t> initialnu_py;
  std::vector<Float_t> initialnu_pz;

  //Final nu
  std::vector<Int_t> finalnu_pdg;
  std::vector<Float_t> finalnu_E;
  std::vector<Float_t> finalnu_px;
  std::vector<Float_t> finalnu_py;
  std::vector<Float_t> finalnu_pz;

  //Cross section
  std::vector<Float_t> xsecs;
  std::vector<Float_t> y;

  //Make TTree
  TFile* file = new TFile("nuescat.root", "RECREATE");
  TTree *nuetree = new TTree("nuescat","Nue Kinematics");

  //Electron
  nuetree->Branch("initiale_E",&initiale_E);
  nuetree->Branch("initiale_px",&initiale_px);
  nuetree->Branch("initiale_py",&initiale_py);
  nuetree->Branch("initiale_pz",&initiale_pz);

  nuetree->Branch("finale_E",&finale_E);
  nuetree->Branch("finale_px",&finale_px);
  nuetree->Branch("finale_py",&finale_py);
  nuetree->Branch("finale_pz",&finale_pz);

  //Nu
  nuetree->Branch("initialnu_pdg",&initialnu_pdg);
  nuetree->Branch("initialnu_E",&initialnu_E);
  nuetree->Branch("initialnu_px",&initialnu_px);
  nuetree->Branch("initialnu_py",&initialnu_py);
  nuetree->Branch("initialnu_pz",&initialnu_pz);

  nuetree->Branch("finalnu_pdg",&finalnu_pdg);
  nuetree->Branch("finalnu_E",&finalnu_E);
  nuetree->Branch("finalnu_px",&finalnu_px);
  nuetree->Branch("finalnu_py",&finalnu_py);
  nuetree->Branch("finalnu_pz",&finalnu_pz);

  //Cross section
  nuetree->Branch("xsecs",&xsecs);
  nuetree->Branch("y",&y);

  TFile infile(dir+filename);
  //TFile infile('gntp.0.ghep.root');

  // Get the GENIE GHEP tree and set its branch address
  TTree * tree = dynamic_cast<TTree*> ( infile.Get("gtree") );

  NtpMCEventRecord * mcrec = 0;
  tree->SetBranchAddress( "gmcrec", &mcrec);

  // Event loop
  for(Long64_t i=0; i < tree->GetEntries(); i++){
    tree->GetEntry(i);

    EventRecord & event = *( mcrec->event );
    
    const Interaction & inter = *( event.Summary() ) ;
    y.push_back(inter.Kine().y());

    const ProcessInfo & proc_info = inter.ProcInfo() ;

    if (proc_info.IsElectronScattering()){
      GHepParticle * electron = event.HitElectron(); //Find initial electron
      const TLorentzVector & pe = *( electron -> P4() ) ; //Get its 4 momentum

      GHepParticle * finale = event.FinalStatePrimaryLepton(); //Find final electron
      const TLorentzVector & finalpe = *(finale -> P4()); //Get its 4 momentum

      int init_nupdg = 0;
      int final_nupdg = 0;

      // GHepStatus_t initialstatus = GHepStatus(0);
      // GHepStatus_t finalstatus = GHepStatus(1);
      // for (int pdg : nu_pdgs){
      //   if (event.FindParticle(pdg,initialstatus,0) != 0){
      //     GHepParticle * nu = event.FindParticle(pdg,0,0);
      //     init_nupdg = pdg;
      //   }
      //   if (event.FindParticle(pdg,finalstatus,0) != 0){
      //     GHepParticle * finalnu = event.FindParticle(pdg,0,0);
      //     final_nupdg = pdg;
      //   }
      GHepParticle * nu = event.Particle(0); //Initial nu
      const TLorentzVector & pnu = *( nu -> P4() ) ; //Get its 4 momentum

      GHepParticle * finalnu = event.Particle(5); //Final nu
      const TLorentzVector & finalpnu = *( finalnu -> P4() ) ; //Get its 4 momentum
      

      //Fill vectors
      initiale_E.push_back(pe.E());
      initiale_px.push_back(pe.Px());
      initiale_py.push_back(pe.Py());
      initiale_pz.push_back(pe.Pz());

      finale_E.push_back(finalpe.E());
      finale_px.push_back(finalpe.Px());
      finale_py.push_back(finalpe.Py());
      finale_pz.push_back(finalpe.Pz());

      //Fill with dumby values for now
      initialnu_pdg.push_back(nu->Pdg());
      initialnu_E.push_back(pnu.E());
      initialnu_px.push_back(pnu.Px());
      initialnu_py.push_back(pnu.Py());
      initialnu_pz.push_back(pnu.Pz());

      finalnu_pdg.push_back(finalnu->Pdg());
      finalnu_E.push_back(finalpnu.E());
      finalnu_px.push_back(finalpnu.Px());
      finalnu_py.push_back(finalpnu.Py());
      finalnu_pz.push_back(finalpnu.Pz());

      //Fill xsec
      xsecs.push_back(event.XSec());

    }//If is scattering
  }//end event loop
  nuetree->Fill();
  file->Write();


}