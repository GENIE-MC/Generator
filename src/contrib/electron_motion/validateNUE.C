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

void make_dists(TString filename = "gntp.0.ghep.root", TString dir="",bool print_to_pdf=false,
  bool print_summary=false){
  TFile infile(dir+filename);
  //TFile infile('gntp.0.ghep.root');

  // Get the GENIE GHEP tree and set its branch address
  TTree * tree = dynamic_cast<TTree*> ( infile.Get("gtree") );

  std::string title_car = "Velocity Distributions";
  const char *char_title_car = title_car.c_str();

  std::string title_ang = "Angular Distributions";
  const char *char_title_ang = title_ang.c_str();

  std::string title_y = "y Distribution";
  const char *char_title_y = title_y.c_str();


  //SEt histograms
  TH1F *vx_dist = new TH1F("vx",char_title_car,20,-0.15,0.15);
  TH1F *vy_dist = new TH1F("vy",char_title_car,20,-0.15,0.15);
  TH1F *vz_dist = new TH1F("vz",char_title_car,20,-0.15,0.15);

  //2D histogram with cos theta and phi of velocities
  TH2F *ang_dist = new TH2F("ang",char_title_ang,20,-1,1,20,-3.1415,3.1415);

  //y distribution
  TH1F *y_dist = new TH1F("y",char_title_y,20,0,1);

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

      //Fill y distribution
      y_dist->Fill(inter.Kine().y(true));

      //Print out values
      if (print_summary){
        std::cout<<"Mass = "<<electron->Mass()<<std::endl;
        std::cout<<"px = "<<pe.Px()<< " py = "<<pe.Py()<<" pz = "<<pe.Pz()<<std::endl;
        std::cout<<"vx = "<<vx<< " vy = "<<vy<<" vz = "<<vz<<std::endl;
        std::cout<<"costheta = "<<costheta<<" phi = "<<phi<<std::endl;
        std::cout<<"y = "<<inter.Kine().y(true)<<std::endl;
      }
      

    }//If is scattering

  }//end event loop

  TFile *file1 = new TFile("Plots/histos.root","recreate");
  TCanvas canvas("canvas");
  vx_dist->Draw("HIST");
  vx_dist->SetLineColor(kBlack);
  //vx_dist->SetLineWidth
  vx_dist->SetFillColorAlpha(kBlue,0.35);
  vx_dist->GetXaxis()->SetTitle("v [c]");
  vx_dist->SetStats(kFALSE);
  vx_dist->Write();

  vy_dist->Draw("SAME");
  vy_dist->SetLineColor(kBlack);
  vy_dist->SetFillColorAlpha(kRed,0.35);
  vy_dist->GetXaxis()->SetTitle("v [c]");
  vy_dist->Write();

  vz_dist->Draw("SAME");
  vz_dist->SetLineColor(kBlack);
  vz_dist->SetFillColorAlpha(kGreen,0.35);
  vy_dist->GetXaxis()->SetTitle("v [c]");
  vz_dist->Write();

  if (print_to_pdf){
    auto legend = new TLegend(0.1,0.7,0.3,0.9);
    //legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
    legend->AddEntry(vx_dist,"vx");
    legend->AddEntry(vy_dist,"vy");
    legend->AddEntry(vz_dist,"vz");
    legend->Draw();

    std::string file_out = "Plots/vdist.pdf";
    const char *fo = file_out.c_str();

    canvas.Print(fo);
  }

  TCanvas canvas2("canvas2");
  ang_dist->Draw("COLZ");
  ang_dist->GetXaxis()->SetTitle("cos#theta");
  ang_dist->GetYaxis()->SetTitle("#phi");
  ang_dist->SetStats(kFALSE);
  ang_dist->Write();

  if (print_to_pdf){
    std::string file_out2 = "Plots/angdist.pdf";
    const char *fo2 = file_out2.c_str();
    canvas2.Print(fo2);
  }

  TCanvas canvas3("canvas3");
  y_dist->Draw("HIST");
  y_dist->GetXaxis()->SetTitle("y");
  y_dist->SetStats(kFALSE);
  y_dist->Write();

  if (print_to_pdf){
    std::string file_out3 = "Plots/ydist.pdf";
    const char *fo3 = file_out3.c_str();
    canvas3.Print(fo3);
  }
  
  file1->Close();
}

void validateNUE(TString filename = "gntp.0.ghep.root",TString dir=""){
  //Make all distributions
  gSystem->Exec("mkdir -p Plots");
  make_dists(filename,dir);
}