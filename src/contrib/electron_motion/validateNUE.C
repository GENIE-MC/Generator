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

  std::string title_mom = "Momentum Distributions";
  const char *char_title_mom = title_mom.c_str();

  std::string title_ang = "Angular Distributions";
  const char *char_title_ang = title_ang.c_str();

  std::string title_y = "y Distribution";
  const char *char_title_y = title_y.c_str();

  std::string title_xsec = "xsec";
  const char *char_title_xsec = title_xsec.c_str();


  //SEt histograms
  TH1F *vx_dist = new TH1F("vx",char_title_car,20,-0.15,0.15);
  TH1F *vy_dist = new TH1F("vy",char_title_car,20,-0.15,0.15);
  TH1F *vz_dist = new TH1F("vz",char_title_car,20,-0.15,0.15);

  TH1F *px_dist = new TH1F("px",char_title_mom,20,-1e-4,1e-4);
  TH1F *py_dist = new TH1F("py",char_title_mom,20,-1e-4,1e-4);
  TH1F *pz_dist = new TH1F("pz",char_title_mom,20,-1e-4,1e-4);

  //2D histogram with cos theta and phi of velocities
  TH2F *ang_dist = new TH2F("ang",char_title_ang,20,-1,1,20,-3.1415,3.1415);

  //y distribution
  TH1F *y_dist = new TH1F("y",char_title_y,20,0,1);

  //xsec distribution
  TH1F *xsec_dist = new TH1F("xsec",char_title_xsec,20,1e-13,1e-16);

  //Scatter plot of E_l theta_e^2 vs theta_l
  std::vector<double> E_l; //final lepton energy
  std::vector<double> Etheta2_l; //final lepton etheta^2
  std::vector<double> theta_l; //final lepton angle

  //Scatter plot of E_l theta_e^2 vs theta_l for final state lepton in electron rest frame
  std::vector<double> E_l_erf; //final lepton energy in electron rest frame
  std::vector<double> Etheta2_l_erf; //final lepton etheta^2 in electron rest frame
  std::vector<double> theta_l_erf; //final lepton angle in electron rest frame

  //Differences between distributions in electron rest frame and lab frame
  std::vector<double> E_l_diff; //final lepton energy
  std::vector<double> Etheta2_l_diff; //final lepton etheta^2
  std::vector<double> theta_l_diff; //final lepton angle

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
      GHepParticle * fs_lepton = event.FinalStatePrimaryLepton(); //Find final state lepton
      GHepParticle * init_neutrino = event.Probe(); //Find initial neutrino

      const TLorentzVector & pe = *( electron -> P4() ) ; //Get its 4 momentum
      const TLorentzVector & ple = *( fs_lepton -> P4() ) ; //Get final state lepton 4 momentum
      
      //Get velocity from four momentum
      double vx = pe.Px()/electron->Energy();
      double vy = pe.Py()/electron->Energy();
      double vz = pe.Pz()/electron->Energy();

      double fs_vx = ple.Px()/fs_lepton->Energy();
      double fs_vy = ple.Py()/fs_lepton->Energy();
      double fs_vz = ple.Pz()/fs_lepton->Energy();

      //Get angles
      double costheta = vz/(sqrt(pow(vx,2)+pow(vy,2)+pow(vz,2)));
      double phi = vy/abs(vy)*acos(vx/(sqrt(pow(vx,2)+pow(vy,2))));

      double fs_costheta = fs_vz/(sqrt(pow(fs_vx,2)+pow(fs_vy,2)+pow(fs_vz,2)));

      //Fill vectors
      E_l.push_back(ple.E());
      theta_l.push_back(acos(fs_costheta));
      Etheta2_l.push_back(pow(acos(fs_costheta),2)*ple.E());

      //Fill histograms
      vx_dist->Fill(vx);
      vy_dist->Fill(vy);
      vz_dist->Fill(vz);

      //Fill momentum distributions
      px_dist->Fill(pe.Px());
      py_dist->Fill(pe.Py());
      pz_dist->Fill(pe.Pz());

      //Fill 2d histogram with angles
      ang_dist->Fill(costheta,phi);

      //Fill y distribution
      y_dist->Fill(inter.Kine().y(true));

      //Fill xsec distribution
      xsec_dist->Fill(event.XSec()/init_neutrino->Energy()); //normalize to neutrino energy

      //Boost final state lepton to initial electron rest frame
      TLorentzVector ple_boosted = ple;
      ple_boosted.Boost(-pe.BoostVector());

      //Extract values for final state lepton in electron rest frame
      double E_l_erf_val = ple_boosted.E();
      double fs_costheta_erf_val = ple_boosted.Pz()/ple_boosted.P();
      double Etheta2_l_erf_val = pow(acos(fs_costheta_erf_val),2)*ple_boosted.E();

      //Fill vectors
      E_l_erf.push_back(E_l_erf_val);
      theta_l_erf.push_back(acos(fs_costheta_erf_val));
      Etheta2_l_erf.push_back(Etheta2_l_erf_val);

      //Extract values for final state lepton in electron rest frame and compare to lab frame
      double E_l_diff_val = E_l_erf_val - ple.E();
      double fs_costheta_diff_val = fs_costheta_erf_val - fs_costheta;
      double Etheta2_l_diff_val = Etheta2_l_erf_val - Etheta2_l[i];

      //Fill vectors
      E_l_diff.push_back(E_l_diff_val);
      theta_l_diff.push_back(fs_costheta_diff_val);
      Etheta2_l_diff.push_back(Etheta2_l_diff_val);

      //Print out values
      if (print_summary){
        std::cout<<"Event "<<i<<std::endl;
        std::cout<<"Mass = "<<electron->Mass()<<std::endl;
        std::cout<<"px = "<<pe.Px()<< " py = "<<pe.Py()<<" pz = "<<pe.Pz()<<std::endl;
        std::cout<<"vx = "<<vx<< " vy = "<<vy<<" vz = "<<vz<<std::endl;
        std::cout<<"costheta = "<<costheta<<" phi = "<<phi<<std::endl;
        std::cout<<"y = "<<inter.Kine().y(true)<<std::endl;
        std::cout<<"E_l = "<<ple.E()<<std::endl;
        std::cout<<"theta_l = "<<acos(fs_costheta)<<std::endl;
        std::cout<<"xsec = "<<event.XSec()<<std::endl;
        std::cout<<"E_nu = "<<init_neutrino->Energy()<<std::endl;
        std::cout<<"E_l_erf = "<<E_l_erf_val<<std::endl;
        std::cout<<"theta_l_erf = "<<acos(fs_costheta_erf_val)<<std::endl;
        std::cout<<"E_l_diff = "<<E_l_diff_val<<std::endl;
        std::cout<<"theta_l_diff = "<<fs_costheta_diff_val<<std::endl;
        std::cout<<"---------------------------------"<<std::endl;
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
  ang_dist->Write();

  if (print_to_pdf){
    std::string file_out2 = "Plots/angdist.pdf";
    const char *fo2 = file_out2.c_str();
    canvas2.Print(fo2);
  }

  TCanvas canvas3("canvas3");
  y_dist->Draw("HIST");
  y_dist->GetXaxis()->SetTitle("y");
  // Set y limit to be from 0 to 600
  y_dist->SetMinimum(0);
  y_dist->SetMaximum(600);
  y_dist->Write();

  if (print_to_pdf){
    std::string file_out3 = "Plots/ydist.pdf";
    const char *fo3 = file_out3.c_str();
    canvas3.Print(fo3);
  }

  //Scatter plot of E_l theta_e^2 vs theta_l
  TCanvas canvas4("canvas4");
  TGraph *graph = new TGraph(theta_l.size(),&theta_l[0],&Etheta2_l[0]);
  graph->SetTitle("Etheta2");
  graph->SetMarkerStyle(20);
  graph->SetMarkerSize(0.5);
  graph->SetLineWidth(0);
  graph->Draw("AP");
  graph->GetXaxis()->SetTitle("#theta_{l}");
  graph->GetYaxis()->SetTitle("E_{l}#theta_{l}^{2}");
  graph->Write();

  if (print_to_pdf){
    std::string file_out4 = "Plots/etheta2.pdf";
    const char *fo4 = file_out4.c_str();
    canvas4.Print(fo4);
  }
  //Scatter plot of E_l vs theta_l
  TCanvas canvas6("canvas6");
  TGraph *g2 = new TGraph(theta_l.size(),&theta_l[0],&E_l[0]);
  g2->SetTitle("El");
  g2->SetMarkerStyle(20);
  g2->SetMarkerSize(0.5);
  g2->SetLineWidth(0);
  g2->Draw("AP");
  g2->GetXaxis()->SetTitle("#theta_{l}");
  g2->GetYaxis()->SetTitle("E_{l}");
  g2->Write();

  if (print_to_pdf){
    std::string file_out6 = "Plots/etheta.pdf";
    const char *fo6 = file_out6.c_str();
    canvas6.Print(fo6);
  }

  //Momentum distribution
  TCanvas canvas5("canvas5");
  px_dist->Draw("HIST");
  px_dist->SetLineColor(kBlack);
  //px_dist->SetLineWidth
  px_dist->SetFillColorAlpha(kBlue,0.35);
  px_dist->GetXaxis()->SetTitle("p [c]");
  px_dist->SetStats(kFALSE);
  px_dist->Write();

  py_dist->Draw("SAME");
  py_dist->SetLineColor(kBlack);
  py_dist->SetFillColorAlpha(kRed,0.35);
  py_dist->GetXaxis()->SetTitle("p [c]");
  py_dist->Write();

  pz_dist->Draw("SAME");
  pz_dist->SetLineColor(kBlack);
  pz_dist->SetFillColorAlpha(kGreen,0.35);
  py_dist->GetXaxis()->SetTitle("p [c]");
  pz_dist->Write();

  if (print_to_pdf){
    auto legend = new TLegend(0.1,0.7,0.3,0.9);
    //legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
    legend->AddEntry(px_dist,"px");
    legend->AddEntry(py_dist,"py");
    legend->AddEntry(pz_dist,"pz");
    legend->Draw();

    std::string file_out5 = "Plots/pdist.pdf";
    const char *fo5 = file_out5.c_str();

    canvas5.Print(fo5);
  }

  TCanvas canvas7("canvas7");
  xsec_dist->Draw("HIST");
  xsec_dist->GetXaxis()->SetTitle("xsec/E_{#nu}");
  xsec_dist->Write();
  if (print_to_pdf){
    std::string file_out7 = "Plots/xsec.pdf";
    const char *fo7 = file_out7.c_str();
    canvas7.Print(fo7);
  }

  //Scatter plot of E_l vs theta_l for final state lepton in electron rest frame
  TCanvas canvas8("canvas8");
  TGraph *g3 = new TGraph(theta_l_erf.size(),&theta_l_erf[0],&E_l_erf[0]);
  g3->SetTitle("El (electron rest frame)");
  g3->SetMarkerStyle(20);
  g3->SetMarkerSize(0.5);
  g3->SetLineWidth(0);
  g3->Draw("AP");
  g3->GetXaxis()->SetTitle("#theta_{l}");
  g3->GetYaxis()->SetTitle("E_{l}");
  g3->Write();

  if (print_to_pdf){
    std::string file_out8 = "Plots/etheta_erf.pdf";
    const char *fo8 = file_out8.c_str();
    canvas8.Print(fo8);
  }

  //Scatter plot of E_l theta_l^2 vs theta_l for final state lepton in electron rest frame
  TCanvas canvas9("canvas9");
  TGraph *g4 = new TGraph(theta_l_erf.size(),&theta_l_erf[0],&Etheta2_l_erf[0]);
  g4->SetTitle("Etheta2 (electron rest frame)");
  g4->SetMarkerStyle(20);
  g4->SetMarkerSize(0.5);
  g4->SetLineWidth(0);
  g4->Draw("AP");
  g4->GetXaxis()->SetTitle("#theta_{l}");
  g4->GetYaxis()->SetTitle("E_{l}#theta_{l}^{2}");
  g4->Write();
  if (print_to_pdf){
    std::string file_out9 = "Plots/etheta2_erf.pdf";
    const char *fo9 = file_out9.c_str();
    canvas9.Print(fo9);
  }

  //Differences between distributions in electron rest frame and lab frame
  TCanvas canvas10("canvas10");
  TGraph *g5 = new TGraph(theta_l_diff.size(),&theta_l_diff[0],&E_l_diff[0]);
  g5->SetTitle("El (electron rest frame - lab frame)");
  g5->SetMarkerStyle(20);
  g5->SetMarkerSize(0.5);
  g5->SetLineWidth(0);
  g5->Draw("AP");
  g5->GetXaxis()->SetTitle("#theta_{l}");
  g5->GetYaxis()->SetTitle("E_{l}");
  g5->Write();
  if (print_to_pdf){
    std::string file_out10 = "Plots/etheta_diff.pdf";
    const char *fo10 = file_out10.c_str();
    canvas10.Print(fo10);
  }

  TCanvas canvas11("canvas11");
  TGraph *g6 = new TGraph(theta_l_diff.size(),&theta_l_diff[0],&Etheta2_l_diff[0]);
  g6->SetTitle("Etheta2 (electron rest frame - lab frame)");
  g6->SetMarkerStyle(20);
  g6->SetMarkerSize(0.5);
  g6->SetLineWidth(0);
  g6->Draw("AP");
  g6->GetXaxis()->SetTitle("#theta_{l}");
  g6->GetYaxis()->SetTitle("E_{l}#theta_{l}^{2}");
  g6->Write();
  if (print_to_pdf){
    std::string file_out11 = "Plots/etheta2_diff.pdf";
    const char *fo11 = file_out11.c_str();
    canvas11.Print(fo11);
  }
  
  file1->Close();
}

void validateNUE(TString filename = "gntp.0.ghep.root",TString dir="",bool print_to_pdf=false,
  bool print_summary=false){
  //Make all distributions
  gSystem->Exec("mkdir -p Plots");
  make_dists(filename,dir,print_to_pdf,print_summary);
}