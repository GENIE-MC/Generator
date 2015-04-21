// *********************************************************************
//  Generate 4D diff-xsec histograms - Martti Nirkko (28th Nov 2014)
//  Compile and run in terminal:     root -l -q d4sigma_hist.C+g
// *********************************************************************

#include "code/singlekaon_xsec.cxx"
#include <TMath.h>
#include <TFile.h>
#include <TH3D.h>

// Compile using:   root -l -b -q d4sigma_hist.C+g
void d4sigma_hist() {
  
  // ***************************
  // *   STEERING PARAMETERS   *
  // ***************************
  const int COMP = 10;                  // time complexity (runs with O(n^4)...)
  const int pLeptonIsUsed = 0;          // use momentum instead of kinetic energy
  const int thetaIsUsed = 0;            // use angle instead of cosine
  const double pi = TMath::Pi();
  
  // NOTE: All results for a changing "pleptonIsUsed" and "thetaIsUsed" give the same results within
  //       about 1% regarding the total cross-section. Here, for making histograms of the differential
  //       1D cross-section, I use the variables I want to actually plot (Tlep and costheta)
  
  int type = 2;                         // lepton:   1=electron, 2=muon, 3=tau
  int reac = 3;                         // reaction: 1=NN, 2=NP, 3=PP
  
  double Enu;                           // neutrino energy [GeV]
  printf("Please enter neutrino energy: ");
  scanf("%lf", &Enu);
  printf("E_nu set to %.3lf GeV...\n", Enu);
  singlekaon_xsec *inst = new singlekaon_xsec();
  inst->init(Enu, type, reac);
  
  double Mlep = 0.;                     // lepton mass [GeV]
  if      (type==1) Mlep = 0.510998928e-3;
  else if (type==2) Mlep = 0.1056583715;
  else if (type==3) Mlep = 1.77682;
  else {std::cout<<"ERROR: Invalid type!"<<std::endl; return;}
  
  double Mkaon = 0.;                    // kaon mass [GeV]
  if      (reac==1) Mkaon = 0.493677;
  else if (reac==2) Mkaon = 0.497614;
  else if (reac==3) Mkaon = 0.493677;
  else {std::cout<<"ERROR: Invalid reaction!"<<std::endl; return;}
  
  // Initialisation of variables
  int i,j,k,l;
  const int nsteps1 = 10*COMP, nsteps2 = 10*COMP, nsteps3 = 10*COMP, nsteps4 = 2*COMP;
  double varmin1, varmin2, varmin3, varmin4, varmax1, varmax2, varmax3, varmax4, temp;
  double binsvar1[nsteps1+1], binsvar2[nsteps2+1], binsvar3[nsteps3+1], binsvar4[nsteps4+1];
  double varvals1[nsteps1],   varvals2[nsteps2],   varvals3[nsteps3],   varvals4[nsteps4];
  
  double Tlep, Elep;                    // LEPTON ENERGY [GeV]
  double Tkaon;                         // KAON ENERGY [GeV]
  double theta;                         // LEPTON ANGLE [rad]
  double phikq;                         // KAON ANGLE [rad]
  double diff4sigma;                    // DIFFERENTIAL CROSS SECTION [cm^2/GeV^2/rad^2]
  
  // Integration ranges
  varmin1 = 0.; varmax1 = Enu-Mkaon-Mlep;     // Tkaon
  if (pLeptonIsUsed) {                        // plep
    varmin2 = 0.;
    varmax2 = sqrt((Enu-Mkaon)*(Enu-Mkaon)-Mlep*Mlep);
  } else {                                    // Tlep
    varmin2 = 0.;
    varmax2 = Enu-Mkaon-Mlep;
  }
  if (thetaIsUsed) {                          // theta
    varmin3 = 0.;
    varmax3 = pi;
  } else {                                    // cos(theta)
    varmin3 = -1.0;
    varmax3 = 1.0;
  }
  varmin4 = -pi; varmax4 = pi;                // phi_kq
  
  // Calculate edges of bins
  for (i=0; i<=nsteps1; i++) {
    binsvar1[i] = varmin1 + i*(varmax1-varmin1)/nsteps1;
  }
  for (j=0; j<=nsteps2; j++) {
    temp = varmin2 + j*(varmax2-varmin2)/nsteps2;   // plep OR Tlep
    if (pLeptonIsUsed) {
      binsvar2[j] = sqrt(temp*temp+Mlep*Mlep)-Mlep; // plep -> Tlep
    } else {
      binsvar2[j] = temp;                           // Tlep
    }
  }
  for (k=0; k<=nsteps3; k++) {
    if (thetaIsUsed) {
      temp = varmax3 - k*(varmax3-varmin3)/nsteps3;   // theta
      binsvar3[k] = cos(temp);                        // cos(theta)
    } else {
      binsvar3[k] = varmin3 + k*(varmax3-varmin3)/nsteps3;
    }
  }
  for (l=0; l<=nsteps4; l++) {
    binsvar4[l] = varmin4 + l*(varmax4-varmin4)/nsteps4;
  }
  
  // Calculate edges of bins (NEW / TEST)
  /*for (i=0; i<=nsteps1; i++) binsvar1[i] = varmin1 + i*(varmax1-varmin1)/nsteps1;
  for (j=0; j<=nsteps2; j++) binsvar2[j] = varmin2 + j*(varmax2-varmin2)/nsteps2;   // plep OR Tlep
  for (k=0; k<=nsteps3; k++) binsvar3[k] = varmin3 + k*(varmax3-varmin3)/nsteps3;   // theta OR costh
  for (l=0; l<=nsteps4; l++) binsvar4[l] = varmin4 + l*(varmax4-varmin4)/nsteps4;*/
  
  // Calculate central bin values (average of left and right edge of bin)
  for (i=0; i<nsteps1; i++) varvals1[i] = 0.5*(binsvar1[i+1]+binsvar1[i]);
  for (j=0; j<nsteps2; j++) varvals2[j] = 0.5*(binsvar2[j+1]+binsvar2[j]);
  for (k=0; k<nsteps3; k++) varvals3[k] = 0.5*(binsvar3[k+1]+binsvar3[k]);
  for (l=0; l<nsteps4; l++) varvals4[l] = 0.5*(binsvar4[l+1]+binsvar4[l]);
  
  TH3D *hist[nsteps4];
  for (l=0; l<nsteps4; l++) {
    hist[l] = new TH3D(Form("hist_%2d",l), "diff4sigma", nsteps1, binsvar1, nsteps2, binsvar2, nsteps3, binsvar3);
  }
  
  // CALCULATE CROSS-SECTION
  // -----------------------
  
  for (i=0; i<nsteps1; i++) {
    if (binsvar1[i+1]==binsvar1[i]) continue;
    Tkaon = varvals1[i];
    
    for (j=0; j<nsteps2; j++) {
      if (binsvar2[j+1]==binsvar2[j]) continue;
      Tlep = varvals2[j];
      Elep = Tlep+Mlep;
      
      for (k=0; k<nsteps3; k++) {
        if (binsvar3[k+1]==binsvar3[k]) continue;
        theta = acos(varvals3[k]);
        
        for (l=0; l<nsteps4; l++) {
          if (binsvar4[l+1]==binsvar4[l]) continue;
          phikq = varvals4[l];
          
          // Calculate 4D differential xsec
          diff4sigma = inst->diffxsec(Tlep,Tkaon,theta,phikq);
          diff4sigma *= 2*pi;
          
          // This Jacobian no longer needed, since binning is adjusted to cos(theta)
          // diff4sigma *= sin(theta);
          
          // Multiplication with Jacobian for transformation plep->Tlep (it is E/p)
          diff4sigma *= Elep/sqrt(Elep*Elep-Mlep*Mlep);
          
          hist[l]->SetBinContent(i+1,j+1,k+1,diff4sigma);
        }
        
      }
      
    }
    
    if ((i+1)%5==0)
      std::cout << "Processing... (" << 100.*(i+1.)/nsteps1 << "%)" << std::endl;
    
  }
  
  std::string fname = Form("data/d4sigma_hist_%3.1lfGeV.root", Enu);
  TFile* outfile = new TFile(fname.c_str(), "RECREATE");
  for (l=0; l<nsteps4; l++) hist[l]->Write(Form("d4sigma_hist_%d",l));
  outfile->Close();
  std::cout << std::endl << "Output written to file: " << fname << std::endl << std::endl;
  
}

