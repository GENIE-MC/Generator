// *********************************************************************
//  Xsec (4D) of single kaon production - Martti Nirkko (28th Nov 2014)
//  Compile and run in terminal:     root -l -q d4sigma_calc.C+g
// *********************************************************************

#include "code/singlekaon_xsec.cxx"
#include <TMath.h>
#include <TString.h>

// Compile using:   root -l -b -q d4sigma_calc.C+g
void d4sigma_calc() {
  
  // ***************************
  // *   STEERING PARAMETERS   *
  // ***************************
  const int verbose = 0;                // debugging mode on/off
  const int pLeptonIsUsed = 1;          // binning for momentum instead of kinetic energy
  const int thetaIsUsed = 1;            // binning for angle instead of cosine
  const double pi = TMath::Pi();
  
  // NOTE: All results for a changing "pleptonIsUsed" and "thetaIsUsed" give the same results within
  //       about 1%. Here, both are set to 1 because this is closest to what the Fortran code does.
  //       For making differential histograms, it may be convenient to use other combinations.
  
  const int type = 2;                   // lepton:   1=electron, 2=muon, 3=tau
  const int reac = 3;                   // reaction: 1=NN, 2=NP, 3=PP
  
  int COMP;                             // time complexity (runs with O(n^4))
  printf("Please enter time complexity: ");
  scanf("%d", &COMP);
  printf("COMP set to %d...\n", COMP);
  
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
  
  double Enu = 0.;                      // neutrino energy [GeV]
  singlekaon_xsec *inst = new singlekaon_xsec();
  inst->init(Enu, type, reac);
  
  double Emin = inst->threshold;
  double Emax = 3.1;
  int Esteps = (int)(20.0*(Emax-Emin));
  
  // Initialisation of variables
  int i,j,k,l,itenu;
  const int nsteps1 = 10*COMP, nsteps2 = 10*COMP, nsteps3 = 10*COMP, nsteps4 = 2*COMP;
  double varmin1, varmin2, varmin3, varmin4, varmax1, varmax2, varmax3, varmax4, temp;
  double binsvar1[nsteps1+1], binsvar2[nsteps2+1], binsvar3[nsteps3+1], binsvar4[nsteps4+1];
  double varvals1[nsteps1],   varvals2[nsteps2],   varvals3[nsteps3],   varvals4[nsteps4];
  
  double Tlep,  Elep;                   // LEPTON ENERGY [GeV]
  double Tkaon, Ekaon;                  // KAON ENERGY [GeV]
  double theta;                         // LEPTON ANGLE [rad]
  double phikq;                         // KAON ANGLE [rad]
  
  double diff4sigma;
  double diff3sigma;
  double diff2sigma;
  double diff1sigma;
  double sigma;                         // CROSS SECTION [cm^2]
  
  std::cout << std::endl;
  std::cout << "  E [GeV]\txsec [cm^2]" << std::endl;
  std::cout << "-----------------------------" << std::endl;
  
  std::string fpath = "./data/";
  std::string fname = "d4sigma_calc.out";
  FILE *outfile = fopen((fpath+fname).c_str(),"w");
  fprintf(outfile, "#  E [GeV]\txsec [cm^2]\n");
  fprintf(outfile, "#-----------------------------\n");
  
  for (itenu=0; itenu<=Esteps; itenu++) {
    
    // Initialisation of energy
    Enu = Emin + itenu*0.05;
    inst->init(Enu, type, reac);
    
    // Integration ranges
    varmin1 = 0.; varmax1 = Enu-Mkaon-Mlep;     // Tkaon
    // varmin2, varmax2 filled in loop below
    if (thetaIsUsed) {                          // theta
      varmin3 = 0.;
      varmax3 = pi;
    } else {                                    // cos(theta)
      varmin3 = -1.0;
      varmax3 = 1.0;
    }
    varmin4 = 0.;  varmax4 = 2.*pi;             // phi_kq
    
    // Calculate edges of bins
    for (i=0; i<=nsteps1; i++) {
      binsvar1[i] = varmin1 + i*(varmax1-varmin1)/nsteps1;
    }
    // binsvar2 filled in loop below
    for (k=0; k<=nsteps3; k++) {
      if (thetaIsUsed) {
        temp = varmax3 - k*(varmax3-varmin3)/nsteps3;   // theta
        binsvar3[k] = cos(temp);                        // cos(theta)
      } else {
        binsvar3[k] = varmin3 + k*(varmax3-varmin3)/nsteps3;
      }
      if (verbose && !itenu && k)
        std::cout << "Bin #" << k << " has cos(theta) = " << 0.5*(binsvar3[k]+binsvar3[k-1])
                  << ",\t width = " << binsvar3[k]-binsvar3[k-1] << std::endl;
    }
    for (l=0; l<=nsteps4; l++) {
      binsvar4[l] = varmin4 + l*(varmax4-varmin4)/nsteps4;
    }
    
    // Calculate central bin values (average of left and right edge of bin)
    for (i=0; i<nsteps1; i++) varvals1[i] = 0.5*(binsvar1[i+1]+binsvar1[i]);
    // binsvar2 filled in loop below
    for (k=0; k<nsteps3; k++) varvals3[k] = 0.5*(binsvar3[k+1]+binsvar3[k]);
    for (l=0; l<nsteps4; l++) varvals4[l] = 0.5*(binsvar4[l+1]+binsvar4[l]);
    
    // CALCULATE CROSS-SECTION
    // -----------------------
    sigma = 0.;  // reset integral over Tkaon
    
    // **********************
    //  Integrate over Tkaon
    // **********************
    for (i=0; i<nsteps1; i++) {
      if (binsvar1[i+1]==binsvar1[i]) continue;
      Tkaon = varvals1[i];
      Ekaon = Tkaon+Mkaon;
      
      // ----------------------------------------------
      // Define binning for second integration variable
      // ----------------------------------------------
      if (pLeptonIsUsed) {                        // plep
        varmin2 = 0.;
        varmax2 = sqrt((Enu-Ekaon)*(Enu-Ekaon)-Mlep*Mlep);
      } else {                                    // Tlep
        varmin2 = 0.;
        varmax2 = Enu-Ekaon-Mlep;
      }
      for (j=0; j<=nsteps2; j++) {
        temp = varmin2 + j*(varmax2-varmin2)/nsteps2;   // plep OR Tlep
        if (pLeptonIsUsed) {
          binsvar2[j] = sqrt(temp*temp+Mlep*Mlep)-Mlep; // plep -> Tlep
        } else {
          binsvar2[j] = temp;                           // Tlep
        }
        if (verbose && !pLeptonIsUsed && !itenu && j)
          std::cout << "Bin edge #" << j << " has Tlep = " << binsvar2[j] 
                    << ",\t width = " << binsvar2[j]-binsvar2[j-1] << std::endl;
      }
      for (j=0; j<nsteps2; j++) varvals2[j] = 0.5*(binsvar2[j+1]+binsvar2[j]);
      // ----------------------------------------------
      
      diff1sigma = 0.;  // reset integral over plep
      
      // *******************************
      //  Integrate over plep (or Tlep)
      // *******************************
      for (j=0; j<nsteps2; j++) {
        if (binsvar2[j+1]==binsvar2[j]) continue;
        Tlep = varvals2[j];
        Elep = Tlep+Mlep;
        
        diff2sigma = 0.;  // reset integral over theta
        
        // ************************************
        //  Integrate over theta (or costheta)
        // ************************************
        for (k=0; k<nsteps3; k++) {
          if (binsvar3[k+1]==binsvar3[k]) continue;
          theta = acos(varvals3[k]);
          
          diff3sigma = 0.;  // reset integral over phikq
          
          // **********************
          //  Integrate over phikq
          // **********************
          for (l=0; l<nsteps4; l++) {
            if (binsvar4[l+1]==binsvar4[l]) continue;
            phikq = varvals4[l];
            
            // Calculate 4D differential xsec
            diff4sigma = inst->diffxsec(Tlep,Tkaon,theta,phikq);
            diff4sigma *= 2*pi;
            
            // Multiply by bin width
            diff3sigma += diff4sigma*(binsvar4[l+1]-binsvar4[l]);
          }
          
          // This Jacobian no longer needed, since binning is adjusted to cos(theta)
          // diff3sigma *= sin(theta);
          
          diff2sigma += diff3sigma*(binsvar3[k+1]-binsvar3[k]);
        }
        
        // Multiplication with Jacobian for transformation plep->Tlep (it is E/p)
        diff2sigma *= Elep/sqrt(Elep*Elep-Mlep*Mlep);
        
        diff1sigma += diff2sigma*(binsvar2[j+1]-binsvar2[j]);
      }
      
      sigma += diff1sigma*(binsvar1[i+1]-binsvar1[i]);
    }
    
    // Total cross section for this energy
    std::cout << Form("%12.6lf",Enu) << "\t" << sigma << std::endl;
    fprintf(outfile, "%lf\t%e\n",Enu,sigma);
  }
  
  fclose(outfile);
  std::cout << std::endl << "Output written to file: " << fpath << fname << std::endl << std::endl;
}

