//____________________________________________________________________________
/*
 Copyright (C) 2003-2006, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Steve Dytman <dytman+@pitt.edu>, Pittsburgh Univ.
         Hugh Gallagher <gallag@minos.phy.tufts.edu>, Tufts Univ.
         Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>, Rutherford Lab.
         October 02, 2006

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <cassert>
#include <string>
#include <iostream>

#include <TSystem.h>
#include <TNtupleD.h>
#include <TTree.h>

#include "HadronTransport/INukeHadroData.h"
#include "Messenger/Messenger.h"
#include "Numerical/Spline.h"

using std::cout;
using std::endl;

using namespace genie;

//____________________________________________________________________________
INukeHadroData * INukeHadroData::fInstance = 0;
//____________________________________________________________________________
INukeHadroData::INukeHadroData()
{
  this->LoadData();
  this->CalcData();

  fInstance = 0;
}
//____________________________________________________________________________
INukeHadroData::~INukeHadroData()
{
  cout << "INukeHadroData singleton dtor: "
                      << "Deleting all hadron cross section splines" << endl;

  //-- delete SAID h+N x-section splines
  delete fXSecPipP_Elas;
  delete fXSecPipP_Reac;
  delete fXSecPipD_Abs;    
  delete fXSecPimP_Elas;
  delete fXSecPimP_Reac;
  delete fXSecPimP_CEx; 
  delete fXSecPP_Elas;      
  delete fXSecPP_Reac;      
  delete fXSecNP_Elas;      
  delete fXSecNP_Reac;      

  //-- delete Mashnik p+Fe x-section splines
  delete fXSecPFe_Elas;
  delete fXSecPFe_Reac;
  delete fXSecPFe_P;   
  delete fXSecPFe_N;   
  delete fXSecPFe_NP;  
  delete fXSecPFe_PP;  
  delete fXSecPFe_NPP; 
  delete fXSecPFe_NNP; 
  delete fXSecPFe_NNPP;
  delete fXSecPFe_Pim; 
  delete fXSecPFe_Pi0; 
  delete fXSecPFe_Pip; 

  //-- delete Mashnik pi+Fe x-section splines
  delete fXSecPiFe_P;    
  delete fXSecPiFe_PP;   
  delete fXSecPiFe_PPP;  
  delete fXSecPiFe_N;    
  delete fXSecPiFe_NN;   
  delete fXSecPiFe_NNN;  
  delete fXSecPiFe_NP;   
  delete fXSecPiFe_NPP;  
  delete fXSecPiFe_NPPP; 
  delete fXSecPiFe_NNP;  
  delete fXSecPiFe_NNPP; 
  delete fXSecPiFe_Pi0;  

  //-- delete x-sections from data (Ash/Carrol)
  delete fXSecAshPiFe_Abs; 
  delete fXSecAshPiFe_Reac;
  delete fXSecCarPiFe_Tot; 

  //-- delete total x-sections 
  delete fXSecPip_Tot;   
  delete fXSecPim_Tot;   
  delete fXSecPi0_Tot;   
  delete fXSecP_Tot;     
  delete fXSecN_Tot;     

  //-- delete x-section fractions
  delete fFracPip_CEx;   
  delete fFracPip_Elas;  
  delete fFracPip_Reac;  
  delete fFracPip_Abs;   
  delete fFracPim_CEx;   
  delete fFracPim_Elas;  
  delete fFracPim_Reac;  
  delete fFracPim_Abs;   
  delete fFracPi0_CEx;    
  delete fFracPi0_Elas;   
  delete fFracPi0_Reac;   
  delete fFracPi0_Abs;    
  delete fFracP_Reac;     
  delete fFracN_Reac;     
  delete fFracPiA_Elas;   
  delete fFracPiA_Inel;   
  delete fFracPiA_CEx;    
  delete fFracPiA_Abs;    
  delete fFracPiA_PP;     
  delete fFracPiA_NPP;    
  delete fFracPiA_NNP;    
  delete fFracPiA_4N4P;   
  delete fFracPiA_PiProd; 
  delete fFracPA_Elas;    
  delete fFracPA_Inel;    
  delete fFracPA_Abs;     
  delete fFracPA_PP;      
  delete fFracPA_NPP;     
  delete fFracPA_NNP;     
  delete fFracPA_4N4P;    
  delete fFracPA_PiProd;  
}
//____________________________________________________________________________
INukeHadroData * INukeHadroData::Instance()
{
  if(fInstance == 0) {
    LOG("INukeData", pINFO) << "INukeHadroData late initialization";
    static INukeHadroData::Cleaner cleaner;
    cleaner.DummyMethodAndSilentCompiler();
    fInstance = new INukeHadroData;
  }
  return fInstance;
}
//____________________________________________________________________________
void INukeHadroData::LoadData(void)
{
// Loads hadronic x-section data

  //-- Get the directory with the SAID hadron cross section data (search for
  //   $GINUKEHADRONDATA or use default: $GENIE/data/hadron_xsec)
  string data_dir = (gSystem->Getenv("GINUKEHADRONDATA")) ?
             string(gSystem->Getenv("GINUKEHADRONDATA")) :
             string(gSystem->Getenv("GENIE")) + string("/data/hadron_xsec");

  LOG("INukeData", pNOTICE)  
      << "Loading INTRANUKE hadron data from: " << data_dir;

  //-- Build filenames

  string fnm_said_pip_p_elas = data_dir + "/pi+p-xselas-said.dat";
  string fnm_said_pip_p_reac = data_dir + "/pi+p-xsreac-said.dat";
  string fnm_said_pim_p_elas = data_dir + "/pi-p-xselas-said.dat";
  string fnm_said_pim_p_reac = data_dir + "/pi-p-xsreac-said.dat";
  string fnm_said_pim_p_cex  = data_dir + "/pi-p-cex-xstot-said.dat";
  string fnm_said_pid_pp_tot = data_dir + "/pid-pp-xstot-said.dat";
  string fnm_said_np_elas    = data_dir + "/np-xselas-said.dat";
  string fnm_said_np_reac    = data_dir + "/np-xsreac-said.dat";
  string fnm_said_pp_elas    = data_dir + "/pp-xselas-said.dat";
  string fnm_said_pp_reac    = data_dir + "/pp-xsreac-said.dat";
  string fnm_mashnik_pFe     = data_dir + "/mashnik-pfe-xs.dat";
  string fnm_mashnik_piFe    = data_dir + "/mashnik-pife-xs.dat";
  string fnm_ash_piFe_abs    = data_dir + "/exp-ash-pife-xsabs.dat";
  string fnm_ash_piFe_reac   = data_dir + "/exp-ash-pife-xsreac.dat";
  string fnm_carrol_piFe_tot = data_dir + "/exp-carrol-pife-xstot.dat";

  //-- Make sure that all data files are available

  assert( ! gSystem->AccessPathName(fnm_said_pip_p_elas.c_str()) );
  assert( ! gSystem->AccessPathName(fnm_said_pip_p_reac.c_str()) );
  assert( ! gSystem->AccessPathName(fnm_said_pim_p_elas.c_str()) );
  assert( ! gSystem->AccessPathName(fnm_said_pim_p_reac.c_str()) );
  assert( ! gSystem->AccessPathName(fnm_said_pim_p_cex.c_str())  );
  assert( ! gSystem->AccessPathName(fnm_said_pid_pp_tot.c_str()) );
  assert( ! gSystem->AccessPathName(fnm_said_np_elas.c_str())    );
  assert( ! gSystem->AccessPathName(fnm_said_np_reac.c_str())    );
  assert( ! gSystem->AccessPathName(fnm_said_pp_elas.c_str())    );
  assert( ! gSystem->AccessPathName(fnm_said_pp_reac.c_str())    );
  assert( ! gSystem->AccessPathName(fnm_mashnik_pFe.c_str())     );
  assert( ! gSystem->AccessPathName(fnm_mashnik_piFe.c_str())    );
  assert( ! gSystem->AccessPathName(fnm_ash_piFe_abs.c_str())    );
  assert( ! gSystem->AccessPathName(fnm_ash_piFe_reac.c_str())   );
  assert( ! gSystem->AccessPathName(fnm_carrol_piFe_tot.c_str()) );

  LOG("INukeData", pNOTICE)  << "Found all necessary data files...";

  //-- Load the SAID data into splines
  
  fXSecPipP_Elas = new Spline ( fnm_said_pip_p_elas );
  fXSecPipP_Reac = new Spline ( fnm_said_pip_p_reac );
  fXSecPimP_Elas = new Spline ( fnm_said_pim_p_elas );
  fXSecPimP_Reac = new Spline ( fnm_said_pim_p_reac );
  fXSecPimP_CEx  = new Spline ( fnm_said_pim_p_cex  );
  fXSecNP_Elas   = new Spline ( fnm_said_np_elas    );
  fXSecNP_Reac   = new Spline ( fnm_said_np_reac    );
  fXSecPP_Elas   = new Spline ( fnm_said_pp_elas    );   
  fXSecPP_Reac   = new Spline ( fnm_said_pp_reac    );   

  // Especially for pion absorption cross section copy S.Dytmnan's correction.
  // His comments: 'model not valid above 600 MeV, slow linear falloff matches data'

  Spline xsabs(fnm_said_pid_pp_tot); 

  double   frac = 1; //?
  int      nk   = xsabs.NKnots();
  double * xsec = new double[nk];
  double * E    = new double[nk];

  for(int i=0; i<nk; i++) {
    xsec[i] = (E[i] < 600.) ? frac * xsabs.Evaluate(E[i]) : 
                              frac * xsabs.Evaluate(600.)*(1800.-E[i])/1200.;
  }
  fXSecPipD_Abs = new Spline(nk,E,xsec); // corrected pi abs xsec

  delete xsec;
  delete E;

  LOG("INukeData", pNOTICE) 
      << "... Done loading SAID hadron cross section data";

  //-- Load the Mashnik data into splines

  TTree tpFe;
  tpFe.ReadFile(fnm_mashnik_pFe.c_str(),
     "ke/D:p/D:np/D:pp/D:npp/D:pim/D:nnpp/D:pi0/D:nnp/D:pip/D:n/D:reac/D:elas/D");

  fXSecPFe_Elas  = new Spline(&tpFe,"ke:elas"); 
  fXSecPFe_Reac  = new Spline(&tpFe,"ke:reac"); 
  fXSecPFe_P     = new Spline(&tpFe,"ke:p"); 
  fXSecPFe_N     = new Spline(&tpFe,"ke:n"); 
  fXSecPFe_NP    = new Spline(&tpFe,"ke:np"); 
  fXSecPFe_PP    = new Spline(&tpFe,"ke:pp"); 
  fXSecPFe_NPP   = new Spline(&tpFe,"ke:npp"); 
  fXSecPFe_NNP   = new Spline(&tpFe,"ke:nnp"); 
  fXSecPFe_NNPP  = new Spline(&tpFe,"ke:nnpp"); 
  fXSecPFe_Pim   = new Spline(&tpFe,"ke:pim"); 
  fXSecPFe_Pi0   = new Spline(&tpFe,"ke:pi0"); 
  fXSecPFe_Pip   = new Spline(&tpFe,"ke:pip"); 

  TTree tpiFe;
  tpiFe.ReadFile(fnm_mashnik_piFe.c_str(),
      "ke/D:p/D:pp/D:ppp/D:n/D:nn/D:nnn/D:np/D:npp/D:nppp/D:nnp/D:nnpp/D:pi0/D");

  fXSecPiFe_P    = new Spline(&tpiFe,"ke:p"); 
  fXSecPiFe_PP   = new Spline(&tpiFe,"ke:pp"); 
  fXSecPiFe_PPP  = new Spline(&tpiFe,"ke:ppp"); 
  fXSecPiFe_N    = new Spline(&tpiFe,"ke:n"); 
  fXSecPiFe_NN   = new Spline(&tpiFe,"ke:nn"); 
  fXSecPiFe_NNN  = new Spline(&tpiFe,"ke:nnn"); 
  fXSecPiFe_NP   = new Spline(&tpiFe,"ke:np"); 
  fXSecPiFe_NPP  = new Spline(&tpiFe,"ke:npp"); 
  fXSecPiFe_NPPP = new Spline(&tpiFe,"ke:nppp"); 
  fXSecPiFe_NNP  = new Spline(&tpiFe,"ke:nnp"); 
  fXSecPiFe_NNPP = new Spline(&tpiFe,"ke:nnpp"); 
  fXSecPiFe_Pi0  = new Spline(&tpiFe,"ke:pi0"); 

  LOG("INukeData", pNOTICE) 
      << "... Done loading Mashnik hadron cross section data";

  fXSecAshPiFe_Abs  = new Spline( fnm_ash_piFe_abs    );
  fXSecAshPiFe_Reac = new Spline( fnm_ash_piFe_reac   );
  fXSecCarPiFe_Tot  = new Spline( fnm_carrol_piFe_tot );

  LOG("INukeData", pNOTICE) 
      << "... Done loading Ash & Carrol hadron cross section data";
};
//____________________________________________________________________________
void INukeHadroData::CalcData(void)
{
// Hopefully this esction would be simplified (be made redundant) by making all 
// the corrections offline and loading the corrected x-sections
//
  LOG("INukeData", pNOTICE) 
      << "Computing missing x-sections & applying corrections...";

  //
  //-- corrections to Mashnik total pi0 and Carrol total x-section
  //

  // Mashnik: divide pi0 to match Bowles data
  fXSecPiFe_Pi0->Divide(2.0); 

  // This extrapolation is empirical, based on total x-sections at similar
  // energies measured by Cluogh et al for lighter targets
/*
  double fac0  = 44.8;
  double fac1  = -0.045;
  double fac2  =  0.000048;
  double expo1 =  0.9715;
  double expo2 = -0.000215;

  loop from bins 16 to 29... and do something...
  Stuff like this make make the code more error prone and obscure that it needs
  to be.
*/

  //
  // -- hN mode: Build total x-sections and fractions splines
  //

  //total x-sections
  TNtupleD * ntxs = new TNtupleD("tot","total x-sec","ke:pip:pim:pi0:p:n");

  //pi0 fractions
  TNtupleD * ntfrpi0 = new TNtupleD("ntfrpi0","","ke:cex:elas:reac:abs");

  int nknots = fXSecPipP_Elas->NKnots();
  for(int i=0; i<nknots; i++) {
    double ke        = fXSecPipP_Elas->GetKnotX(i);  // kinetic energy
    double xspipelas = fXSecPipP_Elas->Evaluate(ke);
    double xspipreac = fXSecPipP_Reac->Evaluate(ke);
    double xspiabs   = fXSecPipD_Abs ->Evaluate(ke);

    LOG("INukeData", pDEBUG) 
      << "pi+ @ ke = " << ke << ": elas =  " <<  xspipelas 
      << ", reac = " << xspipreac << ", abs = " <<  xspiabs;

    double xspimelas = fXSecPimP_Elas->Evaluate(ke);
    double xspimreac = fXSecPimP_Reac->Evaluate(ke);
    double xspimcex  = fXSecPimP_CEx ->Evaluate(ke);

    LOG("INukeData", pDEBUG) 
      << "pi- @ at ke = " << ke << ": elas =  " <<  xspimelas 
      << ", reac = " << xspimreac << ", cex = " <<  xspimcex;

    double xspi0pelas = (5./18.)*xspipelas - (1./6.)*xspimelas + (5./6.)*xspimcex;
    double xspi0nelas = (1./2. )*xspipelas + (1./2.)*xspimelas - (1./2.)*xspimcex;
    double xspi0reac  = (xspimreac+xspipreac)/2.;
    double xspi0cex   =  xspimcex;
    double xspi0abs   =  xspiabs;

    double xsppelas  = fXSecPP_Elas->Evaluate(ke);
    double xsnpelas  = fXSecNP_Elas->Evaluate(ke);
    double xsppreac  = fXSecPP_Reac->Evaluate(ke);
    double xsnpreac  = fXSecNP_Reac->Evaluate(ke);

    LOG("INukeData", pDEBUG) 
      << "pp @ ke = " << ke << ": elas =  " <<  xsppelas 
      << ", reac = " << xsppreac;
    LOG("INukeData", pDEBUG) 
      << "np @ ke = " << ke << ": elas =  " <<  xsnpelas 
      << ", reac = " << xsnpreac;

    double xsplptot  = xspipelas + xspipreac +  xspiabs;
    double xsplntot  = xspimelas + xspimreac +  xspiabs;
    double xspi0ptot = xspi0pelas + xspi0cex + xspi0reac + xspi0abs;
    double xspi0ntot = xspi0nelas + xspi0cex + xspi0reac + xspi0abs;

    double xspip = (xsplptot  + xsplntot)/2.;
    double xspim = (xsplptot  + xsplntot)/2.;
    double xspi0 = (xspi0ptot + xspi0ntot)/2.;
    double xsp   = (xsppelas  + xsnpelas + xsppreac + xsnpreac)/2.;
    double xsn   =  xsnpelas  + xsnpreac;

    ntxs->Fill(ke,xspip,xspim,xspi0,xsp,xsn);

    //-- pi0 fractions
    if(xspi0>0) {
      double frpi0cex  = xspi0cex/xspi0;
      double frpi0elas = 1 - frpi0cex;
      double frpi0reac = frpi0elas - 0.5*(xspi0pelas+xspi0nelas)/xspi0;
      double frpi0abs  = frpi0reac - xspi0reac/xspi0;

      ntfrpi0->Fill(ke,frpi0cex,frpi0elas,frpi0reac,frpi0abs);
    } 
  }

  fXSecPip_Tot    = new Spline(ntxs,"ke:pip");
  fXSecPim_Tot    = new Spline(ntxs,"ke:pim");
  fXSecPi0_Tot    = new Spline(ntxs,"ke:pi0");
  fXSecP_Tot      = new Spline(ntxs,"ke:p");
  fXSecN_Tot      = new Spline(ntxs,"ke:n");

  fFracPi0_CEx    = new Spline(ntfrpi0,"ke:cex");
  fFracPi0_Elas   = new Spline(ntfrpi0,"ke:elas");
  fFracPi0_Reac   = new Spline(ntfrpi0,"ke:reac");
  fFracPi0_Abs    = new Spline(ntfrpi0,"ke:abs");

  delete ntxs;
  delete ntfrpi0;

  // the following 10 splines are not filled-in yet...
  fFracPip_CEx    = new Spline; 
  fFracPip_Elas   = new Spline; 
  fFracPip_Reac   = new Spline; 
  fFracPip_Abs    = new Spline; 
  fFracPim_CEx    = new Spline; 
  fFracPim_Elas   = new Spline; 
  fFracPim_Reac   = new Spline; 
  fFracPim_Abs    = new Spline; 
  fFracP_Reac     = new Spline; 
  fFracN_Reac     = new Spline; 

  //
  //-- hA mode: calculate p+A fractions
  //
  LOG("INukeData", pNOTICE) << "Calculating p+A fractions...";

  TNtupleD ntfrpA("ntfrpA","","ke:elas:inel:abs:pp:npp:nnp:n4p4:piprod"); 
  nknots = fXSecPFe_Elas->NKnots();
  for(int i=0; i<nknots; i++) {
      double ke = fXSecPFe_Elas->GetKnotX(i); // kinetic energy
      ke = TMath::Max(1.,ke);
      double sigreac        = fXSecPFe_Reac->Evaluate(ke);
      double sigelas        = fXSecPFe_Elas->Evaluate(ke);
      double sigtot         = sigreac + sigelas;
      // if(i==0) sigtot=1; // Steve's comment: to avoid NaN !
      if(sigtot==0) sigtot=1;
      double  siginel       = fXSecPFe_P->Evaluate(ke);
      double sigcex         = fXSecPFe_N->Evaluate(ke);
      double sig_pfe_np     = fXSecPFe_NP->Evaluate(ke);
      double sig_pfe_pp     = fXSecPFe_PP->Evaluate(ke);
      double sig_pfe_npp    = fXSecPFe_NPP->Evaluate(ke);
      double sig_pfe_nnp    = fXSecPFe_NNP->Evaluate(ke);
      double sig_pfe_pip    = fXSecPFe_Pip->Evaluate(ke);

      double frac_pa_elas   = TMath::Max(0., 1 - sigcex/sigtot);
      double frac_pa_inel   = TMath::Max(0., frac_pa_elas   - sigelas / sigtot);
      double frac_pa_abs    = TMath::Max(0., frac_pa_inel   - siginel / sigtot);
      double frac_pa_pp     = TMath::Max(0., frac_pa_abs    - sig_pfe_np / sigtot);
      double frac_pa_npp    = TMath::Max(0., frac_pa_pp     - sig_pfe_pp / sigtot);
      double frac_pa_nnp    = TMath::Max(0., frac_pa_npp    - sig_pfe_npp / sigtot);
      double frac_pa_4n4p   = TMath::Max(0., frac_pa_nnp    - sig_pfe_nnp / sigtot);
      double frac_pa_piprod = TMath::Max(0., sig_pfe_pip / sigtot);

      SLOG("INukeData", pDEBUG) 
        << "\n ******* pA @ ke = " << ke << "\n"
        << "- x-sections: \n"
        << "reac = "  << sigreac     << ", elas = " << sigelas 
        << ", tot = " << sigtot      << ", inel = " << siginel  
        << ", cex = " << sigcex      << ", pp = "   << sig_pfe_pp   
        << ", np = "  << sig_pfe_np  << ", nnp = "  << sig_pfe_nnp 
        << ", npp = " << sig_pfe_npp << ", piprod = " << sig_pfe_pip << "\n"
        << "- fractions: \n"
        << "elas = "  << frac_pa_elas  << ". inel = "   << frac_pa_inel     
        << ", abs = " << frac_pa_abs   << ", pp = "     << frac_pa_pp
        << ", npp = "  << frac_pa_npp  << ", nnp = "    << frac_pa_nnp  
        << ", 4n4p = " << frac_pa_4n4p << ", piprod = " << frac_pa_piprod;

      ntfrpA.Fill(ke, frac_pa_elas, frac_pa_inel, frac_pa_abs, frac_pa_pp, 
                  frac_pa_npp, frac_pa_nnp, frac_pa_4n4p, frac_pa_piprod);
  }
  fFracPA_Elas    = new Spline(&ntfrpA, "ke:elas"  );
  fFracPA_Inel    = new Spline(&ntfrpA, "ke:inel"  );
  fFracPA_Abs     = new Spline(&ntfrpA, "ke:abs"   );
  fFracPA_PP      = new Spline(&ntfrpA, "ke:pp"    );
  fFracPA_NPP     = new Spline(&ntfrpA, "ke:npp"   );
  fFracPA_NNP     = new Spline(&ntfrpA, "ke:nnp"   );
  fFracPA_4N4P    = new Spline(&ntfrpA, "ke:n4p4"  );
  fFracPA_PiProd  = new Spline(&ntfrpA, "ke:piprod");

  //
  //-- hA mode: calculate pi+A fractions
  //
  LOG("INukeData", pINFO) << "Calculating pi+A fractions...";

  TNtupleD ntfrpiA("ntfrpiA","","ke:elas:inel:cex:abs:pp:npp:nnp:n4p4:piprod"); 
  nknots = fXSecPiFe_Pi0->NKnots();
  for(int i=0; i<nknots; i++) {
      double ke = fXSecPiFe_NP->GetKnotX(i); // kinetic energy
      ke = TMath::Max(1.,ke);
      double sig_pife_pi0    = fXSecPiFe_Pi0->Evaluate(ke);
      // car x-section goes only up to 320 MeV - take it to be constant above that
      double sigtot          = (ke<320) ? fXSecCarPiFe_Tot->Evaluate(ke) : 
                                         fXSecCarPiFe_Tot->Evaluate(310);
      double sigabs          = fXSecAshPiFe_Abs->Evaluate(ke); 
      double sigreac         = fXSecAshPiFe_Reac->Evaluate(ke);
      double sigelas         = TMath::Max(0., sigtot - sigreac);
      double siginel         = TMath::Max(0., sigreac - sig_pife_pi0 - sigabs);
      // this is an estimate of the single pi0 x-section from Ashery (elas only)
      double sigcex          = 2.2 * fXSecPimP_CEx->Evaluate(ke);
      double sig_cex_elas;
      double sig_cex_inel;
      if(ke<400) {
          sig_cex_elas = sig_pife_pi0;
          sig_cex_inel = 0.;
      } else {
          sig_cex_elas = sigcex;
          sig_cex_inel = TMath::Max(0., sig_pife_pi0 - sigcex);
      }
      double sig_pife_np     = fXSecPiFe_NP->Evaluate(ke);
      double sig_pife_pp     = fXSecPiFe_PP->Evaluate(ke);
      double sig_pife_npp    = fXSecPiFe_NPP->Evaluate(ke);
      double sig_pife_nnp    = fXSecPiFe_NNP->Evaluate(ke);

      double frac_pia_elas   = TMath::Max(0., 1 - sig_cex_elas / sigtot);
      double frac_pia_inel   = TMath::Max(0., frac_pia_elas - sigelas / sigtot);
      double frac_pia_abs    = TMath::Max(0., frac_pia_inel - siginel / sigtot);
      double frac_pia_pp     = TMath::Max(0., frac_pia_abs  - sig_pife_np / sigtot);
      double frac_pia_npp    = TMath::Max(0., frac_pia_pp   - sig_pife_pp / sigtot);
      double frac_pia_nnp    = TMath::Max(0., frac_pia_npp  - sig_pife_npp / sigtot);
      double frac_pia_4n4p   = TMath::Max(0., frac_pia_nnp  - sig_pife_nnp / sigtot);
      double frac_pia_piprod = TMath::Max(0., sig_cex_inel / sigtot);
      double frac_pia_cex    = 0; // ???

      SLOG("INukeData", pDEBUG) 
         << "\n ******* piA @ ke = " << ke << "\n"
         << "- x-sections: \n"
         << "elas = " << sigelas << ", inel = " << siginel << ", tot = " << sigtot
         << ", abs = " << sigabs << ", reac = " << sigreac << ", cex = " << sigcex
         << ", cex/elas = " << sig_cex_elas << ", cex/inel = " << sig_cex_inel
         << ", pi0 = " << sig_pife_pi0 << ", np = " << sig_pife_np << ", pp = " << sig_pife_pp
         << ", npp = " << sig_pife_npp << ", nnp = " << sig_pife_nnp << "\n"
         << "- fractions: \n"
         << "elas = "  << frac_pia_elas  << ". inel = " << frac_pia_inel     
         << ", abs = " << frac_pia_abs   << ", pp = "   << frac_pia_pp
         << ", npp = "  << frac_pia_npp  << ", nnp = "  << frac_pia_nnp  
         << ", 4n4p = " << frac_pia_4n4p << ", piprod = " << frac_pia_piprod;

      ntfrpiA.Fill(ke, frac_pia_elas, frac_pia_inel, frac_pia_cex,
                   frac_pia_abs, frac_pia_pp, frac_pia_npp, frac_pia_nnp,
                   frac_pia_4n4p, frac_pia_piprod); 
  }
  fFracPiA_Elas   = new Spline(&ntfrpiA, "ke:elas"); 
  fFracPiA_Inel   = new Spline(&ntfrpiA, "ke:inel"); 
  fFracPiA_CEx    = new Spline(&ntfrpiA, "ke:cex"); 
  fFracPiA_Abs    = new Spline(&ntfrpiA, "ke:abs"); 
  fFracPiA_PP     = new Spline(&ntfrpiA, "ke:pp"); 
  fFracPiA_NPP    = new Spline(&ntfrpiA, "ke:npp"); 
  fFracPiA_NNP    = new Spline(&ntfrpiA, "ke:nnp"); 
  fFracPiA_4N4P   = new Spline(&ntfrpiA, "ke:n4p4"); 
  fFracPiA_PiProd = new Spline(&ntfrpiA, "ke:piprod"); 

  LOG("INukeData", pNOTICE) 
       << "... Done computing total hadron cross section data";
}
//____________________________________________________________________________
