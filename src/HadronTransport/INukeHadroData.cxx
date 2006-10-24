//____________________________________________________________________________
/*
 Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
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
  delete fXsPipPElas;
  delete fXsPipPReac;
  delete fXsPipDAbs;    
  delete fXsPimPElas;
  delete fXsPimPReac;
  delete fXsPimPCEx; 
  delete fXsPPElas;      
  delete fXsPPReac;      
  delete fXsNPElas;      
  delete fXsNPReac;      

  //-- delete Mashnik p+Fe x-section splines
  delete fXsPFeElas;
  delete fXsPFeReac;
  delete fXsPFeP;   
  delete fXsPFePP;  
  delete fXsPFeNPP; 
  delete fXsPFeNNP; 
  delete fXsPFeNNPP;
  delete fXsPFePim; 
  delete fXsPFePi0; 
  delete fXsPFePip; 

  //-- delete Mashnik pi+Fe x-section splines
  delete fXsPiFeP;    
  delete fXsPiFePP;   
  delete fXsPiFePPP;  
  delete fXsPiFeN;    
  delete fXsPiFeNN;   
  delete fXsPiFeNNN;  
  delete fXsPiFeNP;   
  delete fXsPiFeNPP;  
  delete fXsPiFeNPPP; 
  delete fXsPiFeNNP;  
  delete fXsPiFeNNPP; 
  delete fXsPiFePi0;  

  //-- delete x-sections from data (Ash/Carrol)
  delete fXsAshPiFeAbs; 
  delete fXsAshPiFeReac;
  delete fXsCarPiFeTot; 

  //-- delete total x-sections 
  delete fXsPipTot;   
  delete fXsPimTot;   
  delete fXsPi0Tot;   
  delete fXsPTot;     
  delete fXsNTot;     

  //-- delete x-section fractions
  delete fFrPipCEx;   
  delete fFrPipElas;  
  delete fFrPipReac;  
  delete fFrPipAbs;   
  delete fFrPimCEx;   
  delete fFrPimElas;  
  delete fFrPimReac;  
  delete fFrPimAbs;   
  delete fFrPi0CEx;    
  delete fFrPi0Elas;   
  delete fFrPi0Reac;   
  delete fFrPi0Abs;    
  delete fFrPReac;     
  delete fFrNReac;     
  delete fFrPiAElas;   
  delete fFrPiAInel;   
  delete fFrPiACEx;    
  delete fFrPiAAbs;    
  delete fFrPiAPP;     
  delete fFrPiANPP;    
  delete fFrPiANNP;    
  delete fFrPiA4N4P;   
  delete fFrPiAPiProd; 
  delete fFrPAElas;    
  delete fFrPAInel;    
  delete fFrPAAbs;     
  delete fFrPAPP;      
  delete fFrPANPP;     
  delete fFrPANNP;     
  delete fFrPA4N4P;    
  delete fFrPAPiProd;  
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
  
  fXsPipPElas = new Spline ( fnm_said_pip_p_elas );
  fXsPipPReac = new Spline ( fnm_said_pip_p_reac );
  fXsPimPElas = new Spline ( fnm_said_pim_p_elas );
  fXsPimPReac = new Spline ( fnm_said_pim_p_reac );
  fXsPimPCEx  = new Spline ( fnm_said_pim_p_cex  );
  fXsNPElas   = new Spline ( fnm_said_np_elas    );
  fXsNPReac   = new Spline ( fnm_said_np_reac    );
  fXsPPElas   = new Spline ( fnm_said_pp_elas    );   
  fXsPPReac   = new Spline ( fnm_said_pp_reac    );   

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
  fXsPipDAbs = new Spline(nk,E,xsec); // corrected pi abs xsec

  delete [] xsec;
  delete [] E;

  LOG("INukeData", pNOTICE) 
      << "... Done loading SAID hadron cross section data";

  //-- Load the Mashnik data into splines

  TTree tpFe;
  tpFe.ReadFile(fnm_mashnik_pFe.c_str(),
     "ke/D:p/D:np/D:pp/D:npp/D:pim/D:nnpp/D:pi0/D:nnp/D:pip/D:n/D:reac/D:elas/D");

  fXsPFeElas  = new Spline(&tpFe,"ke:elas"); 
  fXsPFeReac  = new Spline(&tpFe,"ke:reac"); 
  fXsPFeP     = new Spline(&tpFe,"ke:p"); 
  fXsPFePP    = new Spline(&tpFe,"ke:pp"); 
  fXsPFeNPP   = new Spline(&tpFe,"ke:npp"); 
  fXsPFeNNP   = new Spline(&tpFe,"ke:nnp"); 
  fXsPFeNNPP  = new Spline(&tpFe,"ke:nnpp"); 
  fXsPFePim   = new Spline(&tpFe,"ke:pim"); 
  fXsPFePi0   = new Spline(&tpFe,"ke:pi0"); 
  fXsPFePip   = new Spline(&tpFe,"ke:pip"); 

  TTree tpiFe;
  tpiFe.ReadFile(fnm_mashnik_piFe.c_str(),
      "ke/D:p/D:pp/D:ppp/D:n/D:nn/D:nnn/D:np/D:npp/D:nppp/D:nnp/D:nnpp/D:pi0/D");

  fXsPiFeP    = new Spline(&tpiFe,"ke:p"); 
  fXsPiFePP   = new Spline(&tpiFe,"ke:pp"); 
  fXsPiFePPP  = new Spline(&tpiFe,"ke:ppp"); 
  fXsPiFeN    = new Spline(&tpiFe,"ke:n"); 
  fXsPiFeNN   = new Spline(&tpiFe,"ke:nn"); 
  fXsPiFeNNN  = new Spline(&tpiFe,"ke:nnn"); 
  fXsPiFeNP   = new Spline(&tpiFe,"ke:np"); 
  fXsPiFeNPP  = new Spline(&tpiFe,"ke:npp"); 
  fXsPiFeNPPP = new Spline(&tpiFe,"ke:nppp"); 
  fXsPiFeNNP  = new Spline(&tpiFe,"ke:nnp"); 
  fXsPiFeNNPP = new Spline(&tpiFe,"ke:nnpp"); 
  fXsPiFePi0  = new Spline(&tpiFe,"ke:pi0"); 

  LOG("INukeData", pNOTICE) 
      << "... Done loading Mashnik hadron cross section data";

  fXsAshPiFeAbs  = new Spline( fnm_ash_piFe_abs    );
  fXsAshPiFeReac = new Spline( fnm_ash_piFe_reac   );
  fXsCarPiFeTot  = new Spline( fnm_carrol_piFe_tot );

  LOG("INukeData", pNOTICE) 
      << "... Done loading Ash & Carrol hadron cross section data";
};
//____________________________________________________________________________
void INukeHadroData::CalcData(void)
{
// Calculates hadronic x-sections 

  LOG("INukeData", pNOTICE) 
      << "Computing missing x-sections & applying corrections...";

  //total x-sections
  TNtupleD * ntxs = new TNtupleD("tot","total x-sec","ke:pip:pim:pi0:p:n");

  //pi0 fractions
  TNtupleD * ntfrpi0 = new TNtupleD("ntfrpi0","","ke:cex:elas:reac:abs");

  int nknots = fXsPipPElas->NKnots();

  for(int i=0; i<nknots; i++) {

    double ke, tmp;
    fXsPipPElas->GetAsTSpline()->GetKnot(i,ke,tmp);

    double xspipelas = fXsPipPElas->Evaluate(ke);
    double xspipreac = fXsPipPReac->Evaluate(ke);
    double xspiabs   = fXsPipDAbs ->Evaluate(ke);

    LOG("INukeData", pDEBUG) 
      << "pi+ @ ke = " << ke << ": elas =  " <<  xspipelas 
      << ", reac = " << xspipreac << ", abs = " <<  xspiabs;

    double xspimelas = fXsPimPElas->Evaluate(ke);
    double xspimreac = fXsPimPReac->Evaluate(ke);
    double xspimcex  = fXsPimPCEx ->Evaluate(ke);

    LOG("INukeData", pDEBUG) 
      << "pi- @ at ke = " << ke << ": elas =  " <<  xspimelas 
      << ", reac = " << xspimreac << ", cex = " <<  xspimcex;

    double xspi0pelas = (5./18.)*xspipelas - (1./6.)*xspimelas + (5./6.)*xspimcex;
    double xspi0nelas = (1./2. )*xspipelas + (1./2.)*xspimelas - (1./2.)*xspimcex;
    double xspi0reac  = (xspimreac+xspipreac)/2.;
    double xspi0cex   =  xspimcex;
    double xspi0abs   =  xspiabs;

    double xsppelas  = fXsPPElas->Evaluate(ke);
    double xsnpelas  = fXsNPElas->Evaluate(ke);
    double xsppreac  = fXsPPReac->Evaluate(ke);
    double xsnpreac  = fXsNPReac->Evaluate(ke);

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

  fXsPipTot    = new Spline(ntxs,"ke:pip");
  fXsPimTot    = new Spline(ntxs,"ke:pim");
  fXsPi0Tot    = new Spline(ntxs,"ke:pi0");
  fXsPTot      = new Spline(ntxs,"ke:p");
  fXsNTot      = new Spline(ntxs,"ke:n");

  fFrPi0CEx    = new Spline(ntfrpi0,"ke:cex");
  fFrPi0Elas   = new Spline(ntfrpi0,"ke:elas");
  fFrPi0Reac   = new Spline(ntfrpi0,"ke:reac");
  fFrPi0Abs    = new Spline(ntfrpi0,"ke:abs");

  delete ntxs;
  delete ntfrpi0;



  fFrPipCEx    = new Spline; 
  fFrPipElas   = new Spline; 
  fFrPipReac   = new Spline; 
  fFrPipAbs    = new Spline; 
  fFrPimCEx    = new Spline; 
  fFrPimElas   = new Spline; 
  fFrPimReac   = new Spline; 
  fFrPimAbs    = new Spline; 

  fFrPReac     = new Spline; 
  fFrNReac     = new Spline; 
  fFrPiAElas   = new Spline; 
  fFrPiAInel   = new Spline; 
  fFrPiACEx    = new Spline; 
  fFrPiAAbs    = new Spline; 
  fFrPiAPP     = new Spline; 
  fFrPiANPP    = new Spline; 
  fFrPiANNP    = new Spline; 
  fFrPiA4N4P   = new Spline; 
  fFrPiAPiProd = new Spline; 
  fFrPAElas    = new Spline; 
  fFrPAInel    = new Spline; 
  fFrPAAbs     = new Spline; 
  fFrPAPP      = new Spline; 
  fFrPANPP     = new Spline; 
  fFrPANNP     = new Spline; 
  fFrPA4N4P    = new Spline; 
  fFrPAPiProd  = new Spline; 



  LOG("INukeData", pNOTICE) 
       << "... Done computing total hadron cross section data";
}
//____________________________________________________________________________
