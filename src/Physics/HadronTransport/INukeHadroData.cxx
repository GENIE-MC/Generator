//____________________________________________________________________________
/*
 Copyright (C) 2003-2015, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 

 Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>, Rutherford Lab.
         Steve Dytman <dytman+@pitt.edu>, Pittsburgh Univ.
	 Aaron Meyer <asm58@pitt.edu>, Pittsburgh Univ.
	 Alex Bell, Pittsburgh Univ.
         February 01, 2007

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Dec 06, 2008 - CA
   Tweak dtor so as not to clutter the output if GENIE exits in err so as to
   spot the fatal mesg immediately.
 @ Jul 15, 2010 - AM
   MeanFreePath now distinguishes between protons and neutrons. To account for
   this, additional splines were added for each nucleon. Absorption fates
   condensed into a single fate, splines updated accordingly. Added gamma and kaon
   splines.Outdated splines were removed. Function IntBounce implemented to calculate
   a CM scattering angle for given probe, target, product, and fate. AngleAndProduct
   similar to IntBounce, but also determines the target nucleon.
 @ May 01, 2012 - CA
   Pick data from $GENIE/data/evgen/intranuke/

*/
//____________________________________________________________________________

#include <cassert>
#include <string>

#include <TSystem.h>
#include <TNtupleD.h>
#include <TGraph2D.h>
#include <TTree.h>
#include <TMath.h>

#include "Framework/Conventions/Constants.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Physics/HadronTransport/INukeHadroData.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/Numerical/Spline.h"
#include "Framework/ParticleData/PDGCodes.h"

using std::ostringstream;
using std::ios;

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
INukeHadroData * INukeHadroData::fInstance = 0;
//____________________________________________________________________________
double INukeHadroData::fMinKinEnergy   =    1.0; // MeV
double INukeHadroData::fMaxKinEnergyHA =  999.0; // MeV
double INukeHadroData::fMaxKinEnergyHN = 1799.0; // MeV
//____________________________________________________________________________
INukeHadroData::INukeHadroData()
{
  this->LoadCrossSections();
  fInstance = 0;
}
//____________________________________________________________________________
INukeHadroData::~INukeHadroData()
{
  // pi+n/p hA x-section splines
  delete fXSecPipn_Tot;
  delete fXSecPipn_CEx;
  delete fXSecPipn_Elas;
  delete fXSecPipn_Reac;
  delete fXSecPipp_Tot;
  delete fXSecPipp_CEx;
  delete fXSecPipp_Elas;
  delete fXSecPipp_Reac;
  delete fXSecPipd_Abs;

  // pi0n/p hA x-section splines
  delete fXSecPi0n_Tot;
  delete fXSecPi0n_CEx;
  delete fXSecPi0n_Elas;
  delete fXSecPi0n_Reac;
  delete fXSecPi0p_Tot;
  delete fXSecPi0p_CEx;
  delete fXSecPi0p_Elas;
  delete fXSecPi0p_Reac;
  delete fXSecPi0d_Abs;

  // K+N x-section splines (elastic only)
  delete fXSecKpn_Elas;
  delete fXSecKpp_Elas;
  delete fXSecKpN_Abs;
  delete fXSecKpN_Tot;

  // gamma x-section splines (inelastic only)
  delete fXSecGamp_fs;
  delete fXSecGamn_fs;
  delete fXSecGamN_Tot;

  // N+A x-section splines
  delete fFracPA_Tot;
  delete fFracPA_Elas;
  delete fFracPA_Inel;
  delete fFracPA_CEx;
  delete fFracPA_Abs;
  delete fFracPA_Pipro;
  delete fFracNA_Tot;
  delete fFracNA_Elas;
  delete fFracNA_Inel;
  delete fFracNA_CEx;
  delete fFracNA_Abs;
  delete fFracNA_Pipro;

  // pi+A x-section splines
  delete fFracPipA_Tot;
  delete fFracPipA_Elas;
  delete fFracPipA_Inel;
  delete fFracPipA_CEx;
  delete fFracPipA_Abs;
  delete fFracPipA_PiProd;
  delete fFracPimA_Tot;
  delete fFracPimA_Elas;
  delete fFracPimA_Inel;
  delete fFracPimA_CEx;
  delete fFracPimA_Abs;
  delete fFracPimA_PiProd;
  delete fFracPi0A_Tot;
  delete fFracPi0A_Elas;
  delete fFracPi0A_Inel;
  delete fFracPi0A_CEx;
  delete fFracPi0A_Abs;
  delete fFracPi0A_PiProd;

  // hN data
  delete fhN2dXSecPP_Elas;
  delete fhN2dXSecNP_Elas;
  delete fhN2dXSecPipN_Elas;   
  delete fhN2dXSecPi0N_Elas;   
  delete fhN2dXSecPimN_Elas;   
  delete fhN2dXSecKpN_Elas;
  delete fhN2dXSecKpP_Elas;
  delete fhN2dXSecPiN_CEx;
  delete fhN2dXSecPiN_Abs;
  delete fhN2dXSecGamPi0P_Inelas;
  delete fhN2dXSecGamPi0N_Inelas;
  delete fhN2dXSecGamPipN_Inelas;
  delete fhN2dXSecGamPimP_Inelas;
  
  // p/n+p/n hA x-section splines
  delete fXSecPp_Tot;
  delete fXSecPp_Elas;      
  delete fXSecPp_Reac;      
  delete fXSecPn_Tot;     
  delete fXSecPn_Elas;      
  delete fXSecPn_Reac;      
  delete fXSecNn_Tot;     
  delete fXSecNn_Elas;      
  delete fXSecNn_Reac;      


  // K+A x-section fraction splines
  delete fFracKA_Tot;
  delete fFracKA_Elas;
  delete fFracKA_Inel;
  delete fFracKA_Abs;

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
void INukeHadroData::LoadCrossSections(void)
{
// Loads hadronic x-section data

  //-- Get the top-level directory with input hadron cross-section data
  //   (search for $GINUKEHADRONDATA or use default location)
  string data_dir = (gSystem->Getenv("GINUKEHADRONDATA")) ?
             string(gSystem->Getenv("GINUKEHADRONDATA")) :
             string(gSystem->Getenv("GENIE")) + string("/data/evgen/intranuke");

  LOG("INukeData", pINFO)  
      << "Loading INTRANUKE hadron data from: " << data_dir;

  //-- Build filenames

  string datafile_NN   = data_dir + "/tot_xsec/intranuke-xsections-NN.dat";
  string datafile_pipN = data_dir + "/tot_xsec/intranuke-xsections-pi+N.dat";
  string datafile_pi0N = data_dir + "/tot_xsec/intranuke-xsections-pi0N.dat";
  string datafile_NA   = data_dir + "/tot_xsec/intranuke-fractions-NA.dat";
  string datafile_piA  = data_dir + "/tot_xsec/intranuke-fractions-piA.dat";
  string datafile_KA   = data_dir + "/tot_xsec/intranuke-fractions-KA.dat";
  string datafile_gamN = data_dir + "/tot_xsec/intranuke-xsections-gamN.dat";
  string datafile_kN   = data_dir + "/tot_xsec/intranuke-xsections-kaonN.dat";

  //-- Make sure that all data files are available

  assert( ! gSystem->AccessPathName(datafile_NN.  c_str()) );
  assert( ! gSystem->AccessPathName(datafile_pipN.c_str()) );
  assert( ! gSystem->AccessPathName(datafile_pi0N.c_str()) );
  assert( ! gSystem->AccessPathName(datafile_NA.  c_str()) );
  assert( ! gSystem->AccessPathName(datafile_piA. c_str()) );
  assert( ! gSystem->AccessPathName(datafile_KA. c_str())  );
  assert( ! gSystem->AccessPathName(datafile_gamN.c_str())  );
  assert( ! gSystem->AccessPathName(datafile_kN.  c_str())  );

  LOG("INukeData", pINFO)  << "Found all necessary data files...";

  //-- Load data files

  TTree data_NN;
  TTree data_pipN;
  TTree data_pi0N;
  TTree data_NA;
  TTree data_piA;
  TTree data_KA;
  TTree data_gamN; 
  TTree data_kN;

  data_NN.ReadFile(datafile_NN.c_str(),
     "ke/D:pp_tot/D:pp_elas/D:pp_reac/D:pn_tot/D:pn_elas/D:pn_reac/D:nn_tot/D:nn_elas/D:nn_reac/D");
  data_pipN.ReadFile(datafile_pipN.c_str(),
     "ke/D:pipn_tot/D:pipn_cex/D:pipn_elas/D:pipn_reac/D:pipp_tot/D:pipp_cex/D:pipp_elas/D:pipp_reac/D:pipd_abs");
  data_pi0N.ReadFile(datafile_pi0N.c_str(),
     "ke/D:pi0n_tot/D:pi0n_cex/D:pi0n_elas/D:pi0n_reac/D:pi0p_tot/D:pi0p_cex/D:pi0p_elas/D:pi0p_reac/D:pi0d_abs");
  data_NA.ReadFile(datafile_NA.c_str(),
     "ke/D:pA_tot/D:pA_elas/D:pA_inel/D:pA_cex/D:pA_abs/D:pA_pipro/D");
  data_piA.ReadFile(datafile_piA.c_str(),
     "ke/D:piA_tot/D:piA_elas/D:piA_inel/D:piA_cex/D:piA_np/D:piA_pp/D:piA_npp/D:piA_nnp/D:piA_2n2p/D:piA_piprod/D");
  data_gamN.ReadFile(datafile_gamN.c_str(),
    "ke/D:pi0p_tot/D:pipn_tot/D:pimp_tot/D:pi0n_tot/D:gamp_fs/D:gamn_fs/D:gamN_tot/D");
  data_kN.ReadFile(datafile_kN.c_str(),
		   "ke/D:kpn_elas/D:kpp_elas/D:kp_abs/D:kpN_tot/D");  //????
  data_KA.ReadFile(datafile_KA.c_str(),
     "ke/D:KA_tot/D:KA_elas/D:KA_inel/D:KA_abs/D");

  LOG("INukeData", pDEBUG)  << "Number of data rows in NN : "   << data_NN.GetEntries();
  LOG("INukeData", pDEBUG)  << "Number of data rows in pipN : " << data_pipN.GetEntries();
  LOG("INukeData", pDEBUG)  << "Number of data rows in pi0N : " << data_pi0N.GetEntries();
  LOG("INukeData", pDEBUG)  << "Number of data rows in NA  : "  << data_NA.GetEntries();
  LOG("INukeData", pDEBUG)  << "Number of data rows in piA : "  << data_piA.GetEntries();
  LOG("INukeData", pDEBUG)  << "Number of data rows in KA : "   << data_KA.GetEntries();
  LOG("INukeData", pDEBUG)  << "Number of data rows in gamN : " << data_gamN.GetEntries(); 
  LOG("INukeData", pDEBUG)  << "Number of data rows in kN  : "  << data_kN.GetEntries();

  LOG("INukeData", pINFO)  << "Done loading all x-section files...";

  //-- Build x-section splines

  // p/n+p/n hA x-section splines
  fXSecPp_Tot      = new Spline(&data_NN, "ke:pp_tot");
  fXSecPp_Elas     = new Spline(&data_NN, "ke:pp_elas");      
  fXSecPp_Reac     = new Spline(&data_NN, "ke:pp_reac");      
  fXSecPn_Tot      = new Spline(&data_NN, "ke:pn_tot");     
  fXSecPn_Elas     = new Spline(&data_NN, "ke:pn_elas");      
  fXSecPn_Reac     = new Spline(&data_NN, "ke:pn_reac");      
  fXSecNn_Tot      = new Spline(&data_NN, "ke:nn_tot");     
  fXSecNn_Elas     = new Spline(&data_NN, "ke:nn_elas");      
  fXSecNn_Reac     = new Spline(&data_NN, "ke:nn_reac");      

  // pi+n/p hA x-section splines
  fXSecPipn_Tot     = new Spline(&data_pipN, "ke:pipn_tot");    
  fXSecPipn_CEx     = new Spline(&data_pipN, "ke:pipn_cex");    
  fXSecPipn_Elas    = new Spline(&data_pipN, "ke:pipn_elas");    
  fXSecPipn_Reac    = new Spline(&data_pipN, "ke:pipn_reac");    
  fXSecPipp_Tot     = new Spline(&data_pipN, "ke:pipp_tot");    
  fXSecPipp_CEx     = new Spline(&data_pipN, "ke:pipp_cex");    
  fXSecPipp_Elas    = new Spline(&data_pipN, "ke:pipp_elas");    
  fXSecPipp_Reac    = new Spline(&data_pipN, "ke:pipp_reac");    
  fXSecPipd_Abs     = new Spline(&data_pipN, "ke:pipd_abs");    

  // pi0n/p hA x-section splines
  fXSecPi0n_Tot     = new Spline(&data_pi0N, "ke:pi0n_tot");    
  fXSecPi0n_CEx     = new Spline(&data_pi0N, "ke:pi0n_cex");    
  fXSecPi0n_Elas    = new Spline(&data_pi0N, "ke:pi0n_elas");    
  fXSecPi0n_Reac    = new Spline(&data_pi0N, "ke:pi0n_reac");    
  fXSecPi0p_Tot     = new Spline(&data_pi0N, "ke:pi0p_tot");    
  fXSecPi0p_CEx     = new Spline(&data_pi0N, "ke:pi0p_cex");    
  fXSecPi0p_Elas    = new Spline(&data_pi0N, "ke:pi0p_elas");    
  fXSecPi0p_Reac    = new Spline(&data_pi0N, "ke:pi0p_reac");    
  fXSecPi0d_Abs     = new Spline(&data_pi0N, "ke:pi0d_abs");   

   // K+N x-section splines  
  fXSecKpn_Elas   = new Spline(&data_kN,  "ke:kpn_elas");
  fXSecKpp_Elas   = new Spline(&data_kN,  "ke:kpp_elas");
  fXSecKpN_Abs    = new Spline(&data_kN,  "ke:kp_abs");
  fXSecKpN_Tot    = new Spline(&data_kN,  "ke:kpN_tot");

  // gamma x-section splines  
  fXSecGamp_fs     = new Spline(&data_gamN, "ke:gamp_fs");
  fXSecGamn_fs     = new Spline(&data_gamN, "ke:gamn_fs");
  fXSecGamN_Tot    = new Spline(&data_gamN, "ke:gamN_tot");

  // N+A x-section fraction splines
  fFracPA_Tot      = new Spline(&data_NA, "ke:pA_tot");
  fFracPA_Elas     = new Spline(&data_NA, "ke:pA_elas");
  fFracPA_Inel     = new Spline(&data_NA, "ke:pA_inel");   
  fFracPA_CEx      = new Spline(&data_NA, "ke:pA_cex");   
  fFracPA_Abs      = new Spline(&data_NA, "ke:pA_abs");
  fFracPA_Pipro    = new Spline(&data_NA, "ke:pA_pipro");  
  fFracNA_Tot      = new Spline(&data_NA, "ke:pA_tot");  // assuming nA same as pA
  fFracNA_Elas     = new Spline(&data_NA, "ke:pA_elas"); 
  fFracNA_Inel     = new Spline(&data_NA, "ke:pA_inel");   
  fFracNA_CEx      = new Spline(&data_NA, "ke:pA_cex");   
  fFracNA_Abs      = new Spline(&data_NA, "ke:pA_abs");
  fFracNA_Pipro    = new Spline(&data_NA, "ke:pA_pipro");  

  // pi+A x-section splines
  fFracPipA_Tot     = new Spline(&data_piA, "ke:piA_tot");    
  fFracPipA_Elas    = new Spline(&data_piA, "ke:piA_elas");    
  fFracPipA_Inel    = new Spline(&data_piA, "ke:piA_inel");    
  fFracPipA_CEx     = new Spline(&data_piA, "ke:piA_cex");    
  fFracPipA_Abs     = new Spline(&data_piA, "ke:piA_np+piA_pp+piA_npp+piA_nnp+piA_2n2p");
  fFracPipA_PiProd  = new Spline(&data_piA, "ke:piA_piprod");    
  fFracPimA_Tot     = new Spline(&data_piA, "ke:piA_tot");    
  fFracPimA_Elas    = new Spline(&data_piA, "ke:piA_elas");    
  fFracPimA_Inel    = new Spline(&data_piA, "ke:piA_inel");    
  fFracPimA_CEx     = new Spline(&data_piA, "ke:piA_cex");    
  fFracPimA_Abs     = new Spline(&data_piA, "ke:piA_np+piA_pp+piA_npp+piA_nnp+piA_2n2p");
  fFracPimA_PiProd  = new Spline(&data_piA, "ke:piA_piprod");    
  fFracPi0A_Tot     = new Spline(&data_piA, "ke:piA_tot");    
  fFracPi0A_Elas    = new Spline(&data_piA, "ke:piA_elas");    
  fFracPi0A_Inel    = new Spline(&data_piA, "ke:piA_inel");    
  fFracPi0A_CEx     = new Spline(&data_piA, "ke:piA_cex");    
  fFracPi0A_Abs     = new Spline(&data_piA, "ke:piA_np+piA_pp+piA_npp+piA_nnp+piA_2n2p");
  fFracPi0A_PiProd  = new Spline(&data_piA, "ke:piA_piprod");
  // K+A x-section fraction splines
  fFracKA_Tot      = new Spline(&data_KA, "ke:KA_tot");
  fFracKA_Elas     = new Spline(&data_KA, "ke:KA_elas");
  fFracKA_Inel     = new Spline(&data_KA, "ke:KA_inel");
  fFracKA_Abs      = new Spline(&data_KA, "ke:KA_abs");
  //
  // hN stuff
  //


  // 	kIHNFtElas
  //	  pp, nn --> read from pp/pp%.txt
  //	  pn, np --> read from pp/pn%.txt
  //	  pi+ N  --> read from pip/pip%.txt
  //	  pi0 N  --> read from pip/pip%.txt
  //	  pi- N  --> read from pim/pim%.txt
  //      K+  N  --> read from kpn/kpn%.txt
  //      K+  P  --> read from kpp/kpp%.txt
  //    kIHNFtCEx
  //	  pi+, pi0, pi- --> read from pie/pie%.txt (using pip+n->pi0+p data)
  //    kIHNFtAbs
  //      pi+, pi0, pi- --> read from pid2p/pid2p%.txt (using pip+D->2p data)
  //    kIHNFtInelas
  //      gamma p -> p pi0 --> read from gampi0p/%-pi0p.txt
  //      gamma p -> n pi+ --> read from gampi+n/%-pi+n.txt
  //      gamma n -> n pi0 --> read from gampi0n/%-pi0n.txt
  //      gamma n -> p pi- --> read from gampi-p/%-pi-p.txt


  // kIHNFtElas, pp&nn :
  {
    const int hN_ppelas_nfiles = 20;
    const int hN_ppelas_points_per_file = 21;
    const int hN_ppelas_npoints = hN_ppelas_points_per_file * hN_ppelas_nfiles;

    double hN_ppelas_energies[hN_ppelas_nfiles] = {  
        50, 100, 150, 200, 250, 300, 350, 400, 450, 500,
       550, 600, 650, 700, 750, 800, 850, 900, 950, 1000 
    };

    double hN_ppelas_costh [hN_ppelas_points_per_file];
    double hN_ppelas_xsec  [hN_ppelas_npoints];

    int ipoint=0;

    for(int ifile = 0; ifile < hN_ppelas_nfiles; ifile++) {
     // build filename
     ostringstream hN_datafile;
     double ke = hN_ppelas_energies[ifile];
     hN_datafile << data_dir << "/diff_ang/pp/pp" << ke << ".txt";
     // read data
     ReadhNFile(
		hN_datafile.str(), ke, hN_ppelas_points_per_file, 
		ipoint, hN_ppelas_costh, hN_ppelas_xsec,2);
    }//loop over files

    /*double hN_ppelas_costh_cond [hN_ppelas_points_per_file];
    for (int ient = 0; ient < hN_ppelas_points_per_file; ient++) {
      hN_ppelas_costh_cond[ient] = hN_ppelas_costh[ient];
      }*/

    fhN2dXSecPP_Elas = new BLI2DNonUnifGrid(hN_ppelas_nfiles,hN_ppelas_points_per_file,
			   hN_ppelas_energies,hN_ppelas_costh,hN_ppelas_xsec); 
  }

  // kIHNFtElas, pn&np :
  {
    const int hN_npelas_nfiles = 20;
    const int hN_npelas_points_per_file = 21;
    const int hN_npelas_npoints = hN_npelas_points_per_file * hN_npelas_nfiles;

    double hN_npelas_energies[hN_npelas_nfiles] = {  
        50, 100, 150, 200, 250, 300, 350, 400, 450, 500,
       550, 600, 650, 700, 750, 800, 850, 900, 950, 1000 
    };

    double hN_npelas_costh [hN_npelas_points_per_file];
    double hN_npelas_xsec  [hN_npelas_npoints];

    int ipoint=0;

    for(int ifile = 0; ifile < hN_npelas_nfiles; ifile++) {
     // build filename
     ostringstream hN_datafile;
     double ke = hN_npelas_energies[ifile];
     hN_datafile << data_dir << "/diff_ang/pn/pn" << ke << ".txt";
     // read data
     ReadhNFile(
       hN_datafile.str(), ke, hN_npelas_points_per_file, 
       ipoint, hN_npelas_costh, hN_npelas_xsec,2);                
    }//loop over files

    /*double hN_npelas_costh_cond [hN_npelas_points_per_file];
    for (int ient = 0; ient < hN_npelas_points_per_file; ient++) {
      hN_npelas_costh_cond[ient] = hN_npelas_costh[ient];
      }*/

    fhN2dXSecNP_Elas = new BLI2DNonUnifGrid(hN_npelas_nfiles,hN_npelas_points_per_file,
			   hN_npelas_energies,hN_npelas_costh,hN_npelas_xsec); 
  }

  // kIHNFtElas, pipN :
  {
    const int hN_pipNelas_nfiles = 60;
    const int hN_pipNelas_points_per_file = 21;
    const int hN_pipNelas_npoints = hN_pipNelas_points_per_file * hN_pipNelas_nfiles;

    double hN_pipNelas_energies[hN_pipNelas_nfiles] = {  
      10,  20,  30,  40,  50,  60,  70,  80,  90, 
     100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 
     200, 210, 220, 230, 240, 250, 260, 270, 280, 290,
     300, 340, 380, 420, 460, 500, 540, 580, 620, 660, 
     700, 740, 780, 820, 860, 900, 940, 980, 
     1020, 1060, 1100, 1140, 1180, 1220, 1260, 
     1300, 1340, 1380, 1420, 1460, 1500
    };

    double hN_pipNelas_costh [hN_pipNelas_points_per_file];
    double hN_pipNelas_xsec  [hN_pipNelas_npoints];

    int ipoint=0;

    for(int ifile = 0; ifile < hN_pipNelas_nfiles; ifile++) {
     // build filename
     ostringstream hN_datafile;
     double ke = hN_pipNelas_energies[ifile];
     hN_datafile << data_dir << "/diff_ang/pip/pip" << ke << ".txt";
     // read data
     ReadhNFile(
       hN_datafile.str(), ke, hN_pipNelas_points_per_file, 
       ipoint, hN_pipNelas_costh, hN_pipNelas_xsec,2);                
    }//loop over files

    /*double hN_pipNelas_costh_cond [hN_pipNelas_points_per_file];
    for (int ient = 0; ient < hN_pipNelas_points_per_file; ient++) {
      hN_pipNelas_costh_cond[ient] = hN_pipNelas_costh[ient];
      }*/

    fhN2dXSecPipN_Elas = new BLI2DNonUnifGrid(hN_pipNelas_nfiles,hN_pipNelas_points_per_file,
			   hN_pipNelas_energies,hN_pipNelas_costh,hN_pipNelas_xsec); 
  }

  // kIHNFtElas, pi0N :
  {
    const int hN_pi0Nelas_nfiles = 60;
    const int hN_pi0Nelas_points_per_file = 21;
    const int hN_pi0Nelas_npoints = hN_pi0Nelas_points_per_file * hN_pi0Nelas_nfiles;

    double hN_pi0Nelas_energies[hN_pi0Nelas_nfiles] = {  
      10,  20,  30,  40,  50,  60,  70,  80,  90, 
     100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 
     200, 210, 220, 230, 240, 250, 260, 270, 280, 290,
     300, 340, 380, 420, 460, 500, 540, 580, 620, 660, 
     700, 740, 780, 820, 860, 900, 940, 980, 
     1020, 1060, 1100, 1140, 1180, 1220, 1260, 
     1300, 1340, 1380, 1420, 1460, 1500
    };

    double hN_pi0Nelas_costh [hN_pi0Nelas_points_per_file];
    double hN_pi0Nelas_xsec  [hN_pi0Nelas_npoints];

    int ipoint=0;

    for(int ifile = 0; ifile < hN_pi0Nelas_nfiles; ifile++) {
     // build filename
     ostringstream hN_datafile;
     double ke = hN_pi0Nelas_energies[ifile];
     hN_datafile << data_dir << "/diff_ang/pip/pip" << ke << ".txt";
     // read data
     ReadhNFile(
       hN_datafile.str(), ke, hN_pi0Nelas_points_per_file, 
       ipoint, hN_pi0Nelas_costh, hN_pi0Nelas_xsec,2);                
    }//loop over files

    /*double hN_pi0Nelas_costh_cond [hN_pi0Nelas_points_per_file];
    for (int ient = 0; ient < hN_pi0Nelas_points_per_file; ient++) {
      hN_pi0Nelas_costh_cond[ient] = hN_pi0Nelas_costh[ient];
      }*/

    fhN2dXSecPi0N_Elas = new BLI2DNonUnifGrid(hN_pi0Nelas_nfiles,hN_pi0Nelas_points_per_file,
			   hN_pi0Nelas_energies,hN_pi0Nelas_costh,hN_pi0Nelas_xsec); 
  }

  // kIHNFtElas, pimN :
  {
    const int hN_pimNelas_nfiles = 60;
    const int hN_pimNelas_points_per_file = 21;
    const int hN_pimNelas_npoints = hN_pimNelas_points_per_file * hN_pimNelas_nfiles;

    double hN_pimNelas_energies[hN_pimNelas_nfiles] = {  
      10,  20,  30,  40,  50,  60,  70,  80,  90, 
     100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 
     200, 210, 220, 230, 240, 250, 260, 270, 280, 290,
     300, 340, 380, 420, 460, 500, 540, 580, 620, 660, 
     700, 740, 780, 820, 860, 900, 940, 980, 
     1020, 1060, 1100, 1140, 1180, 1220, 1260, 
     1300, 1340, 1380, 1420, 1460, 1500
    };

    double hN_pimNelas_costh [hN_pimNelas_points_per_file];
    double hN_pimNelas_xsec  [hN_pimNelas_npoints];

    int ipoint=0;

    for(int ifile = 0; ifile < hN_pimNelas_nfiles; ifile++) {
     // build filename
     ostringstream hN_datafile;
     double ke = hN_pimNelas_energies[ifile];
     hN_datafile << data_dir << "/diff_ang/pim/pim" << ke << ".txt";
     // read data
     ReadhNFile(
       hN_datafile.str(), ke, hN_pimNelas_points_per_file, 
       ipoint, hN_pimNelas_costh, hN_pimNelas_xsec,2);                
    }//loop over files

    /*double hN_pimNelas_costh_cond [hN_pimNelas_points_per_file];
    for (int ient = 0; ient < hN_pimNelas_points_per_file; ient++) {
      hN_pimNelas_costh_cond[ient] = hN_pimNelas_costh[ient];
      }*/

    fhN2dXSecPimN_Elas = new BLI2DNonUnifGrid(hN_pimNelas_nfiles,hN_pimNelas_points_per_file,
			   hN_pimNelas_energies,hN_pimNelas_costh,hN_pimNelas_xsec); 
  }
 
  // kIHNFtElas, kpn :
  {
    const int hN_kpNelas_nfiles = 18;
    const int hN_kpNelas_points_per_file = 37;
    const int hN_kpNelas_npoints = hN_kpNelas_points_per_file * hN_kpNelas_nfiles;

    double hN_kpNelas_energies[hN_kpNelas_nfiles] = {  
     100, 200, 300, 400, 500, 600, 700, 800, 900, 1000,
     1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800
    };

    double hN_kpNelas_costh [hN_kpNelas_points_per_file];
    double hN_kpNelas_xsec  [hN_kpNelas_npoints];

    int ipoint=0;

    for(int ifile = 0; ifile < hN_kpNelas_nfiles; ifile++) {
     // build filename
     ostringstream hN_datafile;
     double ke = hN_kpNelas_energies[ifile];
     hN_datafile << data_dir << "/diff_ang/kpn/kpn" << ke << ".txt";
     // read data
     ReadhNFile(
       hN_datafile.str(), ke, hN_kpNelas_points_per_file, 
       ipoint, hN_kpNelas_costh, hN_kpNelas_xsec,2);                
    }//loop over files

    /*double hN_kpNelas_costh_cond [hN_kpNelas_points_per_file];
    for (int ient = 0; ient < hN_kpNelas_points_per_file; ient++) {
      hN_kpNelas_costh_cond[ient] = hN_kpNelas_costh[ient];
      }*/

    fhN2dXSecKpN_Elas = new BLI2DNonUnifGrid(hN_kpNelas_nfiles,hN_kpNelas_points_per_file,
			   hN_kpNelas_energies,hN_kpNelas_costh,hN_kpNelas_xsec); 
  }
  
  // kIHNFtElas, kpp :
  {
    const int hN_kpPelas_nfiles = 18;
    const int hN_kpPelas_points_per_file = 37;
    const int hN_kpPelas_npoints = hN_kpPelas_points_per_file * hN_kpPelas_nfiles;

    double hN_kpPelas_energies[hN_kpPelas_nfiles] = {  
     100, 200, 300, 400, 500, 600, 700, 800, 900, 1000,
     1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800
    };

    double hN_kpPelas_costh [hN_kpPelas_points_per_file];
    double hN_kpPelas_xsec  [hN_kpPelas_npoints];

    int ipoint=0;

    for(int ifile = 0; ifile < hN_kpPelas_nfiles; ifile++) {
     // build filename
     ostringstream hN_datafile;
     double ke = hN_kpPelas_energies[ifile];
     hN_datafile << data_dir << "/diff_ang/kpp/kpp" << ke << ".txt";
     // read data
     ReadhNFile(
       hN_datafile.str(), ke, hN_kpPelas_points_per_file, 
       ipoint, hN_kpPelas_costh, hN_kpPelas_xsec,2);                
    }//loop over files

    /*double hN_kpPelas_costh_cond [hN_kpPelas_points_per_file];
    for (int ient = 0; ient < hN_kpPelas_points_per_file; ient++) {
      hN_kpPelas_costh_cond[ient] = hN_kpPelas_costh[ient];
      }*/

    fhN2dXSecKpP_Elas = new BLI2DNonUnifGrid(hN_kpPelas_nfiles,hN_kpPelas_points_per_file,
			   hN_kpPelas_energies,hN_kpPelas_costh,hN_kpPelas_xsec); 
	}

  // kIHNFtCEx, (pi+, pi0, pi-) N 
  {
    const int hN_piNcex_nfiles = 60;
    const int hN_piNcex_points_per_file = 21;
    const int hN_piNcex_npoints = hN_piNcex_points_per_file * hN_piNcex_nfiles;

    double hN_piNcex_energies[hN_piNcex_nfiles] = {  
      10,  20,  30,  40,  50,  60,  70,  80,  90, 
     100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 
     200, 210, 220, 230, 240, 250, 260, 270, 280, 290,
     300, 340, 380, 420, 460, 500, 540, 580, 620, 660, 
     700, 740, 780, 820, 860, 900, 940, 980, 
     1020, 1060, 1100, 1140, 1180, 1220, 1260, 
     1300, 1340, 1380, 1420, 1460, 1500
    };

    double hN_piNcex_costh [hN_piNcex_points_per_file];
    double hN_piNcex_xsec  [hN_piNcex_npoints];

    int ipoint=0;

    for(int ifile = 0; ifile < hN_piNcex_nfiles; ifile++) {
     // build filename
     ostringstream hN_datafile;
     double ke = hN_piNcex_energies[ifile];
     hN_datafile << data_dir << "/diff_ang/pie/pie" << ke << ".txt";
     // read data
     ReadhNFile(
       hN_datafile.str(), ke, hN_piNcex_points_per_file, 
       ipoint, hN_piNcex_costh, hN_piNcex_xsec,2);                
    }//loop over files

    /*double hN_piNcex_costh_cond [hN_piNcex_points_per_file];
    for (int ient = 0; ient < hN_piNcex_points_per_file; ient++) {
      hN_piNcex_costh_cond[ient] = hN_piNcex_costh[ient];
      }*/

    fhN2dXSecPiN_CEx = new BLI2DNonUnifGrid(hN_piNcex_nfiles,hN_piNcex_points_per_file,
			   hN_piNcex_energies,hN_piNcex_costh,hN_piNcex_xsec); 
  }

  // kIHNFtAbs, (pi+, pi0, pi-) N 
  {
    const int hN_piNabs_nfiles = 19;
    const int hN_piNabs_points_per_file = 21;
    const int hN_piNabs_npoints = hN_piNabs_points_per_file * hN_piNabs_nfiles;

    double hN_piNabs_energies[hN_piNabs_nfiles] = {  
      50,  75, 100, 125, 150, 175, 200, 225, 250, 275,  
     300, 325, 350, 375, 400, 425, 450, 475, 500
    };

    double hN_piNabs_costh [hN_piNabs_points_per_file];
    double hN_piNabs_xsec  [hN_piNabs_npoints];

    int ipoint=0;

    for(int ifile = 0; ifile < hN_piNabs_nfiles; ifile++) {
     // build filename
     ostringstream hN_datafile;
     double ke = hN_piNabs_energies[ifile];
     hN_datafile << data_dir << "/diff_ang/pid2p/pid2p" << ke << ".txt";
     // read data
     ReadhNFile(
       hN_datafile.str(), ke, hN_piNabs_points_per_file, 
       ipoint, hN_piNabs_costh, hN_piNabs_xsec,2);                
    }//loop over files

    /*double hN_piNabs_costh_cond [hN_piNabs_points_per_file];
    for (int ient = 0; ient < hN_piNabs_points_per_file; ient++) {
      hN_piNabs_costh_cond[ient] = hN_piNabs_costh[ient];
      }*/

    fhN2dXSecPiN_Abs = new BLI2DNonUnifGrid(hN_piNabs_nfiles,hN_piNabs_points_per_file,
			   hN_piNabs_energies,hN_piNabs_costh,hN_piNabs_xsec);
  }

  // kIHNFtInelas, gamma p -> p pi0
    {
    const int hN_gampi0pInelas_nfiles = 29;
    const int hN_gampi0pInelas_points_per_file = 37;
    const int hN_gampi0pInelas_npoints = hN_gampi0pInelas_points_per_file * hN_gampi0pInelas_nfiles;

    double hN_gampi0pInelas_energies[hN_gampi0pInelas_nfiles] = {  
      160,  180,  200,  220,  240,  260,  280,  300,  320,  340,   
      360,  380,  400,  450,  500,  550,  600,  650,  700,  750,
      800,  850,  900,  950,  1000, 1050, 1100, 1150, 1200
    };

    double hN_gampi0pInelas_costh [hN_gampi0pInelas_points_per_file];
    double hN_gampi0pInelas_xsec  [hN_gampi0pInelas_npoints];

    int ipoint=0;

    for(int ifile = 0; ifile < hN_gampi0pInelas_nfiles; ifile++) {
     // build filename
     ostringstream hN_datafile;
     double ke = hN_gampi0pInelas_energies[ifile];
     hN_datafile << data_dir << "/diff_ang/gampi0p/" << ke << "-pi0p.txt";
     // read data
     ReadhNFile(
       hN_datafile.str(), ke, hN_gampi0pInelas_points_per_file, 
       ipoint, hN_gampi0pInelas_costh, hN_gampi0pInelas_xsec,3);                
    }//loop over files

    /*double hN_gampi0pInelas_costh_cond [hN_gampi0pInelas_points_per_file];
    for (int ient = 0; ient < hN_gampi0pInelas_points_per_file; ient++) {
      hN_gampi0pInelas_costh_cond[ient] = hN_gampi0pInelas_costh[ient];
      }*/

    fhN2dXSecGamPi0P_Inelas = new BLI2DNonUnifGrid(hN_gampi0pInelas_nfiles,hN_gampi0pInelas_points_per_file,
			   hN_gampi0pInelas_energies,hN_gampi0pInelas_costh,hN_gampi0pInelas_xsec);
  }

  // kIHNFtInelas, gamma n -> n pi0
  {
    const int hN_gampi0nInelas_nfiles = 29;
    const int hN_gampi0nInelas_points_per_file = 37;
    const int hN_gampi0nInelas_npoints = hN_gampi0nInelas_points_per_file * hN_gampi0nInelas_nfiles;

    double hN_gampi0nInelas_energies[hN_gampi0nInelas_nfiles] = {  
      160,  180,  200,  220,  240,  260,  280,  300,  320,  340,   
      360,  380,  400,  450,  500,  550,  600,  650,  700,  750,
      800,  850,  900,  950,  1000, 1050, 1100, 1150, 1200
    };

    double hN_gampi0nInelas_costh [hN_gampi0nInelas_points_per_file];
    double hN_gampi0nInelas_xsec  [hN_gampi0nInelas_npoints];
    int ipoint=0;

    for(int ifile = 0; ifile < hN_gampi0nInelas_nfiles; ifile++) {
     // build filename
     ostringstream hN_datafile;
     double ke = hN_gampi0nInelas_energies[ifile];
     hN_datafile << data_dir << "/diff_ang/gampi0n/" << ke << "-pi0n.txt";
     // read data
     ReadhNFile(
       hN_datafile.str(), ke, hN_gampi0nInelas_points_per_file, 
       ipoint, hN_gampi0nInelas_costh, hN_gampi0nInelas_xsec,3);                
    }//loop over files

    /*double hN_gampi0nInelas_costh_cond [hN_gampi0nInelas_points_per_file];
    for (int ient = 0; ient < hN_gampi0nInelas_points_per_file; ient++) {
      hN_gampi0nInelas_costh_cond[ient] = hN_gampi0nInelas_costh[ient];
      }*/

    fhN2dXSecGamPi0N_Inelas = new BLI2DNonUnifGrid(hN_gampi0nInelas_nfiles,hN_gampi0nInelas_points_per_file,
			   hN_gampi0nInelas_energies,hN_gampi0nInelas_costh,hN_gampi0nInelas_xsec);
  }

  // kIHNFtInelas, gamma p -> n pi+
  {
    const int hN_gampipnInelas_nfiles = 29;
    const int hN_gampipnInelas_points_per_file = 37;
    const int hN_gampipnInelas_npoints = hN_gampipnInelas_points_per_file * hN_gampipnInelas_nfiles;

    double hN_gampipnInelas_energies[hN_gampipnInelas_nfiles] = {  
      160,  180,  200,  220,  240,  260,  280,  300,  320,  340,   
      360,  380,  400,  450,  500,  550,  600,  650,  700,  750,
      800,  850,  900,  950,  1000, 1050, 1100, 1150, 1200
    };

    double hN_gampipnInelas_costh [hN_gampipnInelas_points_per_file];
    double hN_gampipnInelas_xsec  [hN_gampipnInelas_npoints];

    int ipoint=0;

    for(int ifile = 0; ifile < hN_gampipnInelas_nfiles; ifile++) {
     // build filename
     ostringstream hN_datafile;
     double ke = hN_gampipnInelas_energies[ifile];
     hN_datafile << data_dir << "/diff_ang/gampi+n/" << ke << "-pi+n.txt";
     // read data
     ReadhNFile(
       hN_datafile.str(), ke, hN_gampipnInelas_points_per_file, 
       ipoint, hN_gampipnInelas_costh, hN_gampipnInelas_xsec,3);                
    }//loop over files

    /*double hN_gampipnInelas_costh_cond [hN_gampipnInelas_points_per_file];
    for (int ient = 0; ient < hN_gampipnInelas_points_per_file; ient++) {
      hN_gampipnInelas_costh_cond[ient] = hN_gampipnInelas_costh[ient];
      }*/

    fhN2dXSecGamPipN_Inelas = new BLI2DNonUnifGrid(hN_gampipnInelas_nfiles,hN_gampipnInelas_points_per_file,
			   hN_gampipnInelas_energies,hN_gampipnInelas_costh,hN_gampipnInelas_xsec);
  }

  // kIHNFtInelas, gamma n -> p pi-
  {
    const int hN_gampimpInelas_nfiles = 29;
    const int hN_gampimpInelas_points_per_file = 37;
    const int hN_gampimpInelas_npoints = hN_gampimpInelas_points_per_file * hN_gampimpInelas_nfiles;

    double hN_gampimpInelas_energies[hN_gampimpInelas_nfiles] = {  
      160,  180,  200,  220,  240,  260,  280,  300,  320,  340,   
      360,  380,  400,  450,  500,  550,  600,  650,  700,  750,
      800,  850,  900,  950,  1000, 1050, 1100, 1150, 1200
    };

    double hN_gampimpInelas_costh [hN_gampimpInelas_points_per_file];
    double hN_gampimpInelas_xsec  [hN_gampimpInelas_npoints];

    int ipoint=0;

    for(int ifile = 0; ifile < hN_gampimpInelas_nfiles; ifile++) {
     // build filename
     ostringstream hN_datafile;
     double ke = hN_gampimpInelas_energies[ifile];
     hN_datafile << data_dir << "/diff_ang/gampi-p/" << ke << "-pi-p.txt";
     // read data
     ReadhNFile(
       hN_datafile.str(), ke, hN_gampimpInelas_points_per_file, 
       ipoint, hN_gampimpInelas_costh, hN_gampimpInelas_xsec,3);                
    }//loop over files

    /*double hN_gampimpInelas_costh_cond [hN_gampimpInelas_points_per_file];
    for (int ient = 0; ient < hN_gampimpInelas_points_per_file; ient++) {
      hN_gampimpInelas_costh_cond[ient] = hN_gampimpInelas_costh[ient];
      }*/

    fhN2dXSecGamPimP_Inelas = new BLI2DNonUnifGrid(hN_gampimpInelas_nfiles,hN_gampimpInelas_points_per_file,
			   hN_gampimpInelas_energies,hN_gampimpInelas_costh,hN_gampimpInelas_xsec);
  }

  LOG("INukeData", pINFO)  << "Done building x-section splines...";
}
//____________________________________________________________________________
void INukeHadroData::ReadhNFile(
  string filename, double ke, int npoints, int & curr_point,
  double * costh_array, double * xsec_array, int cols)
{
  // open 
  std::ifstream hN_stream(filename.c_str(), ios::in);
  if(!hN_stream.good()) {
      LOG("INukeData", pERROR) 
          << "Error reading INTRANUKE/hN data from: " << filename;
      return;
  }

  if(cols<2) {
    LOG("INukeData", pERROR)
      << "Error reading INTRANUKE/hN data from: " << filename;
    LOG("INukeData", pERROR)
      << "Too few columns: " << cols;
    return;
  }

  LOG("INukeData", pINFO)  
     << "Reading INTRANUKE/hN data from: " << filename;

  // skip initial comments
  char cbuf[501];
  hN_stream.getline(cbuf,400);
  hN_stream.getline(cbuf,400);
  hN_stream.getline(cbuf,400);

  // read
  double angle = 0;
  double xsec  = 0;
  double trash = 0;

  for(int ip = 0; ip < npoints; ip++) {
     hN_stream >> angle >> xsec;

     for(int ic = 0; ic < (cols-2); ic++) {
       hN_stream >> trash;
     }

     LOG("INukeData", pDEBUG)  
       << "Adding data point: (KE = " << ke << " MeV, angle = " 
       << angle << ", sigma = " << xsec << " mbarn)";
     costh_array[ip] = TMath::Cos(angle*kPi/180.);
     xsec_array [curr_point] = xsec;
     curr_point++;
  }
}
//____________________________________________________________________________
double INukeHadroData::XSec(
  int hpdgc, int tgtpdgc, int nppdgc, INukeFateHN_t fate, double ke, double costh) const
{
// inputs
//      fate    : h+N fate code
//      hpdgc   : h PDG code
//      tgtpdgc : N PDG code
//      nppdgc  : product N PDG code
//      ke      : kinetic energy (MeV)
//      costh   : cos(scattering angle)
// returns
//      xsec    : mbarn

  double ke_eval    = ke;
  double costh_eval = costh;

  costh_eval = TMath::Min(costh,  1.);
  costh_eval = TMath::Max(costh_eval, -1.);

  if(fate==kIHNFtElas) {

     if( (hpdgc==kPdgProton  && tgtpdgc==kPdgProton) ||
         (hpdgc==kPdgNeutron && tgtpdgc==kPdgNeutron) )
     {
       ke_eval = TMath::Min(ke_eval, 999.);
       ke_eval = TMath::Max(ke_eval,  50.);
       return fhN2dXSecPP_Elas->Evaluate(ke_eval, costh_eval);
     } 
     else 
     if( (hpdgc==kPdgProton  && tgtpdgc==kPdgNeutron) ||
         (hpdgc==kPdgNeutron && tgtpdgc==kPdgProton) )
     {
       ke_eval = TMath::Min(ke_eval, 999.);
       ke_eval = TMath::Max(ke_eval,  50.);
       return fhN2dXSecNP_Elas->Evaluate(ke_eval, costh_eval);
     } 
     else
     if(hpdgc==kPdgPiP) 
     {
       ke_eval = TMath::Min(ke_eval, 1499.);
       ke_eval = TMath::Max(ke_eval,   10.);
       return fhN2dXSecPipN_Elas->Evaluate(ke_eval, costh_eval);
     } 
     else
     if(hpdgc==kPdgPi0) 
     {
       ke_eval = TMath::Min(ke_eval, 1499.);
       ke_eval = TMath::Max(ke_eval,   10.);
       return fhN2dXSecPi0N_Elas->Evaluate(ke_eval, costh_eval);
     } 
     else
     if(hpdgc==kPdgPiM) 
     {
       ke_eval = TMath::Min(ke_eval, 1499.);
       ke_eval = TMath::Max(ke_eval,   10.);
       return fhN2dXSecPimN_Elas->Evaluate(ke_eval, costh_eval);
     }
     else
     if(hpdgc==kPdgKP && tgtpdgc==kPdgNeutron) 
     {
       ke_eval = TMath::Min(ke_eval, 1799.);
       ke_eval = TMath::Max(ke_eval,  100.);
       return fhN2dXSecKpN_Elas->Evaluate(ke_eval, costh_eval);
     }
     else
     if(hpdgc==kPdgKP && tgtpdgc==kPdgProton) 
     {
       ke_eval = TMath::Min(ke_eval, 1799.);
       ke_eval = TMath::Max(ke_eval,  100.);
       return fhN2dXSecKpP_Elas->Evaluate(ke_eval, costh_eval);
     }
  }

  else if(fate == kIHNFtCEx) {
    if( (hpdgc==kPdgPiP || hpdgc==kPdgPi0 || hpdgc==kPdgPiM) &&
         (tgtpdgc==kPdgProton || tgtpdgc==kPdgNeutron) )
     {
        ke_eval = TMath::Min(ke_eval, 1499.);
        ke_eval = TMath::Max(ke_eval,   10.);
        return fhN2dXSecPiN_CEx->Evaluate(ke_eval, costh_eval);
     }
    else if( (hpdgc == kPdgProton && tgtpdgc == kPdgProton) ||
	     (hpdgc == kPdgNeutron && tgtpdgc == kPdgNeutron) )
      {
	LOG("INukeData", pWARN)  << "Inelastic pp does not exist!";
	ke_eval = TMath::Min(ke_eval, 999.);
	ke_eval = TMath::Max(ke_eval,  50.);
	return fhN2dXSecPP_Elas->Evaluate(ke_eval, costh_eval);
      }
    else if( (hpdgc == kPdgProton && tgtpdgc == kPdgNeutron) ||
	     (hpdgc == kPdgNeutron && tgtpdgc == kPdgProton) )
      {
	ke_eval = TMath::Min(ke_eval, 999.);
	ke_eval = TMath::Max(ke_eval,  50.);
	return fhN2dXSecNP_Elas->Evaluate(ke_eval, costh_eval);
      }
  }

  else if(fate == kIHNFtAbs) {
    if( (hpdgc==kPdgPiP || hpdgc==kPdgPi0 || hpdgc==kPdgPiM) &&
         (tgtpdgc==kPdgProton || tgtpdgc==kPdgNeutron) )
     {
        ke_eval = TMath::Min(ke_eval, 499.);
        ke_eval = TMath::Max(ke_eval,  50.);
        return fhN2dXSecPiN_Abs->Evaluate(ke_eval, costh_eval);
     }
    if(hpdgc==kPdgKP) return 1.;  //isotropic since no data ???
  }

  else if(fate == kIHNFtInelas) {
    if( hpdgc==kPdgGamma && tgtpdgc==kPdgProton  &&nppdgc==kPdgProton  )
    {
       ke_eval = TMath::Min(ke_eval, 1199.);
       ke_eval = TMath::Max(ke_eval,  160.);
       return fhN2dXSecGamPi0P_Inelas->Evaluate(ke_eval, costh_eval);
    }
    else
    if( hpdgc==kPdgGamma && tgtpdgc==kPdgProton  && nppdgc==kPdgNeutron )
    {
       ke_eval = TMath::Min(ke_eval, 1199.);
       ke_eval = TMath::Max(ke_eval,  160.);
       return fhN2dXSecGamPipN_Inelas->Evaluate(ke_eval, costh_eval);
    }
    else
    if( hpdgc==kPdgGamma && tgtpdgc==kPdgNeutron && nppdgc==kPdgProton  )
    {
       ke_eval = TMath::Min(ke_eval, 1199.);
       ke_eval = TMath::Max(ke_eval,  160.);
       return fhN2dXSecGamPimP_Inelas->Evaluate(ke_eval, costh_eval);
    }
    else
    if( hpdgc==kPdgGamma && tgtpdgc==kPdgNeutron && nppdgc==kPdgNeutron )
    {
       ke_eval = TMath::Min(ke_eval, 1199.);
       ke_eval = TMath::Max(ke_eval,  160.);
      return fhN2dXSecGamPi0N_Inelas->Evaluate(ke_eval, costh_eval);
    }
  }

  return 0;  
}
//____________________________________________________________________________
double INukeHadroData::Frac(int hpdgc, INukeFateHA_t fate, double ke) const
{
// return the x-section fraction for the input fate for the particle with the input pdg 
// code at the input kinetic energy
//
  ke = TMath::Max(fMinKinEnergy,   ke);
  ke = TMath::Min(fMaxKinEnergyHA, ke);

  LOG("INukeData", pDEBUG)  << "Querying hA cross section at ke = " << ke;

  if(hpdgc == kPdgProton) {
   /* handle protons */
        if (fate == kIHAFtCEx    ) return TMath::Max(0., fFracPA_CEx     -> Evaluate (ke));
   else if (fate == kIHAFtElas   ) return TMath::Max(0., fFracPA_Elas    -> Evaluate (ke));
   else if (fate == kIHAFtInelas ) return TMath::Max(0., fFracPA_Inel    -> Evaluate (ke));
   else if (fate == kIHAFtAbs    ) return TMath::Max(0., fFracPA_Abs     -> Evaluate (ke));
   else if (fate == kIHAFtPiProd ) return TMath::Max(0., fFracPA_Pipro   -> Evaluate (ke));
   else {
     LOG("INukeData", pWARN) 
       << "Protons don't have this fate: " << INukeHadroFates::AsString(fate);
     return 0;
   }

  } else if (hpdgc == kPdgNeutron) {
   /* handle neutrons */
        if (fate == kIHAFtCEx    ) return TMath::Max(0., fFracNA_CEx     -> Evaluate (ke));
   else if (fate == kIHAFtElas   ) return TMath::Max(0., fFracNA_Elas    -> Evaluate (ke));
   else if (fate == kIHAFtInelas ) return TMath::Max(0., fFracNA_Inel    -> Evaluate (ke));
   else if (fate == kIHAFtAbs    ) return TMath::Max(0., fFracNA_Abs     -> Evaluate (ke));
   else if (fate == kIHAFtPiProd ) return TMath::Max(0., fFracNA_Pipro   -> Evaluate (ke));
   else {
     LOG("INukeData", pWARN) 
       << "Neutrons don't have this fate: " << INukeHadroFates::AsString(fate);
     return 0;
   }

  } else if (hpdgc == kPdgPiP) {
   /* handle pi+ */
        if (fate == kIHAFtCEx    ) return TMath::Max(0., fFracPipA_CEx     -> Evaluate (ke));
   else if (fate == kIHAFtElas   ) return TMath::Max(0., fFracPipA_Elas    -> Evaluate (ke));
   else if (fate == kIHAFtInelas ) return TMath::Max(0., fFracPipA_Inel    -> Evaluate (ke));
   else if (fate == kIHAFtAbs    ) return TMath::Max(0., fFracPipA_Abs     -> Evaluate (ke));
	//   else if (fate == kIHAFtPiProd ) return TMath::Max(0., fFracPipA_Pipro   -> Evaluate (ke));
   else if (fate == kIHAFtPiProd ) return TMath::Max(0., fFracPipA_PiProd -> Evaluate (ke));
   else if (fate == kIHAFtPiProd) return TMath::Max(0., fFracPipA_PiProd -> Evaluate (ke));
   else {
     LOG("INukeData", pWARN) 
         << "Pi+'s don't have this fate: " << INukeHadroFates::AsString(fate);
     return 0;
   }

  } else if (hpdgc == kPdgPiM) {
   /* handle pi- */
   if      (fate == kIHAFtCEx    ) return TMath::Max(0., fFracPimA_CEx     -> Evaluate (ke));
   else if (fate == kIHAFtElas   ) return TMath::Max(0., fFracPimA_Elas    -> Evaluate (ke));
   else if (fate == kIHAFtInelas ) return TMath::Max(0., fFracPimA_Inel    -> Evaluate (ke));
   else if (fate == kIHAFtAbs    ) return TMath::Max(0., fFracPimA_Abs     -> Evaluate (ke));
   //   else if (fate == kIHAFtPiProd ) return TMath::Max(0., fFracPimA_Pipro   -> Evaluate (ke));
   else if (fate == kIHAFtPiProd ) return TMath::Max(0., fFracPimA_PiProd -> Evaluate (ke));
   else if (fate == kIHAFtPiProd) return TMath::Max(0., fFracPimA_PiProd -> Evaluate (ke));
   else {
     LOG("INukeData", pWARN) 
        << "Pi-'s don't have this fate: " << INukeHadroFates::AsString(fate);
     return 0;
   }

  } else if (hpdgc == kPdgPi0) {
   /* handle pi0 */
        if (fate == kIHAFtCEx    ) return TMath::Max(0., fFracPi0A_CEx     -> Evaluate (ke));
   else if (fate == kIHAFtElas   ) return TMath::Max(0., fFracPi0A_Elas    -> Evaluate (ke));
   else if (fate == kIHAFtInelas ) return TMath::Max(0., fFracPi0A_Inel    -> Evaluate (ke));
   else if (fate == kIHAFtAbs    ) return TMath::Max(0., fFracPi0A_Abs     -> Evaluate (ke));
	//   else if (fate == kIHAFtPiProd ) return TMath::Max(0., fFracPi0A_Pipro   -> Evaluate (ke));
   else if (fate == kIHAFtPiProd ) return TMath::Max(0., fFracPi0A_PiProd -> Evaluate (ke));
   else if (fate == kIHAFtPiProd) return TMath::Max(0., fFracPi0A_PiProd -> Evaluate (ke));
   else {
     LOG("INukeData", pWARN) 
        << "Pi0's don't have this fate: " << INukeHadroFates::AsString(fate);
       return 0;
   }
  } else if (hpdgc == kPdgKP) {
   /* handle K+ */
  if (fate == kIHAFtInelas ) return TMath::Max(0., fFracKA_Inel    -> Evaluate (ke));
  else if (fate == kIHAFtAbs    ) return TMath::Max(0., fFracKA_Abs     -> Evaluate (ke));
  //  else if (fate == kIHAFtElas   ) return TMath::Max(0., fFracKA_Elas    -> Evaluate (ke));
  else {
    LOG("INukeData", pWARN) 
      << "K+'s don't have this fate: " << INukeHadroFates::AsString(fate);
       return 0;
   }
  }
  LOG("INukeData", pWARN) 
      << "Can't handle particles with pdg code = " << hpdgc;

  return 0;
}
//____________________________________________________________________________
double INukeHadroData::XSec(int hpdgc, INukeFateHN_t fate, double ke, int targA, int targZ) const
{
// return the x-section for the input fate for the particle with the input pdg 
// code at the input kinetic energy
//
  ke = TMath::Max(fMinKinEnergy,   ke);
  ke = TMath::Min(fMaxKinEnergyHN, ke);

  LOG("INukeData", pDEBUG)  << "Querying hN cross section at ke = " << ke;

  double xsec=0;

    if (hpdgc == kPdgPiP) {
    /* handle pi+ */
         if (fate == kIHNFtCEx   ) {xsec = TMath::Max(0., fXSecPipp_CEx  -> Evaluate(ke)) *  targZ;
	                            xsec+= TMath::Max(0., fXSecPipn_CEx  -> Evaluate(ke)) * (targA-targZ);
				    return xsec;}
    else if (fate == kIHNFtElas  ) {xsec = TMath::Max(0., fXSecPipp_Elas -> Evaluate(ke)) *  targZ;
	                            xsec+= TMath::Max(0., fXSecPipn_Elas -> Evaluate(ke)) * (targA-targZ);
				    return xsec;}
    else if (fate == kIHNFtInelas) {xsec = TMath::Max(0., fXSecPipp_Reac -> Evaluate(ke)) *  targZ;
				    xsec+= TMath::Max(0., fXSecPipn_Reac -> Evaluate(ke)) * (targA-targZ);
 				    return xsec;}
    else if (fate == kIHNFtAbs   ) {xsec = TMath::Max(0., fXSecPipd_Abs  -> Evaluate(ke)) *  targA;
				    return xsec;}
    else {
     LOG("INukeData", pWARN) 
        << "Pi+'s don't have this fate: " << INukeHadroFates::AsString(fate);
     return 0;
    }

  } else if (hpdgc == kPdgPiM) {
    /* handle pi- */
         if (fate == kIHNFtCEx   ) {xsec = TMath::Max(0., fXSecPipn_CEx  -> Evaluate(ke)) *  targZ;
	                            xsec+= TMath::Max(0., fXSecPipp_CEx  -> Evaluate(ke)) * (targA-targZ);
				    return xsec;}
    else if (fate == kIHNFtElas  ) {xsec = TMath::Max(0., fXSecPipn_Elas -> Evaluate(ke)) *  targZ;
	                            xsec+= TMath::Max(0., fXSecPipp_Elas -> Evaluate(ke)) * (targA-targZ);
				    return xsec;}
    else if (fate == kIHNFtInelas) {xsec = TMath::Max(0., fXSecPipn_Reac -> Evaluate(ke)) *  targZ;
	                            xsec+= TMath::Max(0., fXSecPipp_Reac -> Evaluate(ke)) * (targA-targZ);
				    return xsec;}
    else if (fate == kIHNFtAbs   ) {xsec = TMath::Max(0., fXSecPipd_Abs  -> Evaluate(ke)) *  targA;
				    return xsec;}
    else {
     LOG("INukeData", pWARN) 
        << "Pi-'s don't have this fate: " << INukeHadroFates::AsString(fate);
     return 0;
    }

  } else if (hpdgc == kPdgPi0) {
    /* handle pi0 */
         if (fate == kIHNFtCEx   ) {xsec = TMath::Max(0., fXSecPi0p_CEx  -> Evaluate(ke)) *  targZ;
	                            xsec+= TMath::Max(0., fXSecPi0n_CEx  -> Evaluate(ke)) * (targA-targZ);
				    return xsec;}
    else if (fate == kIHNFtElas  ) {xsec = TMath::Max(0., fXSecPi0p_Elas -> Evaluate(ke)) *  targZ;
	                            xsec+= TMath::Max(0., fXSecPi0n_Elas -> Evaluate(ke)) * (targA-targZ);
				    return xsec;}
    else if (fate == kIHNFtInelas) {xsec = TMath::Max(0., fXSecPi0p_Reac -> Evaluate(ke)) *  targZ;
	                            xsec+= TMath::Max(0., fXSecPi0n_Reac -> Evaluate(ke)) * (targA-targZ);
				    return xsec;}
    else if (fate == kIHNFtAbs   ) {xsec = TMath::Max(0., fXSecPi0d_Abs  -> Evaluate(ke)) *  targA;
				    return xsec;}
    else {
     LOG("INukeData", pWARN) 
        << "Pi0's don't have this fate: " << INukeHadroFates::AsString(fate);
     return 0;
    }

  } else if (hpdgc == kPdgProton) {
    /* handle protons */
         if (fate == kIHNFtElas  ) {xsec = TMath::Max(0., fXSecPp_Elas -> Evaluate(ke)) *  targZ;
	                            xsec+= TMath::Max(0., fXSecPn_Elas -> Evaluate(ke)) * (targA-targZ);
				    return xsec;}
    else if (fate == kIHNFtInelas) {xsec = TMath::Max(0., fXSecPp_Reac -> Evaluate(ke)) *  targZ;
	                            xsec+= TMath::Max(0., fXSecPn_Reac -> Evaluate(ke)) * (targA-targZ);
				    return xsec;}
    else {
     LOG("INukeData", pWARN) 
        << "Protons don't have this fate: " << INukeHadroFates::AsString(fate);
     return 0;
    }

  } else if (hpdgc == kPdgNeutron) {
    /* handle protons */
         if (fate == kIHNFtElas  ) {xsec = TMath::Max(0., fXSecPn_Elas -> Evaluate(ke)) *  targZ;
	                            xsec+= TMath::Max(0., fXSecNn_Elas -> Evaluate(ke)) * (targA-targZ);
				    return xsec;}
    else if (fate == kIHNFtInelas) {xsec = TMath::Max(0., fXSecPn_Reac -> Evaluate(ke)) *  targZ;
	                            xsec+= TMath::Max(0., fXSecNn_Reac -> Evaluate(ke)) * (targA-targZ);
				    return xsec;}
    else {
     LOG("INukeData", pWARN) 
        << "Neutrons don't have this fate: " << INukeHadroFates::AsString(fate);
     return 0;
    }
	 /*   }  else if (hpdgc == kPdgGamma) {  
    / * handle gamma * /
         if (fate == kIHNFtInelas) {xsec = TMath::Max(0., fXSecGamp_fs   -> Evaluate(ke)) *  targZ;
	                            xsec+= TMath::Max(0., fXSecGamn_fs   -> Evaluate(ke)) * (targA-targZ);
				    return xsec;}
    else {
     LOG("INukeData", pWARN)
        << "Gamma's don't have this fate: " << INukeHadroFates::AsString(fate);
     return 0;
     }*/
   }
  LOG("INukeData", pWARN) 
      << "Can't handle particles with pdg code = " << hpdgc;

  return 0;
}

double INukeHadroData::Frac(int hpdgc, INukeFateHN_t fate, double ke, int targA, int targZ) const
{
// return the x-section fraction for the input fate for the particle with the 
// input pdg code at the input kinetic energy

  ke = TMath::Max(fMinKinEnergy,   ke);
  ke = TMath::Min(fMaxKinEnergyHN, ke);

  // get x-section
  double xsec = this->XSec(hpdgc,fate,ke,targA,targZ);

  // get max x-section
  double xsec_tot = 0;
       if (hpdgc == kPdgPiP    ){xsec_tot = TMath::Max(0., fXSecPipp_Tot  -> Evaluate(ke)) *  targZ;
				 xsec_tot+= TMath::Max(0., fXSecPipn_Tot  -> Evaluate(ke)) * (targA-targZ);}
  else if (hpdgc == kPdgPiM    ){xsec_tot = TMath::Max(0., fXSecPipn_Tot  -> Evaluate(ke)) *  targZ;
                        	 xsec_tot+= TMath::Max(0., fXSecPipp_Tot  -> Evaluate(ke)) * (targA-targZ);}
  else if (hpdgc == kPdgPi0    ){xsec_tot = TMath::Max(0., fXSecPi0p_Tot  -> Evaluate(ke)) *  targZ;
                        	 xsec_tot+= TMath::Max(0., fXSecPi0n_Tot  -> Evaluate(ke)) * (targA-targZ);}
  else if (hpdgc == kPdgProton ){xsec_tot = TMath::Max(0., fXSecPp_Tot    -> Evaluate(ke)) *  targZ;
                        	 xsec_tot+= TMath::Max(0., fXSecPn_Tot    -> Evaluate(ke)) * (targA-targZ);}
  else if (hpdgc == kPdgNeutron){xsec_tot = TMath::Max(0., fXSecPn_Tot    -> Evaluate(ke)) *  targZ;
                        	 xsec_tot+= TMath::Max(0., fXSecNn_Tot    -> Evaluate(ke)) * (targA-targZ);}
  else if (hpdgc == kPdgGamma  ) xsec_tot = TMath::Max(0., fXSecGamN_Tot  -> Evaluate(ke)); 
  else if (hpdgc == kPdgKP     ) xsec_tot = TMath::Max(0., fXSecKpN_Tot   -> Evaluate(ke));

  // compute fraction
  double frac = (xsec_tot>0) ? xsec/xsec_tot : 0.;
  return frac;
}
//____________________________________________________________________________
double INukeHadroData::IntBounce(const GHepParticle* p, int target, int scode, INukeFateHN_t fate)
{
  // This method returns a random cos(ang) according to a distribution
  // based upon the particle and fate. The sampling uses the 
  // Accept/Reject method, whereby a distribution is bounded above by
  // an envelope, or in this case, a number of envelopes, which can be
  // easily sampled (here, we use uniform distributions). 
  // To get a random value, first the envelope is sampled to 
  // obtain an x-coordinate (cos(ang)), and then another random value
  // is obtained uniformally in the range [0,h(j,0)], where h(j,0)
  // is the height of the j-th envelope. If the point is beneath the
  // distribution, the x-coordinate is accepted, otherwise, we try
  // again.

  RandomGen * rnd = RandomGen::Instance();

  // numEnv is the number of envelopes in the total envelope,
  // that is, the number of seperate simple uniform distributions
  // that will be fit against the distribution in question in the
  // Accept/Reject process of sampling
  int    numEnv = 4;
  int    numPoints = 1000;            // The number of points to be evaluated
                                      //   for the purpose of finding the max
                                      //   value of the distribution
  assert((numPoints%numEnv)==0);      // numPoints/numEnv has to be an integer
  double sr = 2.0 / numEnv;           // Subrange, i.e., range of an envelope
  double cstep = 2.0 / (numPoints);   // Magnitude of the step between eval. points
  
  double ke = (p->E() - p->Mass()) * 1000.0; // ke in MeV
  if (TMath::Abs((int)ke-ke)<.01) ke+=.3;    // make sure ke isn't an integer,
                                             // otherwise sometimes gives weird results
                                             // due to ROOT's Interpolate() function
  double avg = 0.0; // average value in envelop

  // Matrices to hold data; buff holds the distribution
  //   data per envelope from which the max value is 
  //   obtained. That value is then recorded in dist, where
  //   the integral of the envelope to that point is 
  //   also recorded
  
  double * buff = new double[numPoints/numEnv + 1];
  double ** dist = new double*[numEnv];
  for(int ih=0;ih<numEnv;ih++)
  {
    dist[ih] = new double[3];
  }

  // Acc-Rej Sampling Method
  // -- Starting at the beginning of each envelope,
  // this loop evaluates (numPoints) amount of points
  // on the distribution and holds them in buff;
  // then takes the max and  places it in the first row
  // of h. The second row of h contains the interval of
  // the total envelope, up to the current envelope.
  // Thus, when properly normalized, the last value
  // in the second row of h should be 1.
  double totxsec = 0.0;
  for(int i=0;i<numEnv;i++)
    {
      double lbound = -1 + i*sr;

      for(int j=0;j<=numPoints / numEnv; j++)
	{
	  buff[j] = this->XSec(p->Pdg(),target,scode,fate,ke,lbound+j*cstep);
	  avg += buff[j];
	}

      totxsec+=avg;
      avg/= (double(numPoints)/double(numEnv));
      dist[i][0] = TMath::MaxElement(numPoints/numEnv+1,buff);
      dist[i][1] = avg;
      dist[i][2] = dist[i][1] + ((i==0)?0.0:dist[i-1][2]);
      avg=0.0;
    }


  delete [] buff;
  
  int iter=1;         // keep track of iterations
  int env=0;          // envelope index
  double rval = 0.0;  // random value
  double val = 0.0;   // angle value
  
  // Get a random point, see if its in the distribution, and if not
  // then try again.

  rval = rnd->RndFsi().Rndm()*dist[numEnv-1][2];

  env=0;
  // Check which envelope it's in, to 
  // get proper height
  while(env<numEnv)
    {
      if(rval<=dist[env][2]) break;
      else env++;
    }
  if(env==numEnv) env=numEnv - 1;

while(iter)
    {

      // Obtain the correct x-coordinate from the random sample
      val = rnd->RndFsi().Rndm()*sr;
      val += sr*env-1;
      rval = rnd->RndFsi().Rndm()*dist[env][0]; 

      // Test to see if point is in distribution, if it is, stop and return
      if(rval < this->XSec(p->Pdg(),target,scode,fate,ke,val)) break;

      // Possibly an extremely long loop, don't want to
      // hold up the program
      if(iter==1000)
	{
          int NUM_POINTS=2000;
	  int pvalues=0;
	  double points[200]={0};
	  for(int k=0;k<NUM_POINTS;k++)
          {
	    points[int(k/10)]=this->XSec(p->Pdg(),target,scode,fate,ke,-1+(2.0/NUM_POINTS)*k);
	    if(points[int(k/10)]>0) pvalues++;
	  }
          if(pvalues<(.05*NUM_POINTS))
	  {
	    // if it reaches here, one more test...if momenta of particle is
            // extremely low, just give it an angle from a uniform distribution
	    if(p->P4()->P()<.005) // 5 MeV
	    {
              val = 2*rnd->RndFsi().Rndm()-1;
	      break;
            }
            else
            {
 	      LOG("Intranuke", pWARN) << "Hung-up in IntBounce method - Exiting";
	      LOG("Intranuke", pWARN) << (*p);
	      LOG("Intranuke", pWARN) << "Target: " << target << ", Scode: " << scode << ", fate: " << INukeHadroFates::AsString(fate);
	      for(int ie=0;ie<200;ie+=10) {
		LOG("Intranuke", pWARN)   << points[ie+0] << ", " << points[ie+1] << ", " << points[ie+2] << ", "
		   << points[ie+3] << ", " << points[ie+4] << ", " << points[ie+5] << ", " << points[ie+6] << ", "
					   << points[ie+7] << ", " << points[ie+8] << ", " << points[ie+9];
	      }
              for(int ih=0;ih<numEnv;ih++)
              {
                delete [] dist[ih];
              }
	      delete [] dist;

	      return -2.;
            }
          }
	}
      iter++;
    }

  for(int ih=0;ih<numEnv;ih++)
  {
    delete [] dist[ih];
  }
  delete [] dist;

  return val;
}
//___________________________________________________________________________
/*int INukeHadroData::AngleAndProduct(const GHepParticle* p, int target, double &costh, INukeFateHN_t fate)
{
  //
  // Adapted from genie::utils::intranuke::IntBounce()
  // by Aaron Meyer (7/10/09)
  //
  // This method returns a random cos(ang) and a final product
  // PDG number according to a distribution
  // based upon the probe particle and fate. The sampling uses the 
  // Accept/Reject method, described in the notes for IntBounce. 
  // The final state particle is chosen by comparing a random value
  // that was obtained and accepted within the range to one of the
  // particle cross section envelops.
  // If the value is less than that of the first particle
  // cross section, the first particle is chosen; if not, the second
  // particle is chosen.
  //
  // The argument costh should be declared before calling. The cosine
  // value is determined in the function and returned via argument by
  // reference.

  RandomGen * rnd = RandomGen::Instance();

  // numEnv is the number of envelopes in the total envelope,
  // that is, the number of seperate simple uniform distributions
  // that will be fit against the distribution in question in the
  // Accept/Reject process of sampling
  
  int    numEnv = 4;
  int    numPoints = 1000;            // The number of points to be evaluated
                                      //   for the purpose of finding the max
                                      //   value of the distribution
  assert((numPoints%numEnv)==0);      // numPoints/numEnv has to be an integer
  double sr = 2.0 / numEnv;           // Subrange, i.e., range of an envelope
  double cstep = 2.0 / (numPoints);   // Magnitude of the step between eval. points
  int fsPart1 = 0;                    // First possible particle final state
  int fsPart2 = 0;                    // Second possible particle final state
  int endPart = 0;                    // Particle to be returned at completion

    if      ((p->Pdg()==kPdgGamma) && (target==kPdgProton)) {
      fsPart1 = kPdgProton;
      fsPart2 = kPdgNeutron;
  } else if ((p->Pdg()==kPdgGamma) && (target==kPdgNeutron)) {
      fsPart1 = kPdgProton;
      fsPart2 = kPdgNeutron;
  } else {
      LOG("INukeHadroData", pERROR)  << "Cannot handle particle: " << p->Pdg()
				<< " and target: "            << target;
      return endPart;
  }

    LOG("INukeHadroData", pNOTICE)  << "Possible final state particles :" << fsPart1 << " or " <<fsPart2;

  
  double ke = (p->E() - p->Mass()) * 1000.0; // ke in MeV
  double avg = 0.0;                          // Temp value used to find average of envelope

  double * buff = new double[numPoints/numEnv + 1];
  double ** dist = new double*[numEnv];
  for(int ih=0;ih<numEnv;ih++)
  {
    dist[ih] = new double[3];
  }

  // Acc-Rej Sampling Method
  //   see IntBounce above for description
  for(int i=0;i<numEnv;i++)
    {
      double lbound = -1 + i*sr;

      for(int j=0;j<=numPoints / numEnv; j++)
	{
	  buff[j] = this->XSec(p->Pdg(),target,fsPart1,fate,ke,lbound+j*cstep);
	  buff[j]+= this->XSec(p->Pdg(),target,fsPart2,fate,ke,lbound+j*cstep);
	  avg += buff[j];
	}
      
      avg /= (double(numPoints)/double(numEnv));
      dist[i][0] = TMath::MaxElement(numPoints/numEnv+1,buff);
      dist[i][1] = avg;
      dist[i][2] = dist[i][1] + ((i==0)?0.0:dist[i-1][2]);
      avg = 0.0;
      LOG("INukeHadroData", pNOTICE)  << "max xsec value for env "<<i<<": " << dist[i][0];
      LOG("INukeHadroData", pNOTICE)  << "xsec avg for env "<<i<<": " << dist[i][1];
      LOG("INukeHadroData", pNOTICE)  << "xsec avg integral at env "<<i<<": " << dist[i][2];

    }

  delete [] buff;

  int iter=1;         // keep track of iterations
  int env=0;          // envelope index
  double rval = 0.0;  // random value
  double val = 0.0;   // costh value


  // Get a random point, see if its in the distribution, and if not
  // then try again.

  rval = rnd->RndFsi().Rndm()*dist[numEnv-1][2];
  LOG("INukeHadroData", pNOTICE)  << "random variable for envelope decision: " << rval;

  env=0;
  // Check which envelope it's in, to 
  // get proper height
  while(env<numEnv)
    {
      if(rval<=dist[env][2]) break;
      else env++;
    }
  if(env==numEnv) env=numEnv - 1;  
  LOG("INukeHadroData", pNOTICE)  << "value in envelope: " << env;
  LOG("INukeHadroData", pNOTICE)  << "weighted envelope bound: " << dist[env][2];

  while(iter)
    {
      // Obtain the correct x-coordinate from the random sample
      val = rnd->RndFsi().Rndm()*sr; 
      val+= (env*sr - 1);
      rval = rnd->RndFsi().Rndm()*dist[env][0]; 

      // Test to see if point is in distribution, if it is, stop and return
      if(rval < (this->XSec(p->Pdg(),target,fsPart1,fate,ke,val) +
		 this->XSec(p->Pdg(),target,fsPart2,fate,ke,val))) 
      {
	// Determine final state particle
	LOG("INukeHadroData", pNOTICE)  << "particle 1 bound: " << this->XSec(p->Pdg(),target,fsPart1,fate,ke,val);
	LOG("INukeHadroData", pNOTICE)  << "xsec bound: " << (this->XSec(p->Pdg(),target,fsPart1,fate,ke,val) +\
							 this->XSec(p->Pdg(),target,fsPart2,fate,ke,val));

	if (rval < this->XSec(p->Pdg(),target,fsPart1,fate,ke,val))
	  {endPart = fsPart1;
	  LOG("INukeHadroData", pNOTICE)  << "particle 1 taken: " << endPart;
	  }
	else 
	  {endPart = fsPart2;
	  LOG("INukeHadroData", pNOTICE)  << "particle 2 taken: " << endPart;
          }
	costh = val;
	break;
      }

      // Possibly an extremely long loop, don't want to
      // hold up the program
      if(iter==1000)
	{
          int NUM_POINTS=2000;
	  int pvalues=0;
	  double points[200]={0};
	  for(int k=0;k<NUM_POINTS;k++)
          {
	    points[int(k/10)]=this->XSec(p->Pdg(),target,fsPart1,fate,ke,-1+(2.0/NUM_POINTS)*k) +
	                      this->XSec(p->Pdg(),target,fsPart2,fate,ke,-1+(2.0/NUM_POINTS)*k);
	    if(points[int(k/10)]>0) pvalues++;
	  }
          if(pvalues<(.05*NUM_POINTS))
	  {
	    // if it reaches here, one more test...if momenta of particle is
            // extremely low, just give it an angle from a uniform distribution
	    if(p->P4()->P()<.005) // 5 MeV
	    {
              val = 2*rnd->RndFsi().Rndm()-1;
	      break;
            }
            else
            {
 	      LOG("INukeHadroData", pWARN) << "Hung-up in AngleAndProduct method - Exiting";
	      LOG("INukeHadroData", pWARN) << (*p);
	      LOG("INukeHadroData", pWARN) << "Target: " << target << ", Potential final state particles: "
				       << fsPart1 << " " << fsPart2 <<", fate: " << INukeHadroFates::AsString(fate);
	      for(int ie=0;ie<200;ie+=10) {
		LOG("INukeHadroData", pWARN)   << points[ie+0] << ", " << points[ie+1] << ", " << points[ie+2] << ", "
		   << points[ie+3] << ", " << points[ie+4] << ", " << points[ie+5] << ", " << points[ie+6] << ", "
					   << points[ie+7] << ", " << points[ie+8] << ", " << points[ie+9];
	      }

              for(int ih=0;ih<numEnv;ih++)
              {
                delete [] dist[ih];
              }
	      delete [] dist;

	      return 0;
            }
          }
	}
      iter++;
    }
  for(int ih=0;ih<numEnv;ih++)
  {
    delete [] dist[ih];
  }
  delete [] dist;

  //  LOG("INukeHadroData", pNOTICE)  << "return value: " << endPart;

  return endPart;
}*/
//___________________________________________________________________________
