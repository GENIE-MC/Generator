//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>, Rutherford Lab.
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
 @ Jan 9, 2015 - SD, Nick Geary, Tomek Golan
   Added 2014 version of INTRANUKE codes for independent development.
 @ Oct, 2015 - TG
   Added 2015 version of INTRANUKE codes for independent development.
   Include Oset data files.
 @ Apr, 2016 - Flor Blasczyk
   Added K+ cex data files
*/
//____________________________________________________________________________

#include <cassert>
#include <string>

#include <TSystem.h>
#include <TNtupleD.h>
#include <TGraph2D.h>
#include <TTree.h>
#include <TMath.h>
#include <TFile.h>
#include <iostream>

#include "Framework/Conventions/Constants.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Physics/HadronTransport/INukeHadroData2018.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/Numerical/Spline.h"
#include "Framework/ParticleData/PDGCodes.h"

using std::ostringstream;
using std::ios;

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
INukeHadroData2018 * INukeHadroData2018::fInstance = 0;
//____________________________________________________________________________
double INukeHadroData2018::fMinKinEnergy   =    1.0; // MeV
double INukeHadroData2018::fMaxKinEnergyHA =  999.0; // MeV
double INukeHadroData2018::fMaxKinEnergyHN = 1799.0; // MeV
//____________________________________________________________________________
INukeHadroData2018::INukeHadroData2018()
{
  this->LoadCrossSections();
  fInstance = 0;
}
//____________________________________________________________________________
INukeHadroData2018::~INukeHadroData2018()
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

  delete fXSecPp_Cmp; //added to fix memory leak; no noticeable effect, but good convention.
  delete fXSecPn_Cmp;
  delete fXSecNn_Cmp;
  delete fXSecPp_Tot;
  delete fXSecPp_Elas;
  delete fXSecPp_Reac;
  delete fXSecPn_Tot;
  delete fXSecPn_Elas;
  delete fXSecPn_Reac;
  delete fXSecNn_Tot;
  delete fXSecNn_Elas;
  delete fXSecNn_Reac;

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

  // K+N x-section splines
  delete fXSecKpn_Elas;
  delete fXSecKpp_Elas;
  delete fXSecKpn_CEx;
  delete fXSecKpN_Abs;
  delete fXSecKpN_Tot;

  // gamma x-section splines (inelastic only)
  delete fXSecGamp_fs;
  delete fXSecGamn_fs;
  delete fXSecGamN_Tot;

  // N+A x-section splines
  delete fFracPA_Tot;
  //  delete fFracPA_Elas;
  delete fFracPA_Inel;
  delete fFracPA_CEx;
  delete fFracPA_Abs;
  delete fFracPA_PiPro;
  delete fFracNA_Tot;
  //  delete fFracNA_Elas;
  delete fFracNA_Inel;
  delete fFracNA_CEx;
  delete fFracNA_Abs;
  delete fFracNA_PiPro;

  delete fFracPA_Cmp;  //  cmp - add support later
  delete fFracNA_Cmp;

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
  delete fhN2dXSecKpN_CEx;
  
//  delete TPipA_Tot; 
  delete TfracPipA_Abs;
  delete TfracPipA_CEx;
//  delete TfracPipA_Elas;
  delete TfracPipA_Inelas;
  delete TfracPipA_PiPro;
  
  // K+A x-section fraction splines
  delete fFracKA_Tot;
  delete fFracKA_Elas;
  delete fFracKA_CEx;
  delete fFracKA_Inel;
  delete fFracKA_Abs;
  
}
//____________________________________________________________________________
INukeHadroData2018 * INukeHadroData2018::Instance()
{
  if(fInstance == 0) {
    LOG("INukeData", pINFO) << "INukeHadroData2018 late initialization";
    static INukeHadroData2018::Cleaner cleaner;
    cleaner.DummyMethodAndSilentCompiler();
    fInstance = new INukeHadroData2018;
  }
  return fInstance;
}
//____________________________________________________________________________
void INukeHadroData2018::LoadCrossSections(void)
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

  string datafile_NN   = data_dir + "/tot_xsec/intranuke-xsections-NN2014.dat";
  string datafile_pipN = data_dir + "/tot_xsec/intranuke-xsections-pi+N.dat";
  string datafile_pi0N = data_dir + "/tot_xsec/intranuke-xsections-pi0N.dat";
  string datafile_NA   = data_dir + "/tot_xsec/intranuke-fractions-NA2016.dat";
  string datafile_KA   = data_dir + "/tot_xsec/intranuke-fractions-KA.dat";
  string datafile_gamN = data_dir + "/tot_xsec/intranuke-xsections-gamN.dat";
  string datafile_kN   = data_dir + "/tot_xsec/intranuke-xsections-kaonN2018.dat";

  //-- Make sure that all data files are available

  assert( ! gSystem->AccessPathName(datafile_NN.  c_str()) );
  assert( ! gSystem->AccessPathName(datafile_pipN.c_str()) );
  assert( ! gSystem->AccessPathName(datafile_pi0N.c_str()) );
  assert( ! gSystem->AccessPathName(datafile_NA.  c_str()) );
  assert( ! gSystem->AccessPathName(datafile_KA. c_str())  );
  assert( ! gSystem->AccessPathName(datafile_gamN.c_str())  );
  assert( ! gSystem->AccessPathName(datafile_kN.  c_str())  );

  LOG("INukeData", pINFO)  << "Found all necessary data files...";

  //-- Load data files

  TTree data_NN;
  TTree data_pipN;
  TTree data_pi0N;
  TTree data_NA;
  TTree data_KA;
  TTree data_gamN;
  TTree data_kN;

  data_NN.ReadFile(datafile_NN.c_str(),"ke/D:pp_tot/D:pp_elas/D:pp_reac/D:pn_tot/D:pn_elas/D:pn_reac/D:nn_tot/D:nn_elas/D:nn_reac/D:pp_cmp/D:pn_cmp/D:nn_cmp/D");
  data_pipN.ReadFile(datafile_pipN.c_str(),
     "ke/D:pipn_tot/D:pipn_cex/D:pipn_elas/D:pipn_reac/D:pipp_tot/D:pipp_cex/D:pipp_elas/D:pipp_reac/D:pipd_abs");
  data_pi0N.ReadFile(datafile_pi0N.c_str(),
     "ke/D:pi0n_tot/D:pi0n_cex/D:pi0n_elas/D:pi0n_reac/D:pi0p_tot/D:pi0p_cex/D:pi0p_elas/D:pi0p_reac/D:pi0d_abs");
  //data_NA.ReadFile(datafile_NA.c_str(),
  //"ke/D:pA_tot/D:pA_elas/D:pA_inel/D:pA_cex/D:pA_abs/D:pA_pipro/D");  // add support for cmp here (?)
  data_NA.ReadFile(datafile_NA.c_str(),
		   "ke/D:pA_tot/D:pA_inel/D:pA_cex/D:pA_abs/D:pA_pipro/D:pA_cmp/D");  // add support for cmp here (?)
  data_gamN.ReadFile(datafile_gamN.c_str(),
    "ke/D:pi0p_tot/D:pipn_tot/D:pimp_tot/D:pi0n_tot/D:gamp_fs/D:gamn_fs/D:gamN_tot/D");
  data_kN.ReadFile(datafile_kN.c_str(),
		   "ke/D:kpp_elas/D:kpn_elas/D:kpn_cex/D:kp_abs/D:kpN_tot/D");
  data_KA.ReadFile(datafile_KA.c_str(),
     "ke/D:KA_tot/D:KA_elas/D:KA_inel/D:KA_abs/D");

  LOG("INukeData", pDEBUG)  << "Number of data rows in NN : "   << data_NN.GetEntries();
  LOG("INukeData", pDEBUG)  << "Number of data rows in pipN : " << data_pipN.GetEntries();
  LOG("INukeData", pDEBUG)  << "Number of data rows in pi0N : " << data_pi0N.GetEntries();
  LOG("INukeData", pDEBUG)  << "Number of data rows in NA  : "  << data_NA.GetEntries();
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
  fXSecPp_Cmp      = new Spline(&data_NN, "ke:pp_cmp"); //for compound nucleus fate
  fXSecPn_Cmp      = new Spline(&data_NN, "ke:pn_cmp");
  fXSecNn_Cmp      = new Spline(&data_NN, "ke:nn_cmp");

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
  fXSecKpn_CEx    = new Spline(&data_kN,  "ke:kpn_cex");
  fXSecKpN_Abs    = 0; // new Spline(&data_kN,  "ke:kp_abs");  why not used?
  fXSecKpN_Tot    = new Spline(&data_kN,  "ke:kpN_tot");

  // gamma x-section splines
  fXSecGamp_fs     = new Spline(&data_gamN, "ke:gamp_fs");
  fXSecGamn_fs     = new Spline(&data_gamN, "ke:gamn_fs");
  fXSecGamN_Tot    = new Spline(&data_gamN, "ke:gamN_tot");

  // N+A x-section fraction splines
  fFracPA_Tot      = new Spline(&data_NA, "ke:pA_tot");
  //  fFracPA_Elas     = new Spline(&data_NA, "ke:pA_elas");
  fFracPA_Inel     = new Spline(&data_NA, "ke:pA_inel");
  fFracPA_CEx      = new Spline(&data_NA, "ke:pA_cex");
  fFracPA_Abs      = new Spline(&data_NA, "ke:pA_abs");
  fFracPA_PiPro    = new Spline(&data_NA, "ke:pA_pipro");
  fFracNA_Tot      = new Spline(&data_NA, "ke:pA_tot");  // assuming nA same as pA
  //  fFracNA_Elas     = new Spline(&data_NA, "ke:pA_elas");
  fFracNA_Inel     = new Spline(&data_NA, "ke:pA_inel");
  fFracNA_CEx      = new Spline(&data_NA, "ke:pA_cex");
  fFracNA_Abs      = new Spline(&data_NA, "ke:pA_abs");
  fFracNA_PiPro    = new Spline(&data_NA, "ke:pA_pipro");

  fFracPA_Cmp      = new Spline(&data_NA, "ke:pA_cmp");  //cmp - add support later
  fFracNA_Cmp      = new Spline(&data_NA, "ke:pA_cmp");

  // K+A x-section fraction splines
  fFracKA_Tot      = new Spline(&data_KA, "ke:KA_tot");
  fFracKA_Elas     = new Spline(&data_KA, "ke:KA_elas");
  fFracKA_CEx      = 0; // new Spline(&data_KA, "ke:KA_cex"); //Added, but needs to be computed
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

    fhN2dXSecKpN_Elas = new BLI2DNonUnifGrid(hN_kpNelas_nfiles,hN_kpNelas_points_per_file,
			   hN_kpNelas_energies,hN_kpNelas_costh,hN_kpNelas_xsec);
  }
  // kIHNFtCEx, kpn :
  {
    const int hN_kpNcex_nfiles = 18;
    const int hN_kpNcex_points_per_file = 37;
    const int hN_kpNcex_npoints = hN_kpNcex_points_per_file * hN_kpNcex_nfiles;

    double hN_kpNcex_energies[hN_kpNcex_nfiles] = {
     100, 200, 300, 400, 500, 600, 700, 800, 900, 1000,
     1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800
    };

    double hN_kpNcex_costh [hN_kpNcex_points_per_file];
    double hN_kpNcex_xsec  [hN_kpNcex_npoints];

    int ipoint=0;

    for(int ifile = 0; ifile < hN_kpNcex_nfiles; ifile++) {
     // build filename
     ostringstream hN_datafile;
     double ke = hN_kpNcex_energies[ifile];
     hN_datafile << data_dir << "/diff_ang/kpncex/kpcex" << ke << ".txt";
     // read data
     ReadhNFile(
       hN_datafile.str(), ke, hN_kpNcex_points_per_file,
       ipoint, hN_kpNcex_costh, hN_kpNcex_xsec,2);
    }//loop over files

    /*double hN_kpNcex_costh_cond [hN_kpNcex_points_per_file];
    for (int ient = 0; ient < hN_kpNcex_points_per_file; ient++) {
      hN_kpNcex_costh_cond[ient] = hN_kpNcex_costh[ient];
      }*/

    fhN2dXSecKpN_CEx = new BLI2DNonUnifGrid(hN_kpNcex_nfiles,hN_kpNcex_points_per_file,
			   hN_kpNcex_energies,hN_kpNcex_costh,hN_kpNcex_xsec);
  }
//----------------------------------------------------------------------------------------


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

    fhN2dXSecGamPimP_Inelas = new BLI2DNonUnifGrid(hN_gampimpInelas_nfiles,hN_gampimpInelas_points_per_file,
			   hN_gampimpInelas_energies,hN_gampimpInelas_costh,hN_gampimpInelas_xsec);
  }




  TFile TGraphs_file;
  bool saveTGraphsToFile = false;  //true;

  if (saveTGraphsToFile) {
    string filename = "TGraphs.root";
    LOG("INukeHadroData2018", pNOTICE) << "Saving INTRANUKE hadron x-section data to ROOT file: " << filename;
    TGraphs_file.Open(filename.c_str(), "RECREATE");
  }

  /*
  // kIHNFtTot,   pip + A                                            PipA_Tot
  {
    const int pipATot_nfiles = 22;
    const int pipATot_nuclei[pipATot_nfiles] = {1, 2, 3, 4, 6, 7, 9, 12, 16, 27, 28, 32, 40, 48, 56, 58, 63, 93, 120, 165, 181, 209};
    const int pipATot_npoints = 203;

    TPipA_Tot = new TGraph2D(pipATot_npoints);
    TPipA_Tot->SetNameTitle("TPipA_Tot","TPipA_Tot");
    TPipA_Tot->SetDirectory(0);

    int ipoint=0;
    double x, y;

    for(int ifile=0; ifile < pipATot_nfiles; ifile++) {
      ostringstream ADep_datafile;
      int nucleus = pipATot_nuclei[ifile];
      ADep_datafile << data_dir << "/tot_xsec/pipA_tot/pip" << nucleus << "_tot.txt";
      TGraph * buff = new TGraph(ADep_datafile.str().c_str());
      buff->SetNameTitle("buff","buff");
      for(int i=0; i < buff->GetN(); i++) {
	buff -> GetPoint(i,x,y);
	TPipA_Tot -> SetPoint(ipoint,(double)nucleus,x,y);
	ipoint++;
      }
      delete buff;
    }

    if (saveTGraphsToFile) {
      TPipA_Tot -> Write("TPipA_Tot"); // TPipA_Tot will be _key_ name
    }
  }
  */

  // kIHNFtAbs, pip + A                                                            PipA_Abs_frac
  {
    const int pipAAbs_f_nfiles = 18;
    const int pipAAbs_f_nuclei[pipAAbs_f_nfiles] = {1, 2, 3, 4, 7, 9, 12, 16, 27, 48, 56, 58, 63, 93, 120, 165, 181, 209};
    const int pipAAbs_f_npoints = 111;

    TfracPipA_Abs = new TGraph2D(pipAAbs_f_npoints);
    TfracPipA_Abs->SetNameTitle("TfracPipA_Abs","TfracPipA_Abs");
    TfracPipA_Abs->SetDirectory(0);

    int ipoint=0;
    double x, y;
    for(int ifile=0; ifile < pipAAbs_f_nfiles; ifile++) {
      ostringstream ADep_datafile;
      int nucleus = pipAAbs_f_nuclei[ifile];
      ADep_datafile << data_dir << "/tot_xsec/pipA_abs_frac/pip" << nucleus << "_abs_frac.txt";
      TGraph * buff = new TGraph(ADep_datafile.str().c_str());
      buff->SetNameTitle("buff","buff");
      for(int i=0; i < buff->GetN(); i++) {
	buff -> GetPoint(i,x,y);
	TfracPipA_Abs -> SetPoint(ipoint,(double)nucleus,x,y);
	ipoint++;
      }
      delete buff;
    }
    if (saveTGraphsToFile) {
      TfracPipA_Abs -> Write("TfracPipA_Abs");
    }

  }


  // kIHNFtCEx, pip + A      PipA_CEx_frac
  {
    const int pipACEx_f_nfiles = 18;
    const int pipACEx_f_nuclei[pipACEx_f_nfiles] = {1, 2, 3, 4, 7, 9, 12, 16, 27, 48, 56, 58, 63, 93, 120, 165, 181, 209};
    const int pipACEx_f_npoints = 129;

    TfracPipA_CEx = new TGraph2D(pipACEx_f_npoints);
    TfracPipA_CEx->SetNameTitle("TfracPipA_CEx","TfracPipA_CEx");
    TfracPipA_CEx->SetDirectory(0);

    int ipoint=0;
    double x, y;

    for(int ifile=0; ifile < pipACEx_f_nfiles; ifile++) {
      ostringstream ADep_datafile;
      int nucleus = pipACEx_f_nuclei[ifile];
      ADep_datafile << data_dir << "/tot_xsec/pipA_cex_frac/pip" << nucleus << "_cex_frac.txt";
      TGraph * buff = new TGraph(ADep_datafile.str().c_str());
      buff->SetNameTitle("buff","buff");
      for(int i=0; i < buff->GetN(); i++) {
	buff -> GetPoint(i,x,y);
	TfracPipA_CEx -> SetPoint(ipoint,(double)nucleus,x,y);
	ipoint++;
      }
      delete buff;
    }

    if (saveTGraphsToFile) {
      TfracPipA_CEx -> Write("TfracPipA_CEx");
    }

  }



  // kIHNFtCEx, pip + A                                                            PipA_CEx (just for developmental purposes)
  {
    TGraph2D * TPipA_CEx;

    const int pipACEx_nfiles = 18;
    const int pipACEx_nuclei[pipACEx_nfiles] = {1, 2, 3, 4, 7, 9, 12, 16, 27, 48, 56, 58, 63, 93, 120, 165, 181, 209};
    const int pipACEx_npoints = 129;

    TPipA_CEx = new TGraph2D(pipACEx_npoints);
    TPipA_CEx->SetNameTitle("TPipA_CEx","TPipA_CEx");
    TPipA_CEx->SetDirectory(0);

    int ipoint=0;
    double x, y;

    for(int ifile=0; ifile < pipACEx_nfiles; ifile++) {
      ostringstream ADep_datafile;
      int nucleus = pipACEx_nuclei[ifile];
      ADep_datafile << data_dir << "/tot_xsec/pipA_cex/pip" << nucleus << "_cex.txt";
      TGraph * buff = new TGraph(ADep_datafile.str().c_str());
      buff->SetNameTitle("buff","buff");
      for(int i=0; i < buff->GetN(); i++) {
	buff -> GetPoint(i,x,y);
	TPipA_CEx -> SetPoint(ipoint,(double)nucleus,x,y);
	ipoint++;
      }
      delete buff;
    }

    if (saveTGraphsToFile) {
      TPipA_CEx -> Write("TPipA_CEx");
    }
    delete TPipA_CEx;
   
  }

  // kIHNFtAbs, pip + A                                                            PipA_Abs (just for developmental purposes)
  {
    TGraph2D * TPipA_Abs;

    const int pipAAbs_nfiles = 18;
    const int pipAAbs_nuclei[pipAAbs_nfiles] = {1, 2, 3, 4, 7, 9, 12, 16, 27, 48, 56, 58, 63, 93, 120, 165, 181, 209};
    const int pipAAbs_npoints = 111;

    TPipA_Abs = new TGraph2D(pipAAbs_npoints);
    TPipA_Abs->SetNameTitle("TPipA_Abs","TPipA_Abs");
    TPipA_Abs->SetDirectory(0);

    int ipoint=0;
    double x, y;

    for(int ifile=0; ifile < pipAAbs_nfiles; ifile++) {
      ostringstream ADep_datafile;
      int nucleus = pipAAbs_nuclei[ifile];
      ADep_datafile << data_dir << "/tot_xsec/pipA_abs/pip" << nucleus << "_abs.txt";
      TGraph * buff = new TGraph(ADep_datafile.str().c_str());
      buff->SetNameTitle("buff","buff");
      for(int i=0; i < buff->GetN(); i++) {
	buff -> GetPoint(i,x,y);
	TPipA_Abs -> SetPoint(ipoint,(double)nucleus,x,y);
	ipoint++;
      }
      delete buff;
    }

    if (saveTGraphsToFile) {
      TPipA_Abs -> Write("TPipA_Abs");
    }
    delete TPipA_Abs;
   
  }

  // kIHNFtElas, pip + A                                                            PipA_Elas (just for developmental purposes)
  {
    TGraph2D * TPipA_Elas;

    const int pipAElas_nfiles = 18;
    const int pipAElas_nuclei[pipAElas_nfiles] = {1, 2, 3, 4, 7, 9, 12, 16, 27, 48, 56, 58, 63, 93, 120, 165, 181, 209};
    const int pipAElas_npoints = 125;

    TPipA_Elas = new TGraph2D(pipAElas_npoints);
    TPipA_Elas->SetNameTitle("TPipA_Elas","TPipA_Elas");
    TPipA_Elas->SetDirectory(0);

    int ipoint=0;
    double x, y;

    for(int ifile=0; ifile < pipAElas_nfiles; ifile++) {
      ostringstream ADep_datafile;
      int nucleus = pipAElas_nuclei[ifile];
      ADep_datafile << data_dir << "/tot_xsec/pipA_elas/pip" << nucleus << "_elas.txt";
      TGraph * buff = new TGraph(ADep_datafile.str().c_str());
      buff->SetNameTitle("buff","buff");
      for(int i=0; i < buff->GetN(); i++) {
	buff -> GetPoint(i,x,y);
	TPipA_Elas -> SetPoint(ipoint,(double)nucleus,x,y);
	ipoint++;
      }
      delete buff;
    }

    if (saveTGraphsToFile) {
      TPipA_Elas -> Write("TPipA_Elas");
    }
    delete TPipA_Elas; 
  }

  // kIHNFtInelas, pip + A                                                            PipA_Inelas (just for developmental purposes)
  {
    TGraph2D * TPipA_Inelas;

    const int pipAInelas_nfiles = 20;
    const int pipAInelas_nuclei[pipAInelas_nfiles] = {1, 2, 3, 4, 7, 9, 12, 16, 27, 40, 48, 56, 58, 63, 93, 120, 165, 181, 208, 209};
    const int pipAInelas_npoints = 118;

    TPipA_Inelas = new TGraph2D(pipAInelas_npoints);
    TPipA_Inelas->SetNameTitle("TPipA_Inelas","TPipA_Inelas");
    TPipA_Inelas->SetDirectory(0);

    int ipoint=0;
    double x, y;

    for(int ifile=0; ifile < pipAInelas_nfiles; ifile++) {
      ostringstream ADep_datafile;
      int nucleus = pipAInelas_nuclei[ifile];
      ADep_datafile << data_dir << "/tot_xsec/pipA_inelas/pip" << nucleus << "_inelas.txt";
      TGraph * buff = new TGraph(ADep_datafile.str().c_str());
      buff->SetNameTitle("buff","buff");
      for(int i=0; i < buff->GetN(); i++) {
	buff -> GetPoint(i,x,y);
	TPipA_Inelas -> SetPoint(ipoint,(double)nucleus,x,y);
	ipoint++;
      }
      delete buff;
    }

    if (saveTGraphsToFile) {
      TPipA_Inelas -> Write("TPipA_Inelas");
    }
    delete TPipA_Inelas;
  }


  /*
  // kIHNFtElas, pip + A                                                            PipA_Elas_frac
  {
    const int pipAElas_f_nfiles = 18;
    const int pipAElas_f_nuclei[pipAElas_f_nfiles] = {1, 2, 3, 4, 7, 9, 12, 16, 27, 48, 56, 58, 63, 93, 120, 165, 181, 209};
    const int pipAElas_f_npoints = 125;

    TfracPipA_Elas = new TGraph2D(pipAElas_f_npoints);
    TfracPipA_Elas->SetNameTitle("TfracPipA_Elas","TfracPipA_Elas");
    TfracPipA_Elas->SetDirectory(0);

    int ipoint=0;
    double x, y;

    for(int ifile=0; ifile < pipAElas_f_nfiles; ifile++) {
      ostringstream ADep_datafile;
      int nucleus = pipAElas_f_nuclei[ifile];
      ADep_datafile << data_dir << "/tot_xsec/pipA_elas_frac/pip" << nucleus << "_elas_frac.txt";
      TGraph * buff = new TGraph(ADep_datafile.str().c_str());
      buff->SetNameTitle("buff","buff");
      for(int i=0; i < buff->GetN(); i++) {
	buff -> GetPoint(i,x,y);
	TfracPipA_Elas -> SetPoint(ipoint,(double)nucleus,x,y);
	ipoint++;
      }
      delete buff;
    }

    if (saveTGraphsToFile) {
      TfracPipA_Elas -> Write("TfracPipA_Elas");
    }

  }
  */

  // kIHNFtInelas, pip + A                                                            PipA_Inelas_frac
  {
    const int pipAInelas_f_nfiles = 20;
    const int pipAInelas_f_nuclei[pipAInelas_f_nfiles] = {1, 2, 3, 4, 7, 9, 12, 16, 27, 40, 48, 56, 58, 63, 93, 120, 165, 181, 208, 209};
    const int pipAInelas_f_npoints = 118;

    TfracPipA_Inelas = new TGraph2D(pipAInelas_f_npoints);
    TfracPipA_Inelas->SetNameTitle("TfracPipA_Inelas","TfracPipA_Inelas");
    TfracPipA_Inelas->SetDirectory(0);

    int ipoint=0;
    double x, y;

    for(int ifile=0; ifile < pipAInelas_f_nfiles; ifile++) {
      ostringstream ADep_datafile;
      int nucleus = pipAInelas_f_nuclei[ifile];
      ADep_datafile << data_dir << "/tot_xsec/pipA_inelas_frac/pip" << nucleus << "_inelas_frac.txt";
      TGraph * buff = new TGraph(ADep_datafile.str().c_str());
      buff->SetNameTitle("buff","buff");
      for(int i=0; i < buff->GetN(); i++) {
	buff -> GetPoint(i,x,y);
	TfracPipA_Inelas -> SetPoint(ipoint,(double)nucleus,x,y);
	ipoint++;
      }
      delete buff;
    }

    if (saveTGraphsToFile) {
      TfracPipA_Inelas -> Write("TfracPipA_Inelas");
    }

  }


  // kIHNFtPiPro, pip + A                                                            PipA_PiPro_frac
   {
    const int pipAPiPro_f_nfiles = 17;
    const int pipAPiPro_f_nuclei[pipAPiPro_f_nfiles] = {1, 2, 3, 4, 7, 9, 12, 16, 48, 56, 58, 63, 93, 120, 165, 181, 209};
    const int pipAPiPro_f_npoints = 76;

    TfracPipA_PiPro = new TGraph2D(pipAPiPro_f_npoints);
    TfracPipA_PiPro->SetNameTitle("TfracPipA_PiPro","TfracPipA_PiPro");
    TfracPipA_PiPro->SetDirectory(0);

    int ipoint=0;
    double x, y;

    for(int ifile=0; ifile < pipAPiPro_f_nfiles; ifile++) {
      ostringstream ADep_datafile;
      int nucleus = pipAPiPro_f_nuclei[ifile];
      ADep_datafile << data_dir << "/tot_xsec/pipA_pipro_frac/pip" << nucleus << "_pipro_frac.txt";
      TGraph * buff = new TGraph(ADep_datafile.str().c_str());
      buff->SetNameTitle("buff","buff");
      for(int i=0; i < buff->GetN(); i++) {
	buff -> GetPoint(i,x,y);
	TfracPipA_PiPro -> SetPoint(ipoint,(double)nucleus,x,y);
	ipoint++;
      }
      delete buff;
    }

    if (saveTGraphsToFile) {
      TfracPipA_PiPro -> Write("TfracPipA_PiPro");
    }
   }

   TGraphs_file.Close();

   LOG("INukeData", pINFO)  << "Done building x-section splines...";
   
}
//____________________________________________________________________________
void INukeHadroData2018::ReadhNFile(
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
double INukeHadroData2018::XSec(
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
    else if(hpdgc == kPdgKP && tgtpdgc == kPdgNeutron) {
    	ke_eval = TMath::Min(ke_eval, 1799.);
    	ke_eval = TMath::Max(ke_eval,  100.);
    	return fhN2dXSecKpN_CEx->Evaluate(ke_eval, costh_eval);
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
double INukeHadroData2018::FracADep(int hpdgc, INukeFateHA_t fate, double ke, int targA) const
{
  // return the x-section fraction for the input fate for the particle with the input pdg
  // code and the target with the input mass number at the input kinetic energy

  ke = TMath::Max(fMinKinEnergy,   ke);  // ke >= 1 MeV
  ke = TMath::Min(fMaxKinEnergyHA, ke);  // ke <= 999 MeV

  targA = TMath::Min(208, targA);  // A <= 208

  LOG("INukeData", pDEBUG)  << "Querying hA cross section at ke  = " << ke << " and target " << targA;

  // Handle pions (currently the same cross sections are used for pi+, pi-, and pi0)
  if ( hpdgc == kPdgPiP || hpdgc == kPdgPiM || hpdgc == kPdgPi0 ) {

    double frac_cex = TfracPipA_CEx->Interpolate(targA, ke);
    //double frac_elas = TfracPipA_Elas->Interpolate(targA, ke);
    double frac_inelas = TfracPipA_Inelas->Interpolate(targA, ke);
    double frac_abs = TfracPipA_Abs->Interpolate(targA, ke);
    double frac_pipro = TfracPipA_PiPro->Interpolate(targA, ke);

    // Protect against unitarity violation due to interpolation problems
    // by renormalizing all available fate fractions to unity.
    double total = frac_cex + frac_inelas + frac_abs + frac_pipro; // + frac_elas

    if ( fate == kIHAFtCEx ) return frac_cex / total;
  //else if ( fate == kIHAFtElas   ) return frac_elas / total;
    else if ( fate == kIHAFtInelas ) return frac_inelas / total;
    else if ( fate == kIHAFtAbs    ) return frac_abs / total;
    else if ( fate == kIHAFtPiProd ) return frac_pipro / total;
    else {
      std::string sign("+");
      if ( hpdgc == kPdgPiM ) sign = "-";
      else if ( hpdgc == kPdgPi0 ) sign = "0";
      LOG("INukeData", pWARN) << "Pi" << sign << "'s don't have this fate: " << INukeHadroFates::AsString(fate);
      return 0.;
    }
  }

  LOG("INukeData", pWARN) << "Can't handle particles with pdg code = " << hpdgc;
  return 0.;
}
//____________________________________________________________________________
double INukeHadroData2018::FracAIndep(int hpdgc, INukeFateHA_t fate, double ke) const
{
  // return the x-section fraction for the input fate for the particle with the input pdg
  // code at the input kinetic energy
  ke = TMath::Max(fMinKinEnergy,   ke);
  ke = TMath::Min(fMaxKinEnergyHA, ke);

  LOG("INukeData", pDEBUG)  << "Querying hA cross section at ke = " << ke;

  // TODO: reduce code duplication here
  if (hpdgc == kPdgProton) {
    // handle protons
    double frac_cex = fFracPA_CEx->Evaluate(ke);
    double frac_inelas = fFracPA_Inel->Evaluate(ke);
    double frac_abs = fFracPA_Abs->Evaluate(ke);
    double frac_pipro = fFracPA_PiPro->Evaluate(ke);
    double frac_comp = fFracPA_Cmp->Evaluate(ke);

    // Protect against unitarity violation due to interpolation problems
    // by renormalizing all available fate fractions to unity.
    double total = frac_cex + frac_inelas + frac_abs + frac_pipro + frac_comp; // + frac_elas

    if ( fate == kIHAFtCEx ) return frac_cex / total;
  //else if ( fate == kIHAFtElas   ) return frac_elas / total;
    else if ( fate == kIHAFtInelas ) return frac_inelas / total;
    else if ( fate == kIHAFtAbs    ) return frac_abs / total;
    else if ( fate == kIHAFtPiProd ) return frac_pipro / total;
    else if ( fate == kIHAFtCmp    ) return frac_comp / total; // cmp - add support for this later
    else {
      LOG("INukeData", pWARN)
        << "Protons don't have this fate: " << INukeHadroFates::AsString(fate);
      return 0;
    }
  }
  else if (hpdgc == kPdgNeutron) {
    // handle neutrons
    double frac_cex = fFracNA_CEx->Evaluate(ke);
    double frac_inelas = fFracNA_Inel->Evaluate(ke);
    double frac_abs = fFracNA_Abs->Evaluate(ke);
    double frac_pipro = fFracNA_PiPro->Evaluate(ke);
    double frac_comp = fFracNA_Cmp->Evaluate(ke);

    // Protect against unitarity violation due to interpolation problems
    // by renormalizing all available fate fractions to unity.
    double total = frac_cex + frac_inelas + frac_abs + frac_pipro + frac_comp; // + frac_elas

    if ( fate == kIHAFtCEx ) return frac_cex / total;
  //else if ( fate == kIHAFtElas   ) return frac_elas / total;
    else if ( fate == kIHAFtInelas ) return frac_inelas / total;
    else if ( fate == kIHAFtAbs    ) return frac_abs / total;
    else if ( fate == kIHAFtPiProd ) return frac_pipro / total;
    else if ( fate == kIHAFtCmp    ) return frac_comp / total; // cmp - add support for this later
    else {
      LOG("INukeData", pWARN)
        << "Neutrons don't have this fate: " << INukeHadroFates::AsString(fate);
      return 0;
    }
  }
  else if (hpdgc == kPdgKP) {
    // handle K+
    double frac_inelas = fFracKA_Inel->Evaluate(ke);
  //double frac_elas = fFracKA_Elas->Evaluate(ke);
    double frac_abs = fFracKA_Abs->Evaluate(ke);

    // Protect against unitarity violation due to interpolation problems
    // by renormalizing all available fate fractions to unity.
    double total = frac_inelas + frac_abs; // + frac_elas

    if ( fate == kIHAFtInelas ) return frac_inelas / total;
    else if ( fate == kIHAFtAbs ) return frac_abs / total;
    else {
      LOG("INukeData", pWARN)
        << "K+'s don't have this fate: " << INukeHadroFates::AsString(fate);
      return 0.;
    }
  }
  LOG("INukeData", pWARN) << "Can't handle particles with pdg code = " << hpdgc;
  return 0.;
}
//____________________________________________________________________________
double INukeHadroData2018::XSec(int hpdgc, INukeFateHN_t fate, double ke, int targA, int targZ) const
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
    else if (fate == kIHNFtCmp) {xsec = TMath::Max(0., fXSecPp_Cmp -> Evaluate(ke)) *  targZ;
                                    xsec+= TMath::Max(0., fXSecPn_Cmp -> Evaluate(ke)) * (targA-targZ);
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
    else if (fate == kIHNFtCmp) {xsec = TMath::Max(0., fXSecPp_Cmp -> Evaluate(ke)) *  targZ;
                                    xsec+= TMath::Max(0., fXSecPn_Cmp -> Evaluate(ke)) * (targA-targZ);
                                    return xsec;}
    else {
     LOG("INukeData", pWARN)
        << "Neutrons don't have this fate: " << INukeHadroFates::AsString(fate);
     return 0;
    }
    //Adding here kaons, why elastic only on protons? hA or hN? No _Reac for kaons...
    } else if (hpdgc == kPdgKP) {
    /* handle K+ */
    	if (fate == kIHNFtCEx   ) {xsec = TMath::Max(0., fXSecKpn_CEx  -> Evaluate(ke)) *  targZ;
	                            xsec+= TMath::Max(0., fXSecKpn_CEx  -> Evaluate(ke)) * (targA-targZ);
				    return xsec;}
    	else if (fate == kIHNFtElas  ) {xsec = TMath::Max(0., fXSecKpn_Elas -> Evaluate(ke)) *  targZ;
	                            xsec+= TMath::Max(0., fXSecKpn_Elas -> Evaluate(ke)) * (targA-targZ);
				    return xsec;}
    	/*else if (fate == kIHNFtAbs   ) {xsec = TMath::Max(0., fXSecKpd_Abs  -> Evaluate(ke)) *  targA;
				    return xsec;}*/
    	else {
    		LOG("INukeData", pWARN)
        	<< "K+'s don't have this fate: " << INukeHadroFates::AsString(fate);
     	return 0;
    }
    //------------------------------------------------
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

double INukeHadroData2018::Frac(int hpdgc, INukeFateHN_t fate, double ke, int targA, int targZ) const
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
double INukeHadroData2018::IntBounce(const GHepParticle* p, int target, int scode, INukeFateHN_t fate)
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
