//____________________________________________________________________________
/*
 Copyright (C) 2003-2006, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>, Rutherford Lab.
         Steve Dytman <dytman+@pitt.edu>, Pittsburgh Univ.
         February 01, 2007

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
#include "PDG/PDGCodes.h"

using std::cout;
using std::endl;

using namespace genie;

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
  cout << "INukeHadroData singleton dtor: "
                      << "Deleting all hadron cross section splines" << endl;

  // N+N x-section splines
  delete fXSecPN_Tot;
  delete fXSecPN_Elas;
  delete fXSecPN_Reac; 
  delete fXSecNN_Tot;     
  delete fXSecNN_Elas;   
  delete fXSecNN_Reac; 

  // pi+N x-section splines
  delete fXSecPipN_Tot;
  delete fXSecPipN_CEx;
  delete fXSecPipN_Elas;
  delete fXSecPipN_Reac;
  delete fXSecPipN_Abs;
  delete fXSecPimN_Tot;
  delete fXSecPimN_CEx;
  delete fXSecPimN_Elas;
  delete fXSecPimN_Reac;
  delete fXSecPimN_Abs;
  delete fXSecPi0N_Tot;
  delete fXSecPi0N_CEx;
  delete fXSecPi0N_Elas;
  delete fXSecPi0N_Reac;
  delete fXSecPi0N_Abs;

  // N+A x-section splines
  delete fXSecPA_Tot;
  delete fXSecPA_Elas;
  delete fXSecPA_Inel;
  delete fXSecPA_CEx;
  delete fXSecPA_NP;
  delete fXSecPA_PP;
  delete fXSecPA_NPP;
  delete fXSecPA_NNP;
  delete fXSecPA_NNPPP;
  delete fXSecPA_NPip;
  delete fXSecPA_NPipPi0;
  delete fXSecNA_Tot;
  delete fXSecNA_Elas;
  delete fXSecNA_Inel;
  delete fXSecNA_CEx;
  delete fXSecNA_NP;
  delete fXSecNA_PP;
  delete fXSecNA_NPP;
  delete fXSecNA_NNP;
  delete fXSecNA_NNPPP;
  delete fXSecNA_NPip;
  delete fXSecNA_NPipPi0;

  // pi+A x-section splines
  delete fXSecPipA_Tot;
  delete fXSecPipA_Elas;
  delete fXSecPipA_Inel;
  delete fXSecPipA_CEx;
  delete fXSecPipA_NP;
  delete fXSecPipA_PP;
  delete fXSecPipA_NPP;
  delete fXSecPipA_NNP;
  delete fXSecPipA_NNPP;
  delete fXSecPipA_NPipPi0;
  delete fXSecPimA_Tot;
  delete fXSecPimA_Elas;
  delete fXSecPimA_Inel;
  delete fXSecPimA_CEx;
  delete fXSecPimA_NP;
  delete fXSecPimA_PP;
  delete fXSecPimA_NPP;
  delete fXSecPimA_NNP;
  delete fXSecPimA_NNPP;
  delete fXSecPimA_NPipPi0;
  delete fXSecPi0A_Tot;
  delete fXSecPi0A_Elas;
  delete fXSecPi0A_Inel;
  delete fXSecPi0A_CEx;
  delete fXSecPi0A_NP;
  delete fXSecPi0A_PP;
  delete fXSecPi0A_NPP;
  delete fXSecPi0A_NNP;
  delete fXSecPi0A_NNPP;
  delete fXSecPi0A_NPipPi0;
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

  //-- Get the directory with the SAID hadron cross section data (search for
  //   $GINUKEHADRONDATA or use default: $GENIE/data/hadron_xsec)
  string data_dir = (gSystem->Getenv("GINUKEHADRONDATA")) ?
             string(gSystem->Getenv("GINUKEHADRONDATA")) :
             string(gSystem->Getenv("GENIE")) + string("/data/hadron_xsec");

  LOG("INukeData", pNOTICE)  
      << "Loading INTRANUKE hadron data from: " << data_dir;

  //-- Build filenames

  string datafile_NN  = data_dir + "/intranuke-xsections-NN.dat";
  string datafile_piN = data_dir + "/intranuke-xsections-piN.dat";
  string datafile_NA  = data_dir + "/intranuke-xsections-NA.dat";
  string datafile_piA = data_dir + "/intranuke-xsections-piA.dat";

  //-- Make sure that all data files are available

  assert( ! gSystem->AccessPathName(datafile_NN. c_str()) );
  assert( ! gSystem->AccessPathName(datafile_piN.c_str()) );
  assert( ! gSystem->AccessPathName(datafile_NA. c_str()) );
  assert( ! gSystem->AccessPathName(datafile_piA.c_str()) );

  LOG("INukeData", pNOTICE)  << "Found all necessary data files...";

  //-- Load data files

  TTree data_NN;
  TTree data_piN;
  TTree data_NA;
  TTree data_piA;

  data_NN.ReadFile(datafile_NN.c_str(),
     "ke/D:pN_tot/D:pN_elas/D:pN_reac/D:nN_tot/D:nN_elas/D:nN_reac/D");
  data_piN.ReadFile(datafile_piN.c_str(),
     "ke/D:piN_tot/D:piN_cex/D:piN_elas/D:piN_reac/D:piN_abs/D:pi0N_tot/D:pi0N_cex/D:pi0N_elas/D:pi0N_reac/D:pi0N_abs/D");
  data_NA.ReadFile(datafile_NA.c_str(),
     "ke/D:pA_tot/D:pA_elas/D:pA_inel/D:pA_cex/D:pA_np/D:pA_pp/D:pA_npp/D:pA_nnp/D:pA_2n3p/D:pA_npip/D:pA_npippi0/D");
  data_piA.ReadFile(datafile_piA.c_str(),
     "ke/D:piA_tot/D:piA_elas/D:piA_inel/D:piA_cex/D:piA_np/D:piA_pp/D:piA_npp/D:piA_nnp/D:piA_2n2p/D:piA_npippi0/D");

  LOG("INukeData", pDEBUG)  << "Number of data rows in NN  : " << data_NN.GetEntries();
  LOG("INukeData", pDEBUG)  << "Number of data rows in piN : " << data_piN.GetEntries();
  LOG("INukeData", pDEBUG)  << "Number of data rows in NA  : " << data_NA.GetEntries();
  LOG("INukeData", pDEBUG)  << "Number of data rows in piA : " << data_piA.GetEntries();

  LOG("INukeData", pNOTICE)  << "Done loading all x-section files...";

  //-- Build x-section splines

  // N+N x-section splines
  fXSecPN_Tot      = new Spline(&data_NN, "ke:pN_tot");     
  fXSecPN_Elas     = new Spline(&data_NN, "ke:pN_elas");      
  fXSecPN_Reac     = new Spline(&data_NN, "ke:pN_reac");      
  fXSecNN_Tot      = new Spline(&data_NN, "ke:nN_tot");     
  fXSecNN_Elas     = new Spline(&data_NN, "ke:nN_elas");      
  fXSecNN_Reac     = new Spline(&data_NN, "ke:nN_reac");      

  // pi+N x-section splines
  fXSecPipN_Tot    = new Spline(&data_piN, "ke:piN_tot");
  fXSecPipN_CEx    = new Spline(&data_piN, "ke:piN_cex");
  fXSecPipN_Elas   = new Spline(&data_piN, "ke:piN_elas");
  fXSecPipN_Reac   = new Spline(&data_piN, "ke:piN_reac");
  fXSecPipN_Abs    = new Spline(&data_piN, "ke:piN_abs");
  fXSecPimN_Tot    = new Spline(&data_piN, "ke:piN_tot");  // same as for pi+
  fXSecPimN_CEx    = new Spline(&data_piN, "ke:piN_cex");  
  fXSecPimN_Elas   = new Spline(&data_piN, "ke:piN_elas"); 
  fXSecPimN_Reac   = new Spline(&data_piN, "ke:piN_reac"); 
  fXSecPimN_Abs    = new Spline(&data_piN, "ke:piN_abs");  
  fXSecPi0N_Tot    = new Spline(&data_piN, "ke:pi0N_tot");
  fXSecPi0N_CEx    = new Spline(&data_piN, "ke:pi0N_cex");
  fXSecPi0N_Elas   = new Spline(&data_piN, "ke:pi0N_elas");
  fXSecPi0N_Reac   = new Spline(&data_piN, "ke:pi0N_reac");
  fXSecPi0N_Abs    = new Spline(&data_piN, "ke:pi0N_abs");

  // N+A x-section splines
  fXSecPA_Tot      = new Spline(&data_NA, "ke:pA_tot");
  fXSecPA_Elas     = new Spline(&data_NA, "ke:pA_elas");
  fXSecPA_Inel     = new Spline(&data_NA, "ke:pA_inel");   
  fXSecPA_CEx      = new Spline(&data_NA, "ke:pA_cex");   
  fXSecPA_NP       = new Spline(&data_NA, "ke:pA_np");  
  fXSecPA_PP       = new Spline(&data_NA, "ke:pA_pp");  
  fXSecPA_NPP      = new Spline(&data_NA, "ke:pA_npp");  
  fXSecPA_NNP      = new Spline(&data_NA, "ke:pA_nnp");  
  fXSecPA_NNPPP    = new Spline(&data_NA, "ke:pA_2n3p");  
  fXSecPA_NPip     = new Spline(&data_NA, "ke:pA_npip");  
  fXSecPA_NPipPi0  = new Spline(&data_NA, "ke:pA_npippi0");  
  fXSecNA_Tot      = new Spline(&data_NA, "ke:pA_tot");  // assuming nA same as pA
  fXSecNA_Elas     = new Spline(&data_NA, "ke:pA_elas"); 
  fXSecNA_Inel     = new Spline(&data_NA, "ke:pA_inel");   
  fXSecNA_CEx      = new Spline(&data_NA, "ke:pA_cex");   
  fXSecNA_NP       = new Spline(&data_NA, "ke:pA_np");  
  fXSecNA_PP       = new Spline(&data_NA, "ke:pA_pp");  
  fXSecNA_NPP      = new Spline(&data_NA, "ke:pA_npp");  
  fXSecNA_NNP      = new Spline(&data_NA, "ke:pA_nnp");  
  fXSecNA_NNPPP    = new Spline(&data_NA, "ke:pA_2n3p");  
  fXSecNA_NPip     = new Spline(&data_NA, "ke:pA_npip");  
  fXSecNA_NPipPi0  = new Spline(&data_NA, "ke:pA_npippi0");  

  // pi+A x-section splines
  fXSecPipA_Tot     = new Spline(&data_piA, "ke:piA_tot");    
  fXSecPipA_Elas    = new Spline(&data_piA, "ke:piA_elas");    
  fXSecPipA_Inel    = new Spline(&data_piA, "ke:piA_inel");    
  fXSecPipA_CEx     = new Spline(&data_piA, "ke:piA_cex");    
  fXSecPipA_NP      = new Spline(&data_piA, "ke:piA_np");    
  fXSecPipA_PP      = new Spline(&data_piA, "ke:piA_pp");    
  fXSecPipA_NPP     = new Spline(&data_piA, "ke:piA_npp");    
  fXSecPipA_NNP     = new Spline(&data_piA, "ke:piA_nnp");    
  fXSecPipA_NNPP    = new Spline(&data_piA, "ke:piA_2n2p");    
  fXSecPipA_NPipPi0 = new Spline(&data_piA, "ke:piA_npippi0");    
  fXSecPimA_Tot     = new Spline(&data_piA, "ke:piA_tot");    
  fXSecPimA_Elas    = new Spline(&data_piA, "ke:piA_elas");    
  fXSecPimA_Inel    = new Spline(&data_piA, "ke:piA_inel");    
  fXSecPimA_CEx     = new Spline(&data_piA, "ke:piA_cex");    
  fXSecPimA_NP      = new Spline(&data_piA, "ke:piA_np");    
  fXSecPimA_PP      = new Spline(&data_piA, "ke:piA_pp");    
  fXSecPimA_NPP     = new Spline(&data_piA, "ke:piA_npp");    
  fXSecPimA_NNP     = new Spline(&data_piA, "ke:piA_nnp");    
  fXSecPimA_NNPP    = new Spline(&data_piA, "ke:piA_2n2p");    
  fXSecPimA_NPipPi0 = new Spline(&data_piA, "ke:piA_npippi0");    
  fXSecPi0A_Tot     = new Spline(&data_piA, "ke:piA_tot");    
  fXSecPi0A_Elas    = new Spline(&data_piA, "ke:piA_elas");    
  fXSecPi0A_Inel    = new Spline(&data_piA, "ke:piA_inel");    
  fXSecPi0A_CEx     = new Spline(&data_piA, "ke:piA_cex");    
  fXSecPi0A_NP      = new Spline(&data_piA, "ke:piA_np");    
  fXSecPi0A_PP      = new Spline(&data_piA, "ke:piA_pp");    
  fXSecPi0A_NPP     = new Spline(&data_piA, "ke:piA_npp");    
  fXSecPi0A_NNP     = new Spline(&data_piA, "ke:piA_nnp");    
  fXSecPi0A_NNPP    = new Spline(&data_piA, "ke:piA_2n2p");    
  fXSecPi0A_NPipPi0 = new Spline(&data_piA, "ke:piA_npippi0");    

  LOG("INukeData", pNOTICE)  << "Done building x-section splines...";
}
//____________________________________________________________________________
double INukeHadroData::XSec(int hpdgc, INukeFateHA_t fate, double ke) const
{
// return the x-section for the input fate for the particle with the input pdg 
// code at the input kinetic energy
//
  ke = TMath::Max(fMinKinEnergy,   ke);
  ke = TMath::Min(fMaxKinEnergyHA, ke);

  if(hpdgc == kPdgProton) {
   /* handle protons */
        if (fate == kIHAFtCEx    ) return TMath::Max(0., fXSecPA_CEx     -> Evaluate (ke));
   else if (fate == kIHAFtElas   ) return TMath::Max(0., fXSecPA_Elas    -> Evaluate (ke));
   else if (fate == kIHAFtInelas ) return TMath::Max(0., fXSecPA_Inel    -> Evaluate (ke));
   else if (fate == kIHAFtAbsNP  ) return TMath::Max(0., fXSecPA_NP      -> Evaluate (ke));
   else if (fate == kIHAFtAbsPP  ) return TMath::Max(0., fXSecPA_PP      -> Evaluate (ke));
   else if (fate == kIHAFtAbsNPP ) return TMath::Max(0., fXSecPA_NPP     -> Evaluate (ke));
   else if (fate == kIHAFtAbsNNP ) return TMath::Max(0., fXSecPA_NNP     -> Evaluate (ke));
   else if (fate == kIHAFtAbs2N3P) return TMath::Max(0., fXSecPA_NNPPP   -> Evaluate (ke));
   else if (fate == kIHAFtNPip   ) return TMath::Max(0., fXSecPA_NPip    -> Evaluate (ke));
   else if (fate == kIHAFtNPipPi0) return TMath::Max(0., fXSecPA_NPipPi0 -> Evaluate (ke));
   else {
     LOG("INukeData", pWARN) 
       << "Protons don't have this fate: " << INukeHadroFates::AsString(fate);
     return 0;
   }

  } else if (hpdgc == kPdgNeutron) {
   /* handle neutrons */
        if (fate == kIHAFtCEx    ) return TMath::Max(0., fXSecNA_CEx     -> Evaluate (ke));
   else if (fate == kIHAFtElas   ) return TMath::Max(0., fXSecNA_Elas    -> Evaluate (ke));
   else if (fate == kIHAFtInelas ) return TMath::Max(0., fXSecNA_Inel    -> Evaluate (ke));
   else if (fate == kIHAFtAbsNP  ) return TMath::Max(0., fXSecNA_NP      -> Evaluate (ke));
   else if (fate == kIHAFtAbsPP  ) return TMath::Max(0., fXSecNA_PP      -> Evaluate (ke));
   else if (fate == kIHAFtAbsNPP ) return TMath::Max(0., fXSecNA_NPP     -> Evaluate (ke));
   else if (fate == kIHAFtAbsNNP ) return TMath::Max(0., fXSecNA_NNP     -> Evaluate (ke));
   else if (fate == kIHAFtAbs2N3P) return TMath::Max(0., fXSecNA_NNPPP   -> Evaluate (ke));
   else if (fate == kIHAFtNPip   ) return TMath::Max(0., fXSecNA_NPip    -> Evaluate (ke));
   else if (fate == kIHAFtNPipPi0) return TMath::Max(0., fXSecNA_NPipPi0 -> Evaluate (ke));
   else {
     LOG("INukeData", pWARN) 
       << "Neutrons don't have this fate: " << INukeHadroFates::AsString(fate);
     return 0;
   }

  } else if (hpdgc == kPdgPiP) {
   /* handle pi+ */
        if (fate == kIHAFtCEx    ) return TMath::Max(0., fXSecPipA_CEx     -> Evaluate (ke));
   else if (fate == kIHAFtElas   ) return TMath::Max(0., fXSecPipA_Elas    -> Evaluate (ke));
   else if (fate == kIHAFtInelas ) return TMath::Max(0., fXSecPipA_Inel    -> Evaluate (ke));
   else if (fate == kIHAFtAbsNP  ) return TMath::Max(0., fXSecPipA_NP      -> Evaluate (ke));
   else if (fate == kIHAFtAbsPP  ) return TMath::Max(0., fXSecPipA_PP      -> Evaluate (ke));
   else if (fate == kIHAFtAbsNPP ) return TMath::Max(0., fXSecPipA_NPP     -> Evaluate (ke));
   else if (fate == kIHAFtAbsNNP ) return TMath::Max(0., fXSecPipA_NNP     -> Evaluate (ke));
   else if (fate == kIHAFtAbs2N2P) return TMath::Max(0., fXSecPipA_NNPP    -> Evaluate (ke));
   else if (fate == kIHAFtNPipPi0) return TMath::Max(0., fXSecPipA_NPipPi0 -> Evaluate (ke));
   else {
     LOG("INukeData", pWARN) 
         << "Pi+'s don't have this fate: " << INukeHadroFates::AsString(fate);
     return 0;
   }

  } else if (hpdgc == kPdgPiM) {
   /* handle pi- */
   if      (fate == kIHAFtCEx    ) return TMath::Max(0., fXSecPimA_CEx     -> Evaluate (ke));
   else if (fate == kIHAFtElas   ) return TMath::Max(0., fXSecPimA_Elas    -> Evaluate (ke));
   else if (fate == kIHAFtInelas ) return TMath::Max(0., fXSecPimA_Inel    -> Evaluate (ke));
   else if (fate == kIHAFtAbsNP  ) return TMath::Max(0., fXSecPimA_NP      -> Evaluate (ke));
   else if (fate == kIHAFtAbsPP  ) return TMath::Max(0., fXSecPimA_PP      -> Evaluate (ke));
   else if (fate == kIHAFtAbsNPP ) return TMath::Max(0., fXSecPimA_NPP     -> Evaluate (ke));
   else if (fate == kIHAFtAbsNNP ) return TMath::Max(0., fXSecPimA_NNP     -> Evaluate (ke));
   else if (fate == kIHAFtAbs2N2P) return TMath::Max(0., fXSecPimA_NNPP    -> Evaluate (ke));
   else if (fate == kIHAFtNPipPi0) return TMath::Max(0., fXSecPimA_NPipPi0 -> Evaluate (ke));
   else {
     LOG("INukeData", pWARN) 
        << "Pi-'s don't have this fate: " << INukeHadroFates::AsString(fate);
     return 0;
   }

  } else if (hpdgc == kPdgPi0) {
   /* handle pi0 */
        if (fate == kIHAFtCEx    ) return TMath::Max(0., fXSecPi0A_CEx     -> Evaluate (ke));
   else if (fate == kIHAFtElas   ) return TMath::Max(0., fXSecPi0A_Elas    -> Evaluate (ke));
   else if (fate == kIHAFtInelas ) return TMath::Max(0., fXSecPi0A_Inel    -> Evaluate (ke));
   else if (fate == kIHAFtAbsNP  ) return TMath::Max(0., fXSecPi0A_NP      -> Evaluate (ke));
   else if (fate == kIHAFtAbsPP  ) return TMath::Max(0., fXSecPi0A_PP      -> Evaluate (ke));
   else if (fate == kIHAFtAbsNPP ) return TMath::Max(0., fXSecPi0A_NPP     -> Evaluate (ke));
   else if (fate == kIHAFtAbsNNP ) return TMath::Max(0., fXSecPi0A_NNP     -> Evaluate (ke));
   else if (fate == kIHAFtAbs2N2P) return TMath::Max(0., fXSecPi0A_NNPP    -> Evaluate (ke));
   else if (fate == kIHAFtNPipPi0) return TMath::Max(0., fXSecPi0A_NPipPi0 -> Evaluate (ke));
   else {
     LOG("INukeData", pWARN) 
        << "Pi0's don't have this fate: " << INukeHadroFates::AsString(fate);
       return 0;
   }
  }
  LOG("INukeData", pWARN) 
      << "Can't handle particles with pdg code = " << hpdgc;

  return 0;
}
//____________________________________________________________________________
double INukeHadroData::XSec(int hpdgc, INukeFateHN_t fate, double ke) const
{
// return the x-section for the input fate for the particle with the input pdg 
// code at the input kinetic energy
//
  ke = TMath::Max(fMinKinEnergy,   ke);
  ke = TMath::Min(fMaxKinEnergyHN, ke);

  if (hpdgc == kPdgProton) {  
    /* handle protons */
         if (fate == kIHNFtElas  ) return TMath::Max(0., fXSecPN_Elas -> Evaluate(ke));
    else if (fate == kIHNFtInelas) return TMath::Max(0., fXSecPN_Reac -> Evaluate(ke));
    else {
     LOG("INukeData", pWARN) 
        << "Protons don't have this fate: " << INukeHadroFates::AsString(fate);
     return 0;
    }

  } else if (hpdgc == kPdgNeutron) {
    /* handle neutrons */
         if (fate == kIHNFtElas  ) return TMath::Max(0., fXSecNN_Elas -> Evaluate(ke));
    else if (fate == kIHNFtInelas) return TMath::Max(0., fXSecNN_Reac -> Evaluate(ke));
    else {
     LOG("INukeData", pWARN) 
        << "Neutrons don't have this fate: " << INukeHadroFates::AsString(fate);
     return 0;
    }

  } else if (hpdgc == kPdgPiP) {
    /* handle pi+ */
         if (fate == kIHNFtCEx   ) return TMath::Max(0., fXSecPipN_CEx  -> Evaluate(ke));
    else if (fate == kIHNFtElas  ) return TMath::Max(0., fXSecPipN_Elas -> Evaluate(ke));
    else if (fate == kIHNFtInelas) return TMath::Max(0., fXSecPipN_Reac -> Evaluate(ke));
    else if (fate == kIHNFtAbs   ) return TMath::Max(0., fXSecPipN_Abs  -> Evaluate(ke));
    else {
     LOG("INukeData", pWARN) 
        << "Pi+'s don't have this fate: " << INukeHadroFates::AsString(fate);
     return 0;
    }

  } else if (hpdgc == kPdgPiM) {
    /* handle pi- */
         if (fate == kIHNFtCEx   ) return TMath::Max(0., fXSecPipN_CEx  -> Evaluate(ke));
    else if (fate == kIHNFtElas  ) return TMath::Max(0., fXSecPipN_Elas -> Evaluate(ke));
    else if (fate == kIHNFtInelas) return TMath::Max(0., fXSecPipN_Reac -> Evaluate(ke));
    else if (fate == kIHNFtAbs   ) return TMath::Max(0., fXSecPipN_Abs  -> Evaluate(ke));
    else {
     LOG("INukeData", pWARN) 
        << "Pi-'s don't have this fate: " << INukeHadroFates::AsString(fate);
     return 0;
    }

  } else if (hpdgc == kPdgPi0) {
    /* handle pi0 */
         if (fate == kIHNFtCEx   ) return TMath::Max(0., fXSecPipN_CEx  -> Evaluate(ke));
    else if (fate == kIHNFtElas  ) return TMath::Max(0., fXSecPipN_Elas -> Evaluate(ke));
    else if (fate == kIHNFtInelas) return TMath::Max(0., fXSecPipN_Reac -> Evaluate(ke));
    else if (fate == kIHNFtAbs   ) return TMath::Max(0., fXSecPipN_Abs  -> Evaluate(ke));
    else {
     LOG("INukeData", pWARN) 
        << "Pi0's don't have this fate: " << INukeHadroFates::AsString(fate);
     return 0;
    }
  }
  LOG("INukeData", pWARN) 
      << "Can't handle particles with pdg code = " << hpdgc;

  return 0;
}
//____________________________________________________________________________
double INukeHadroData::Frac(int hpdgc, INukeFateHA_t fate, double ke) const
{
// return the x-section fraction for the input fate for the particle with the 
// input pdg code at the input kinetic energy

  ke = TMath::Max(fMinKinEnergy,   ke);
  ke = TMath::Min(fMaxKinEnergyHA, ke);

  // get x-section
  double xsec = this->XSec(hpdgc,fate,ke);

  // get max x-section
  double xsec_tot = 0;
       if (hpdgc == kPdgProton ) xsec_tot = TMath::Max(0., fXSecPA_Tot   -> Evaluate (ke));
  else if (hpdgc == kPdgNeutron) xsec_tot = TMath::Max(0., fXSecNA_Tot   -> Evaluate (ke));
  else if (hpdgc == kPdgPiP    ) xsec_tot = TMath::Max(0., fXSecPipA_Tot -> Evaluate (ke));
  else if (hpdgc == kPdgPiM    ) xsec_tot = TMath::Max(0., fXSecPimA_Tot -> Evaluate (ke));
  else if (hpdgc == kPdgPi0    ) xsec_tot = TMath::Max(0., fXSecPi0A_Tot -> Evaluate (ke));

  // compute fraction
  double frac = (xsec_tot>0) ? xsec/xsec_tot : 0.;
  return frac;
}
//____________________________________________________________________________
double INukeHadroData::Frac(int hpdgc, INukeFateHN_t fate, double ke) const
{
// return the x-section fraction for the input fate for the particle with the 
// input pdg code at the input kinetic energy

  ke = TMath::Max(fMinKinEnergy,   ke);
  ke = TMath::Min(fMaxKinEnergyHN, ke);

  // get x-section
  double xsec = this->XSec(hpdgc,fate,ke);

  // get max x-section
  double xsec_tot = 0;
       if (hpdgc == kPdgProton ) xsec_tot = TMath::Max(0., fXSecPN_Tot   -> Evaluate (ke));
  else if (hpdgc == kPdgNeutron) xsec_tot = TMath::Max(0., fXSecNN_Tot   -> Evaluate (ke));
  else if (hpdgc == kPdgPiP    ) xsec_tot = TMath::Max(0., fXSecPipN_Tot -> Evaluate (ke));
  else if (hpdgc == kPdgPiM    ) xsec_tot = TMath::Max(0., fXSecPimN_Tot -> Evaluate (ke));
  else if (hpdgc == kPdgPi0    ) xsec_tot = TMath::Max(0., fXSecPi0N_Tot -> Evaluate (ke));

  // compute fraction
  double frac = (xsec_tot>0) ? xsec/xsec_tot : 0.;
  return frac;
}
//____________________________________________________________________________


