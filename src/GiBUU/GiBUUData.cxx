//____________________________________________________________________________
/*
 Copyright (C) 2003-2015, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         May 30, 2009

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ May 30, 2009 - CA
   Was first added in v2.5.1
 @ Sep 29, 2009 - CA
   Add FormFactors nested class & re-organize code
 @ May 01, 2012 - CA
   Pick data from new location ($GENIE/data/evgen/gibuu/)
*/
//____________________________________________________________________________

#include <cassert>
#include <sstream>
#include <string>

#include <TSystem.h>
#include <TTree.h>
#include <TMath.h>

#include "BaryonResonance/BaryonResUtils.h"
#include "GiBUU/GiBUUData.h"
#include "Messenger/Messenger.h"
#include "Numerical/Spline.h"
#include "PDG/PDGCodes.h"

using std::ostringstream;
using std::istream;
using std::ios;
using std::cout;
using std::endl;

using namespace genie;
using namespace genie::utils;

//____________________________________________________________________________
GiBUUData * GiBUUData::fInstance = 0;
//____________________________________________________________________________
GiBUUData::GiBUUData()
{
  this->LoadTables();
  fInstance = 0;
}
//____________________________________________________________________________
GiBUUData::~GiBUUData()
{
  delete fFormFactors;
  fFormFactors = 0;
}
//____________________________________________________________________________
GiBUUData * GiBUUData::Instance()
{
  if(fInstance == 0) {
    LOG("GiBUUData", pINFO) << "GiBUUData late initialization";
    static GiBUUData::Cleaner cleaner;
    cleaner.DummyMethodAndSilentCompiler();
    fInstance = new GiBUUData;
  }
  return fInstance;
}
//____________________________________________________________________________
void GiBUUData::LoadTables(void)
{
  fFormFactors = new FormFactors;
  fFormFactors->LoadTables();
}
//____________________________________________________________________________
const GiBUUData::FormFactors & GiBUUData::FF(void) const
{
  return *fFormFactors;
}
//____________________________________________________________________________
//
// FormFactors nested class
//
//____________________________________________________________________________
double GiBUUData::FormFactors::fMinQ2 = 0.0; // GeV^2
double GiBUUData::FormFactors::fMaxQ2 = 4.0; // GeV^2
//____________________________________________________________________________
GiBUUData::FormFactors::FormFactors(void)
{

}
//____________________________________________________________________________
GiBUUData::FormFactors::~FormFactors(void)
{
  if(!gAbortingInErr) {
    cout << "GiBUUData singleton dtor: Deleting all f/f splines" << endl;
  }

  // resonance form factor splines
  for(int r=0; r<kNRes; r++) {
    for(int i=0; i<kNCurr; i++) {
      for(int j=0; j<kNHitNuc; j++) {
        for(int k=0; k<kNFFRes; k++) {
           if (fFFRes[r][i][j][k]) {
             delete fFFRes[r][i][j][k];
             fFFRes[r][i][j][k] = 0;
           }
        }//k
      }//j
    }//i   
  }//r 
}
//____________________________________________________________________________
void GiBUUData::FormFactors::LoadTables(void)
{
// Loads hadronic x-section data

  for(int r=0; r<kNRes; r++) {
    for(int i=0; i<kNCurr; i++) {
      for(int j=0; j<kNHitNuc; j++) {
        for(int k=0; k<kNFFRes; k++) {
             fFFRes[r][i][j][k] = 0;
        }//k
      }//j
    }//i   
  }//r 

  string data_dir = string(gSystem->Getenv("GENIE")) + 
                    string("/data/evgen/gibuu");

  LOG("GiBUUData", pNOTICE) << "Loading GiBUU data from: " << data_dir;

  //
  // load resonance form factor data
  //
  for(int r=0; r<kNRes; r++) {
    for(int i=0; i<kNCurr; i++) {
      for(int j=0; j<kNHitNuc; j++) {

        bool skip = false;

        Resonance_t resonance = (Resonance_t)r;

        ostringstream datafile;
        datafile << data_dir << "/form_factors/";

        switch(r){
          case ( 0): datafile << "P33_1232"; break;
          case ( 1): datafile << "S11_1535"; break;
          case ( 2): datafile << "D13_1520"; break;
          case ( 3): datafile << "S11_1650"; break;
          case ( 5): datafile << "D15_1675"; break;
          case ( 6): datafile << "S31_1620"; break;
          case ( 7): datafile << "D33_1700"; break;
          case ( 8): datafile << "P11_1440"; break;
          case (10): datafile << "P13_1720"; break;
          case (11): datafile << "F15_1680"; break;
          case (12): datafile << "P31_1910"; break;
          case (14): datafile << "F35_1905"; break;
          case (15): datafile << "F37_1950"; break;
          default  : skip = true;

        }
        switch(i){
          case (0): datafile << "_CC"; break;
          case (1): datafile << "_NC"; break;
          case (2): datafile << "_EM"; break;
          default : skip = true;
        }
        switch(j){
          case (0): datafile << "_neutron"; break;
          case (1): datafile << "_proton";  break;
          default : skip = true;
        }
        datafile << "_FFres.dat";

        if(skip) continue;

        //-- Make sure that all data file exists
        assert( ! gSystem->AccessPathName(datafile.str().c_str()) );

        //-- Read the data and fill a tree

        //   The GiBUU files have the data organized in 9 columns:
        //              0    1      2      3      4      5      6      7      8
        //   I=3/2 res: Qs   C_3^V  C_4^V  C_5^V  C_6^V  C_3^A  C_4^A  C_5^A  C_6^A
        //   I=1/2 res: Qs   F_1^V  F_2^V  -----  -----  F_A    F_P    -----  ----

        TTree data_ffres;
        data_ffres.ReadFile(datafile.str().c_str(), 
                            "Q2/D:f1/D:f2/D:f3/D:f4/D:f5/D:f6/D:f7/D:f8/D");

        LOG("GiBUUData", pDEBUG)  
           << "Number of data rows: " << data_ffres.GetEntries();

        //
        // I=3/2 resonances
        //
        if(res::IsDelta(resonance)) {
           fFFRes[r][i][j][0]  = new Spline(&data_ffres, "Q2:f1"); // F1V = f(Q2)
           fFFRes[r][i][j][1]  = new Spline(&data_ffres, "Q2:f2"); // F2V = f(Q2)
           fFFRes[r][i][j][2]  = new Spline(&data_ffres, "Q2:f5"); // FA  = f(Q2)
           fFFRes[r][i][j][3]  = new Spline(&data_ffres, "Q2:f6"); // FP  = f(Q2)
        } // Delta res
        else 
        //
        // I=1/2 resonances
        //
        if(res::IsN(resonance)) {
           fFFRes[r][i][j][4]  = new Spline(&data_ffres, "Q2:f1"); // C3V = f(Q2)
           fFFRes[r][i][j][5]  = new Spline(&data_ffres, "Q2:f2"); // C4V = f(Q2)
           fFFRes[r][i][j][6]  = new Spline(&data_ffres, "Q2:f3"); // C5V = f(Q2)
           fFFRes[r][i][j][7]  = new Spline(&data_ffres, "Q2:f4"); // C6V = f(Q2)
           fFFRes[r][i][j][8]  = new Spline(&data_ffres, "Q2:f5"); // C3A = f(Q2)
           fFFRes[r][i][j][9]  = new Spline(&data_ffres, "Q2:f6"); // C4A = f(Q2)
           fFFRes[r][i][j][10] = new Spline(&data_ffres, "Q2:f7"); // C5A = f(Q2)
           fFFRes[r][i][j][11] = new Spline(&data_ffres, "Q2:f8"); // C6A = f(Q2)
        } //N res

      }//j
    }//i   
  }//r 


  for(int r=0; r<kNRes; r++) {
    for(int i=0; i<kNCurr; i++) {
      for(int j=0; j<kNHitNuc; j++) {
        for(int k=0; k<kNFFRes; k++) {
             if(fFFRes[r][i][j][k]) fFFRes[r][i][j][k]->YCanBeNegative(true);
        }//k
      }//j
    }//i   
  }//r 

  LOG("GiBUUData", pINFO)  
     << "Done loading all resonance form factor files...";
}
//____________________________________________________________________________
double GiBUUData::FormFactors::C3V(
  double Q2, Resonance_t res, int hit_nucleon_pdg, InteractionType_t it) const
{
  if(!res::IsN(res)) return 0.;
  return this->FFRes(Q2,res,hit_nucleon_pdg,it,4);
}
//____________________________________________________________________________
double GiBUUData::FormFactors::C4V(
  double Q2, Resonance_t res, int hit_nucleon_pdg, InteractionType_t it) const
{
  if(!res::IsN(res)) return 0.;
  return this->FFRes(Q2,res,hit_nucleon_pdg,it,5);
}
//____________________________________________________________________________
double GiBUUData::FormFactors::C5V(
  double Q2, Resonance_t res, int hit_nucleon_pdg, InteractionType_t it) const
{
  if(!res::IsN(res)) return 0.;
  return this->FFRes(Q2,res,hit_nucleon_pdg,it,6);
}
//____________________________________________________________________________
double GiBUUData::FormFactors::C6V(
  double Q2, Resonance_t res, int hit_nucleon_pdg, InteractionType_t it) const
{
  if(!res::IsN(res)) return 0.;
  return this->FFRes(Q2,res,hit_nucleon_pdg,it,7);
}
//____________________________________________________________________________
double GiBUUData::FormFactors::C3A(
  double Q2, Resonance_t res, int hit_nucleon_pdg, InteractionType_t it) const
{
  if(!res::IsN(res)) return 0.;
  return this->FFRes(Q2,res,hit_nucleon_pdg,it,8);
}
//____________________________________________________________________________
double GiBUUData::FormFactors::C4A(
  double Q2, Resonance_t res, int hit_nucleon_pdg, InteractionType_t it) const
{
  if(!res::IsN(res)) return 0.;
  return this->FFRes(Q2,res,hit_nucleon_pdg,it,9);
}
//____________________________________________________________________________
double GiBUUData::FormFactors::C5A(
  double Q2, Resonance_t res, int hit_nucleon_pdg, InteractionType_t it) const
{
  if(!res::IsN(res)) return 0.;
  return this->FFRes(Q2,res,hit_nucleon_pdg,it,10);
}
//____________________________________________________________________________
double GiBUUData::FormFactors::C6A(
  double Q2, Resonance_t res, int hit_nucleon_pdg, InteractionType_t it) const
{
  if(!res::IsN(res)) return 0.;
  return this->FFRes(Q2,res,hit_nucleon_pdg,it,11);
}
//____________________________________________________________________________
double GiBUUData::FormFactors::F1V(
  double Q2, Resonance_t res, int hit_nucleon_pdg, InteractionType_t it) const
{
  if(!res::IsDelta(res)) return 0.;
  return this->FFRes(Q2,res,hit_nucleon_pdg,it,0);
}
//____________________________________________________________________________
double GiBUUData::FormFactors::F2V(
  double Q2, Resonance_t res, int hit_nucleon_pdg, InteractionType_t it) const
{
  if(!res::IsDelta(res)) return 0.;
  return this->FFRes(Q2,res,hit_nucleon_pdg,it,1);
}
//____________________________________________________________________________
double GiBUUData::FormFactors::FA(
  double Q2, Resonance_t res, int hit_nucleon_pdg, InteractionType_t it) const
{
  if(!res::IsDelta(res)) return 0.;
  return this->FFRes(Q2,res,hit_nucleon_pdg,it,2);
}
//____________________________________________________________________________
double GiBUUData::FormFactors::FP(
  double Q2, Resonance_t res, int hit_nucleon_pdg, InteractionType_t it) const
{
  if(!res::IsDelta(res)) return 0.;
  return this->FFRes(Q2,res,hit_nucleon_pdg,it,3);
}
//____________________________________________________________________________
double GiBUUData::FormFactors::FFRes ( 
    double Q2, Resonance_t res, int hit_nucleon_pdg, InteractionType_t it, 
    int ffresid) const
{
  if(Q2 < fMinQ2 || Q2 > fMaxQ2) return 0.;

  int r = -1, i = -1, j = -1;

  if(ffresid<0 || ffresid >= kNFFRes) return 0.;

  r = (int)res;
  if(r<0 || r >= kNRes) return 0.;

  if      (it == kIntWeakCC) { i = 0; }
  else if (it == kIntWeakNC) { i = 1; }
  else if (it == kIntEM)     { i = 2; }

  if      (hit_nucleon_pdg == kPdgNeutron) { j = 0; }
  else if (hit_nucleon_pdg == kPdgProton ) { j = 1; }

  const Spline * spl = fFFRes[r][i][j][ffresid];

  if(!spl) return 0;
  else {
   return spl->Evaluate(Q2);
  }
}
//____________________________________________________________________________
