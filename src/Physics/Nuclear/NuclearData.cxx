//____________________________________________________________________________
/*
 Copyright (C) 2003-2015, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Nov 30, 2009 - CA
   Was first added in v2.5.1
 @ May 01, 2012 - CA
   Pick-up data from new location ($GENIE/data/evgen/nucl/)
*/
//____________________________________________________________________________

#include <cassert>
#include <iostream>

#include <TSystem.h>
#include <TNtupleD.h>
#include <TTree.h>

#include "Messenger/Messenger.h"
#include "Nuclear/NuclearData.h"
#include "Numerical/Spline.h"
#include "PDG/PDGCodes.h"

using std::cout;
using std::endl;

using namespace genie;

//____________________________________________________________________________
NuclearData * NuclearData::fInstance = 0;
//____________________________________________________________________________
NuclearData::NuclearData()
{
  this->Load();
  fInstance = 0;
}
//____________________________________________________________________________
NuclearData::~NuclearData()
{
  if(!gAbortingInErr) {
    cout << "NuclearData singleton dtor: Deleting inputs... " << endl;
  }

  delete fNuclSupprD2;
}
//____________________________________________________________________________
NuclearData * NuclearData::Instance()
{
  if(fInstance == 0) {
    LOG("NuclData", pINFO) << "NuclearData late initialization";
    static NuclearData::Cleaner cleaner;
    cleaner.DummyMethodAndSilentCompiler();
    fInstance = new NuclearData;
  }
  return fInstance;
}
//____________________________________________________________________________
double NuclearData::DeuteriumSuppressionFactor(double Q2)
{
  if(Q2 > 0.20) return 1.; // no suppression

  if(Q2 < 3E-5) { Q2 = 3E-5; }

  double R = fNuclSupprD2->Evaluate(Q2);
  return R;
}
//____________________________________________________________________________
void NuclearData::Load(void)
{
  fNuclSupprD2 = 0;

  string data_dir = 
     string(gSystem->Getenv("GENIE")) + 
     string("/data/evgen/nucl");
  LOG("NuclData", pINFO)  
     << "Loading nuclear data from: " << data_dir;

  // Load D2 nuclear suppression factor

  string nuclsupprd2_file = data_dir + "/D2sup.data";
  assert( ! gSystem->AccessPathName(nuclsupprd2_file.c_str()) );

  TTree nuclsupprd2_tree;
  nuclsupprd2_tree.ReadFile(nuclsupprd2_file.c_str(), "Q2/D:R/D");

  fNuclSupprD2 = new Spline(&nuclsupprd2_tree, "Q2:R");

  LOG("NuclData", pINFO)  << "Done loading all data";
}
//____________________________________________________________________________
