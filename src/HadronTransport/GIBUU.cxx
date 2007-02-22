//____________________________________________________________________________
/*
 Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author:   
    Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk> CCLRC, Rutherford Lab
    Tina Leitner <Tina.J.Leitner@theo.physik.uni-giessen.de> Giessen Univ.
    February 22, 2007

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <cstdlib>

#include <TLorentzVector.h>
#include <TVector3.h>
#include <TSystem.h>

#include "EVGCore/EVGThreadException.h"
#include "GHEP/GHepRecord.h"
#include "GHEP/GHepParticle.h"
#include "GHEP/GHepStatus.h"
#include "HadronTransport/GIBUU.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"

using namespace genie;

//___________________________________________________________________________
GIBUU::GIBUU() :
EventRecordVisitorI("genie::GIBUU")
{
  this->CheckInstallation();
}
//___________________________________________________________________________
GIBUU::GIBUU(string config) :
EventRecordVisitorI("genie::GIBUU", config)
{
  this->CheckInstallation();
}
//___________________________________________________________________________
GIBUU::~GIBUU()
{

}
//___________________________________________________________________________
void GIBUU::ProcessEventRecord(GHepRecord * event) const
{
  if(!fIsEnabled) {
    LOG("GIBUU", pFATAL) 
       << "\n"
       << "\n****************************************************"
       << "\n*** YOU HAVE NOT INSTALLED GIBUU!                ***"
       << "\n*** Please obtain the actual GiBUU code from:    ***"
       << "\n*** http://tp8.physik.uni-giessen.de:8080/GiBUU/ ***"
       << "\n****************************************************"
       << "\n";
    exit(1);
  }

  //-- Check that we have an interaction with a nuclear target. If not skip...
  GHepParticle * nucltgt = event->TargetNucleus();
  if (!nucltgt) {
    LOG("GIBUU", pINFO)
       << "No nuclear target found - GIBUU will not be called....";
    return;
  }

  //-- Translate GENIE GHepRecord to whatever FLUKA needs

  LOG("GIBUU", pDEBUG) << "Translating: GENIE GHepRecord ---> GIBUU input";

  // set nuclear target
  int target_pdg = nucltgt->Pdg();
  SetNucleus(&target_pdg);

  // set hadrons to be transported out
  int ihad=-1;
  int ipos=-1;
  map<int,int> mother_pos; // keep track of ith input hadron's position

  TObjArrayIter piter(event);
  GHepParticle * p = 0;
  while( (p = (GHepParticle *) piter.Next()) )
  {
    ipos++;
    GHepStatus_t ist = p->Status();
    if(ist != kIStHadronInTheNucleus) continue;

    ihad++;
    int   inp_pdg = p->Pdg();
    float inp_px  = (float) p->Px();
    float inp_py  = (float) p->Py();
    float inp_pz  = (float) p->Pz();
    float inp_E   = (float) p->E();

    AddHadron(&ihad, &inp_pdg, &inp_px, &inp_py, &inp_pz, &inp_E);

    mother_pos[ihad] = ipos;
  }  
  int n_inp_had = mother_pos.size();
    
  //-- Run GIBUU (use GiBUU's GENIE hook)
  LOG("GIBUU", pNOTICE) << "************ Running GIBUU ************";

  HadronTransportForGENIE();

  //-- Get FLUKA output & add it to GENIE GHepRecord
  LOG("GIBUU", pDEBUG) << "Copying: GIBUU output ---> GENIE GHepRecord";

  GHepStatus_t out_ist = kIStStableFinalState;
  for(int inp_had=0; inp_had < n_inp_had; inp_had++) {

      int n_daughters = NumOfFinalStateHadrons(&inp_had);
      int mom         = mother_pos[inp_had];

      for(int out_had=0; out_had < n_daughters; out_had++) {

        int   out_pdg = FinalStateHadronPDG (&inp_had, &out_had);
        float out_px  = FinalStateHadronPX  (&inp_had, &out_had);
        float out_py  = FinalStateHadronPY  (&inp_had, &out_had);
        float out_pz  = FinalStateHadronPZ  (&inp_had, &out_had);
        float out_e   = FinalStateHadronE   (&inp_had, &out_had);
        float out_x   = FinalStateHadronX   (&inp_had, &out_had);
        float out_y   = FinalStateHadronY   (&inp_had, &out_had);
        float out_z   = FinalStateHadronZ   (&inp_had, &out_had);
        float out_t   = FinalStateHadronT   (&inp_had, &out_had);

        GHepParticle new_particle(
             out_pdg, out_ist, mom,-1,-1,-1, 
                     out_px,out_py,out_pz,out_e, out_x,out_y,out_z,out_t);
  
        event->AddParticle(new_particle);

      } // input hadron daughters (output hadrons)
  }//input hadrons
}
//___________________________________________________________________________
void GIBUU::CheckInstallation(void)
{
  // check that GiBUU was enabled
  bool enabled = false;
  if(gSystem->Getenv("GOPT_GIBUU_ENABLE")) {
     string gibuu_enable = string(gSystem->Getenv("GOPT_GIBUU_ENABLE"));
     enabled = (gibuu_enable=="YES") || (gibuu_enable=="yes");
  }
  if(!enabled) {
    fIsEnabled = false;
    return;
  }

  // check that the GiBUU library exists
  bool lib_ok = false;
  if(gSystem->Getenv("GIBUU_LIB")) {
    string gibuu_lib = string(gSystem->Getenv("GIBUU_LIB")) + "/libGiBUU.a";
    lib_ok  = ! (gSystem->AccessPathName(gibuu_lib.c_str()));
  }
  if(!lib_ok) {
    fIsEnabled = false;
    return;
  }

  fIsEnabled = true;
}
//___________________________________________________________________________
void GIBUU::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//___________________________________________________________________________
void GIBUU::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//___________________________________________________________________________
void GIBUU::LoadConfig (void)
{
// Access this module's configuration options from its designated Registry
// and pass them to the actual GIBUU code

}
//___________________________________________________________________________
