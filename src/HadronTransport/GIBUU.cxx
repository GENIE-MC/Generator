//____________________________________________________________________________
/*
 Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author:   
    Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk> STFC, Rutherford Lab
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

#include "Algorithm/AlgConfigPool.h"
#include "Conventions/Constants.h"
#include "EVGCore/EVGThreadException.h"
#include "GHEP/GHepRecord.h"
#include "GHEP/GHepParticle.h"
#include "GHEP/GHepStatus.h"
#include "HadronTransport/GIBUU.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"
#include "Utils/PrintUtils.h"

using namespace genie;
using namespace genie::utils;
using namespace genie::constants;

//___________________________________________________________________________
GIBUU::GIBUU() :
EventRecordVisitorI("genie::GIBUU")
{

}
//___________________________________________________________________________
GIBUU::GIBUU(string config) :
EventRecordVisitorI("genie::GIBUU", config)
{

}
//___________________________________________________________________________
GIBUU::~GIBUU()
{

}
//___________________________________________________________________________
void GIBUU::ProcessEventRecord(GHepRecord * event) const
{
  //-- Check that we have an interaction with a nuclear target. If not skip...
  GHepParticle * nucltgt = event->TargetNucleus();
  if (!nucltgt) {
    LOG("GIBUU", pINFO)
       << "No nuclear target found - GIBUU will not be called....";
    return;
  }

#ifdef __GENIE_GIBUU_ENABLED__

  //-- Preparing event before  passing it to GiBUU 

  // Set intranuke-style formation zones. Within GiBUU formation zones are 
  // taken from the JETSET model but, when GiBUU is called from GENIE, that
  // step is bypassed
  this->SetFormationZones(event);

  //-- Translate GENIE GHepRecord to whatever GIBUU needs

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
    int   in_pdg = p->Pdg();
    float in_px  = (float) p->Px();
    float in_py  = (float) p->Py();
    float in_pz  = (float) p->Pz();
    float in_E   = (float) p->E();
    float in_x   = (float) p->Vx();
    float in_y   = (float) p->Vy();
    float in_z   = (float) p->Vz();
    float in_t   = (float) p->Vt();

    AddHadron(&ihad, &in_pdg, &in_px, &in_py, &in_pz, &in_E, &in_x, &in_y, &in_z, &in_t);

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
  
        // inhibit bindong energy removal from latter steps - that has been
        // taken into account internally at GiBUU
        new_particle->SetRemovalEnergy(0.);

        event->AddParticle(new_particle);

      } // input hadron daughters (output hadrons)
  }//input hadrons

#else
  LOG("GIBUU", pFATAL) 
       << "\n"
       << "\n****************************************************"
       << "\n*** YOU HAVE NOT INSTALLED GIBUU!                ***"
       << "\n*** Please obtain the actual GiBUU code from:    ***"
       << "\n*** http://tp8.physik.uni-giessen.de:8080/GiBUU/ ***"
       << "\n****************************************************"
       << "\n";
  exit(1);
#endif
}
//___________________________________________________________________________
void GIBUU::SetFormationZones(GHepRecord * event) const
{
// Set intranuke-style formation zones to 'initial' hadrons

  // Get hadronic system's 3-momentum
  GHepParticle * hadronic_system = event->FinalStateHadronicSystem();
  TVector3 p3hadr = hadronic_system->P4()->Vect(); // (px,py,pz)

  TObjArrayIter piter(event);
  GHepParticle * p = 0;
  while( (p = (GHepParticle *) piter.Next()) )
  {
     // Use only particles marked as hadrons in the nucleus
     GHepStatus_t ist = p->Status();
     if(ist != kIStHadronInTheNucleus) continue;

     // Compute formation zone
     TVector3 p3  = p->P4()->Vect();      // hadron's: p (px,py,pz)
     double   m   = p->Mass();            //           m
     double   m2  = m*m;                  //           m^2
     double   P   = p->P4()->P();         //           |p|
     double   Pt  = p3.Pt(p3hadr);        //           pT
     double   Pt2 = Pt*Pt;                //           pT^2
     double   fz  = P*fct0*m/(m2+fK*Pt2); //           formation zone, in m

     LOG("GIBUU", pINFO)
          << "|P| = " << P << " GeV, Pt = " << Pt
                              << " GeV, Formation Zone = " << fz << " m";

     // Step particle 
     TVector3 dr = p->P4()->Vect().Unit();          // unit vector along its direction
     double c  = kLightSpeed / (units::m/units::s); // c in m/sec
     dr.SetMag(fz);                                 // spatial step size
     double dt = fz/c;                              // temporal step:
     TLorentzVector dx4(dr,dt);                     // 4-vector step
     TLorentzVector x4new = *(p->X4()) + dx4;       // new position

     LOG("GIBUU", pDEBUG)
        << "\n Init direction = " << print::Vec3AsString(&dr) 
        << "\n Init position (in m,sec) = " << print::X4AsString(p->X4())
        << "\n Fin  position (in m,sec) = " << print::X4AsString(&x4new);

     p->SetPosition(x4new);
  }
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

  AlgConfigPool * confp = AlgConfigPool::Instance();
  const Registry * gc = confp->GlobalParameterList();
 
  fct0     = fConfig->GetDoubleDef ("ct0",  gc->GetDouble("INUKE-FormationZone")); // fm
  fK       = fConfig->GetDoubleDef ("Kpt2", gc->GetDouble("INUKE-KPt2"));
}
//___________________________________________________________________________

