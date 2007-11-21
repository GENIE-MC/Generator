//____________________________________________________________________________
/*
 Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         STFC, Rutherford Appleton Laboratory - October 03, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include "Algorithm/AlgConfigPool.h"
#include "EVGModules/CorrelatedNucleonGenerator.h"
#include "GHEP/GHepStatus.h"
#include "GHEP/GHepParticle.h"
#include "GHEP/GHepRecord.h"
#include "Messenger/Messenger.h"
#include "Nuclear/FermiMomentumTablePool.h"
#include "Nuclear/FermiMomentumTable.h"

using namespace genie;

//___________________________________________________________________________
CorrelatedNucleonGenerator::CorrelatedNucleonGenerator() :
EventRecordVisitorI("genie::CorrelatedNucleonGenerator")
{

}
//___________________________________________________________________________
CorrelatedNucleonGenerator::CorrelatedNucleonGenerator(string config) :
EventRecordVisitorI("genie::CorrelatedNucleonGenerator", config)
{

}
//___________________________________________________________________________
CorrelatedNucleonGenerator::~CorrelatedNucleonGenerator()
{

}
//___________________________________________________________________________
void CorrelatedNucleonGenerator::ProcessEventRecord(GHepRecord * evrec) const
{
  if(!fSimulateCorrelN) return;

//  LOG("NNCorrel", pNOTICE) << "Simulating NN correlation";
//  LOG("NNCorrel", pNOTICE) << *evrec;

  GHepParticle * nucltgt = evrec->TargetNucleus();
  if (!nucltgt) {
    LOG("NNCorrel", pINFO) << "No nuclear target found - Doing nothing";
    return;
  }

  GHepParticle * hitnucl = evrec->HitNucleon();
  if(!hitnucl) {
    LOG("NNCorrel", pINFO) << "No hit nucleon found - Doing nothing";
    return;
  }

  // hit nuclear target & nucleon pdg codes
  int target_pdgc  = nucltgt->Pdg(); 
  int nucleon_pdgc = hitnucl->Pdg();

  // check the actual hit nucleon momentum momentum 
  double Pn = hitnucl->P4()->Mag();

  // get kF
  string fKFTable = "Default";
  FermiMomentumTablePool * kftp = FermiMomentumTablePool::Instance();
  const FermiMomentumTable * kft  = kftp->GetTable(fKFTable);
  double kF = kft->FindClosestKF(target_pdgc, nucleon_pdgc);

  // decide whether to eject extra init state nucleon
  bool eject = (Pn > kF);

  // add an extra nucleon
  if(eject) {

    LOG("NNCorrel", pNOTICE) << "Ejecting extra nucleon";

    double px = -1 * hitnucl->Px();
    double py = -1 * hitnucl->Py();
    double pz = -1 * hitnucl->Pz();
    double E  =      hitnucl->E();
//  double w  =      hitnucl->RemovalEnergy();

    TLorentzVector p4(px,py,pz,E);
    TLorentzVector v4(0,0,0,0);

    evrec->AddParticle(
       nucleon_pdgc, kIStCorrelatedNucleon, 1, 1, -1, -1, p4, v4);

//    LOG("NNCorrel", pNOTICE) << *evrec;
  }
}
//___________________________________________________________________________
void CorrelatedNucleonGenerator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void CorrelatedNucleonGenerator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void CorrelatedNucleonGenerator::LoadConfig(void)
{
  AlgConfigPool * confp = AlgConfigPool::Instance();
  const Registry * gc = confp->GlobalParameterList();

  fSimulateCorrelN = fConfig->GetBoolDef(
              "Enable", gc->GetBool("CorrelNN-Enable"));
/*
  fMomentumThr  = fConfig->GetDoubleDef(
           "MomentumThreshold", gc->GetDouble("CorrelNN-MomentumThreshold"));
*/
}
//____________________________________________________________________________

