//____________________________________________________________________________
/*
 Copyright (c) 2003-2016, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab 

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Feb 09, 2009 - CA
   Moved into the new Coherent package from its previous location  (EVGModules 
   package)

*/
//____________________________________________________________________________

#include <TMath.h>
#include <TVector3.h>

#include "Conventions/Constants.h"
#include "Coherent/COHPrimaryLeptonGenerator.h"
#include "GHEP/GHepParticle.h"
#include "GHEP/GHepRecord.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"
#include "Numerical/RandomGen.h"
#include "Base/XSecAlgorithmI.h"
#include "EVGCore/EVGThreadException.h"
#include "EVGCore/EventGeneratorI.h"
#include "EVGCore/RunningThreadInfo.h"

using namespace genie;
using namespace genie::constants;

//___________________________________________________________________________
COHPrimaryLeptonGenerator::COHPrimaryLeptonGenerator() :
  PrimaryLeptonGenerator("genie::COHPrimaryLeptonGenerator")
{

}
//___________________________________________________________________________
COHPrimaryLeptonGenerator::COHPrimaryLeptonGenerator(string config) :
  PrimaryLeptonGenerator("genie::COHPrimaryLeptonGenerator", config)
{

}
//___________________________________________________________________________
COHPrimaryLeptonGenerator::~COHPrimaryLeptonGenerator()
{

}
//___________________________________________________________________________
void COHPrimaryLeptonGenerator::ProcessEventRecord(GHepRecord * evrec) const
{
  //-- Access cross section algorithm for running thread
  RunningThreadInfo * rtinfo = RunningThreadInfo::Instance();
  const EventGeneratorI * evg = rtinfo->RunningThread();
  const XSecAlgorithmI *fXSecModel = evg->CrossSectionAlg();

  // In Rein and Berger-Sehgal, no modification is required to the standard impl.
  if (fXSecModel->Id().Name() == "genie::ReinSehgalCOHPiPXSec") {
    PrimaryLeptonGenerator::ProcessEventRecord(evrec);
  }
  else if ((fXSecModel->Id().Name() == "genie::BergerSehgalCOHPiPXSec2015")) {
    PrimaryLeptonGenerator::ProcessEventRecord(evrec);
  }
  else if ((fXSecModel->Id().Name() == "genie::BergerSehgalFMCOHPiPXSec2015")) {
    PrimaryLeptonGenerator::ProcessEventRecord(evrec);
  }
  else if ((fXSecModel->Id().Name() == "genie::AlvarezRusoCOHPiPXSec")) {
    CalculatePrimaryLepton_AlvarezRuso(evrec);
  }
  else {
    LOG("COHPrimaryLeptonGenerator",pFATAL) <<
      "ProcessEventRecord >> Cannot calculate primary lepton for " <<
      fXSecModel->Id().Name();
  }
}
//___________________________________________________________________________
void COHPrimaryLeptonGenerator::CalculatePrimaryLepton_AlvarezRuso(GHepRecord * evrec) const
{
  Interaction * interaction = evrec->Summary();
  const Kinematics & kinematics = interaction->Kine();
  TLorentzVector p4l = kinematics.FSLeptonP4();
  int pdgc = interaction->FSPrimLepton()->PdgCode();
  this->AddToEventRecord(evrec, pdgc, p4l);
}
//___________________________________________________________________________
