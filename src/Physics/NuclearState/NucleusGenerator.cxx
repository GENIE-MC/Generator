//____________________________________________________________________________
/*!

\class    genie::NucleusGenerator

\brief    It visits the event record & computes a Fermi motion momentum for
          initial state nucleons bound in nuclei.
          Is a concrete implementation of the EventRecordVisitorI interface.

\author   Liang Liu <liangliu \at fnal.gov>
          Fermi National Accelerator Laboratory

\created  October 17, 2024

\cpright  Copyright (c) 2003-2024, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          
*/
//____________________________________________________________________________

#include <cstdlib>

#include <TLorentzVector.h>
#include <TVector3.h>
#include <TParticlePDG.h>
#include <TMath.h>

#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Units.h"
#include "Physics/NuclearState/NucleusGenerator.h"

using namespace genie;
using namespace genie::constants;

//___________________________________________________________________________
NucleusGenerator::NucleusGenerator() :
EventRecordVisitorI("genie::NucleusGenerator")
{

}
//___________________________________________________________________________
NucleusGenerator::NucleusGenerator(string config) :
EventRecordVisitorI("genie::NucleusGenerator", config)
{

}
//___________________________________________________________________________
NucleusGenerator::~NucleusGenerator()
{

}

//___________________________________________________________________________
void NucleusGenerator::ProcessEventRecord(GHepRecord * evrec) const
{
  // skip if not a nuclear target
  //if(! evrec->Summary()->InitState().Tgt().IsNucleus()) return;
  LOG("NucleusGenerator", pINFO) << "Adding final state nucleus";
  fNucleusGen->ProcessEventRecord(evrec);
}

//___________________________________________________________________________
void NucleusGenerator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void NucleusGenerator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void NucleusGenerator::LoadConfig(void)
{
  RgKey nuclkey = "NucleusGenerator";
  fNucleusGen = nullptr;
  //fNuclModel = dynamic_cast<const NucleusGenerator *> (this->SubAlg(nuclkey));
  //
  LOG("NucleusGenerator", pINFO) << "Hello world!";
  fNucleusGen = dynamic_cast<const EventRecordVisitorI *> (this->SubAlg("Nuclear-Model"));
}
//____________________________________________________________________________

