//____________________________________________________________________________
/*!

\class    genie::NucleusGenTraditional

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
#include "Physics/NuclearState/NucleusGenTraditional.h"

#include "Physics/NuclearState/NuclearModel.h"
#include "Physics/NuclearState/NuclearModelI.h"
#include "Framework/EventGen/EVGThreadException.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepStatus.h"
#include "Framework/GHEP/GHepFlags.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Messenger/Messenger.h"
#include "Physics/NuclearState/FermiMomentumTablePool.h"
#include "Physics/NuclearState/FermiMomentumTable.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/Utils/KineUtils.h"
#include "Physics/NuclearState/NuclearUtils.h"
#include "Physics/NuclearState/INCLNucleus.h"

using namespace genie;
using namespace genie::constants;

//___________________________________________________________________________
NucleusGenTraditional::NucleusGenTraditional() :
EventRecordVisitorI("genie::NucleusGenTraditional")
{

}
//___________________________________________________________________________
NucleusGenTraditional::NucleusGenTraditional(string config) :
EventRecordVisitorI("genie::NucleusGenTraditional", config)
{

}
//___________________________________________________________________________
NucleusGenTraditional::~NucleusGenTraditional()
{

}

//___________________________________________________________________________
void NucleusGenTraditional::ProcessEventRecord(GHepRecord * evrec) const
{
  // skip if not a nuclear target
  if(! evrec->Summary()->InitState().Tgt().IsNucleus()) return;

  LOG("NucleusGenTraditional", pINFO) << "Adding final state nucleus";

  fVertexGenerator->ProcessEventRecord(evrec);
  fFermiMover->ProcessEventRecord(evrec);

}

//___________________________________________________________________________
void NucleusGenTraditional::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void NucleusGenTraditional::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void NucleusGenTraditional::LoadConfig(void)
{


  fFermiMover = nullptr;
  fVertexGenerator = nullptr;
  fFermiMover = dynamic_cast<const EventRecordVisitorI *> (this->SubAlg("FermiMover"));
  fVertexGenerator = dynamic_cast<const EventRecordVisitorI *> (this->SubAlg("VertexGenerator"));

}
//____________________________________________________________________________

