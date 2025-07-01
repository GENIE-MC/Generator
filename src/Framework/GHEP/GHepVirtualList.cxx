//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Costas Andreopoulos <c.andreopoulos \at cern.ch>
 University of Liverpool
*/
//____________________________________________________________________________

#include <cassert>
#include <string>

#include <TLorentzVector.h>

#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepVirtualList.h"
#include "Framework/GHEP/GHepStatus.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGUtils.h"

using std::string;
using namespace genie;

ClassImp(GHepVirtualList)

//___________________________________________________________________________
GHepVirtualList::GHepVirtualList() :
TClonesArray("genie::GHepParticle")
{
  this->SetOwner(false);
}
//___________________________________________________________________________
GHepVirtualList::GHepVirtualList(int size) :
TClonesArray("genie::GHepParticle", size)
{
  this->SetOwner(false);
}
//___________________________________________________________________________
GHepVirtualList::GHepVirtualList(const GHepVirtualList & vlist) :
TClonesArray("genie::GHepParticle", vlist.GetEntries())
{
  this->SetOwner(false);

  // clean up
  this->Clear();

  // adjust size
  this->Expand(vlist.GetEntries());

  // copy event record entries
  unsigned int ientry = 0;
  GHepParticle * p = 0;
  TIter ghepiter(&vlist);
  while ( (p = (GHepParticle *) ghepiter.Next()) ) (*this)[ientry++] = p;
}
//___________________________________________________________________________
GHepVirtualList::~GHepVirtualList()
{
  this->Clear();
}
//___________________________________________________________________________
