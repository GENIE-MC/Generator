//____________________________________________________________________________
/*
 Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         STFC, Rutherford Appleton Laboratory - July 16, 2005

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <cassert>
#include <string>

#include <TLorentzVector.h>

#include "GHEP/GHepParticle.h"
#include "GHEP/GHepVirtualList.h"
#include "GHEP/GHepStatus.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGUtils.h"

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
