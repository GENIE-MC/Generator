//____________________________________________________________________________
/*!

\class   genie::GHepVirtualList

\brief   A GHepVirtualList is a 'virtual' collection of GHepParticles.
         Is virtual because it does not own but only points to GHepParticles
         owned by the generated GHepRecord.
         Use it if in your event generation algorithm you need to define & use
         a GHepRecord subset (without duplicating the GHepParticle entries)
         All 'named' lists are managed by the GHepVirtualListFolder singleton
         and get cleared after the generation of each event is completed.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created July 16, 2005

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
