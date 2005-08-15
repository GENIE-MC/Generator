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

#ifndef _GHEP_VIRTUAL_LIST_H_
#define _GHEP_VIRTUAL_LIST_H_

#include <TClonesArray.h>

class TLorentzVector;

namespace genie {

class GHepParticle;

class GHepVirtualList : public TClonesArray {

public :

  GHepVirtualList();
  GHepVirtualList(int size);
  GHepVirtualList(const GHepVirtualList & vlist);
  ~GHepVirtualList();

private :

ClassDef(GHepVirtualList, 1)

};

}      // genie namespace

#endif // _GHEP_VIRTUAL_LIST_H_
