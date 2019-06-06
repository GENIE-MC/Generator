//____________________________________________________________________________
/*!

\class    genie::GHepVirtualList

\brief    A GHepVirtualList is a 'virtual' collection of GHepParticles.
          Is virtual because it does not own but only points to GHepParticles
          owned by the generated GHepRecord. 
          Use it if in your event generation algorithm you need to define & use
          a GHepRecord subset (without duplicating the GHepParticle entries)
          All 'named' lists are managed by the GHepVirtualListFolder singleton
          and get cleared after the generation of each event is completed.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  July 16, 2005

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
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
