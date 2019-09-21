//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Shivesh Mandalia <s.p.mandalia@qmul.ac.uk>
         Queen Mary University of London
*/
//____________________________________________________________________________
#include <TPrimary.h>

#include "Physics/Hadronization/Pythia8Singleton.h"

#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Messenger/Messenger.h"

using namespace genie;
using namespace genie::constants;

ClassImp(Pythia8Singleton)

Pythia8Singleton * Pythia8Singleton::fgInstance = 0;

//____________________________________________________________________________
Pythia8Singleton::Pythia8Singleton()
{
    // Constructor
    if (fgInstance) {
      // TODO raise assertion?
      LOG("Pythia8Singleton", pERROR) <<
          "Instance of Pythia8Singleton already exists";
      return;
    }

#ifdef __GENIE_PYTHIA8_ENABLED__
    fPythia = new Pythia8::Pythia();
#endif
}
//____________________________________________________________________________
Pythia8Singleton::~Pythia8Singleton()
{
    // Destructor
    if (fgInstance) {
        fgInstance->Delete();
        delete fgInstance;
        fgInstance = 0;
    }
#ifdef __GENIE_PYTHIA8_ENABLED__
    delete fPythia;
#endif
}
//____________________________________________________________________________
Pythia8Singleton* Pythia8Singleton::Instance() 
{
    return fgInstance ? fgInstance : (fgInstance = new Pythia8Singleton());
}

