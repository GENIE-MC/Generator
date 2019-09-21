//____________________________________________________________________________
/*!

\class    Pythia8Singleton

\brief    Provides access to the PYTHIA8 instance.

\author   Shivesh Mandalia <s.p.mandalia@qmul.ac.uk>
          Queen Mary University of London

\created  September 21, 2019

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________
#ifndef _PYTHIA8_SINGLETON_H_
#define _PYTHIA8_SINGLETON_H_

// Avoid the inclusion of dlfcn.h by Pythia.h that CINT is not able to process
#ifdef __CINT__
#define _DLFCN_H_
#define _DLFCN_H
#endif

#include <TObject.h>
#include <TPrimary.h>

#include "Framework/Conventions/GBuild.h"

#ifdef __GENIE_PYTHIA8_ENABLED__
#include "Pythia8/Pythia.h"

class Pythia;
#endif

class Pythia8Singleton : public TObject{

public:
    Pythia8Singleton();
    virtual ~Pythia8Singleton();

    static Pythia8Singleton * Instance ();
#ifdef __GENIE_PYTHIA8_ENABLED__
    Pythia8::Pythia * Pythia8  () {return fPythia;}
#endif

private:
    static Pythia8Singleton * fgInstance;  ///< singleton instance
#ifdef __GENIE_PYTHIA8_ENABLED__
    Pythia8::Pythia * fPythia;    ///< PYTHIA8 instance
#endif

ClassDef(Pythia8Singleton,1)
};

#endif    // _PYTHIA8_SINGLETON__H_
