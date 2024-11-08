//____________________________________________________________________________
/*
 Copyright (c) 2003-2024, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

*/
//____________________________________________________________________________

#include <sstream>
#include <iostream>

#include "Framework/Messenger/Messenger.h"
#include "Framework/Utils/Pythia8Singleton.h"

namespace genie {

//____________________________________________________________________________
//ostream & operator << (ostream & stream, const Pythia8Singleton & cache)
//{
//  cache.Print(stream);
//  return stream;
//}
//____________________________________________________________________________
Pythia8Singleton * Pythia8Singleton::fInstance = 0;
//____________________________________________________________________________
Pythia8Singleton::Pythia8Singleton()
{
  fInstance = 0;
  fPythia   = 0;

}
//____________________________________________________________________________
Pythia8Singleton::~Pythia8Singleton()
{
  if (fPythia) {
    delete fPythia;
  }
  fInstance = 0;
}
//____________________________________________________________________________
Pythia8Singleton * Pythia8Singleton::Instance()
{
  if (fInstance == 0) {
    static Pythia8Singleton::Cleaner cleaner;
    cleaner.DummyMethodAndSilentCompiler();

    fInstance = new Pythia8Singleton;

    // actually create the one to be held
    fInstance->fPythia = new Pythia8::Pythia();
  }
  return fInstance;
}
//____________________________________________________________________________
//void Pythia8Singleton::Print(ostream & stream) const
//{
//  stream << "\n [-] GENIE Pythia8Singleton:";
//  stream << "\n ";
//  stream << "\n";
//}
//___________________________________________________________________________

} // genie namespace
