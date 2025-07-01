//____________________________________________________________________________
/*!

\class    genie::Pythia8Singleton

\brief    Manage a single instance of pythia8

\author   Robert Hatcher <rhatcher \at fnal.gov>
          Fermilab

\created  May 15, 2024

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _PYTHIA8SINGLETON_H_
#define _PYTHIA8SINGLETON_H_

#ifdef __GENIE_PYTHIA8_ENABLED__
#include "Pythia8/Pythia.h"
#endif // __GENIE_PYTHIA8_ENABLED__

namespace genie {

class Pythia8Singleton;

//ostream & operator << (ostream & stream, const Pythia8Singleton & pythia8_1);

class Pythia8Singleton
{
public:

  static Pythia8Singleton * Instance(void);

#ifdef __GENIE_PYTHIA8_ENABLED__
  Pythia8::Pythia * Pythia8() { return fPythia; }
#endif // __GENIE_PYTHIA8_ENABLED__

  //! print
  //void   Print (ostream & stream) const;
  //friend ostream & operator << (ostream & stream, const Pythia8Singleton & pythia8_1);

private:

  //! singleton instance
  static Pythia8Singleton * fInstance;

#ifdef __GENIE_PYTHIA8_ENABLED__
  mutable Pythia8::Pythia * fPythia; ///< actual Pythia8 instance
#endif // __GENIE_PYTHIA8_ENABLED__

  //! singleton class: constructors are private
  Pythia8Singleton();
  Pythia8Singleton(const Pythia8Singleton & cache);
  virtual ~Pythia8Singleton();

  //! proper de-allocation of the singleton object
  struct Cleaner {
      void DummyMethodAndSilentCompiler() { }
      ~Cleaner() {
         if (Pythia8Singleton::fInstance !=0) {
            delete Pythia8Singleton::fInstance;
            Pythia8Singleton::fInstance = 0;
         }
      }
  };
  friend struct Cleaner;
};

}      // genie namespace

#endif // _PYTHIA8SINGLETON_H_
