//____________________________________________________________________________
/*!

\class    genie::PDGLibrary

\brief    Singleton class to load & serve a TDatabasePDG.

\author   Costas Andreopoulos <c.andreopoulos \at cern.ch>
          University of Liverpool

          Changes required to implement the GENIE Boosted Dark Matter module
          were installed by Josh Berger (Univ. of Wisconsin)

          Changes required to implement the GENIE Dark Neutrino module
          were installed by Iker de Icaza (Univ. of Sussex)

	  Changes required to implement the GENIE BeamHNL module
	  were installed by John Plows (Univ. of Oxford)

\created  May 06, 2004

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _PDG_LIBRARY_H_
#define _PDG_LIBRARY_H_

#include <TDatabasePDG.h>
#include <TParticlePDG.h>

namespace genie {

class PDGLibrary
{
public:

  static PDGLibrary * Instance(void);

  TDatabasePDG * DBase (void);
  TParticlePDG * Find  (int pdgc, bool must_exist = true );
  void           ReloadDBase (void);

  // Add dark matter and mediator with parameters from Boosted Dark Matter app configuration
  // Ideally, this code should be in the Dark Matter app, not here.
  // But presently there is no way to edit the PDGLibrary after it has been created.
  void AddDarkMatter  (double mass, double med_ratio);

private:

  PDGLibrary();
  PDGLibrary(const PDGLibrary & config_pool);
  virtual ~PDGLibrary();

  bool LoadDBase(void);
  bool AddDarkSector ();
  bool AddHNL  (void);

  static PDGLibrary * fInstance;
  TDatabasePDG      * fDatabasePDG;

  struct Cleaner {
      void DummyMethodAndSilentCompiler() { }
      ~Cleaner() {
         if (PDGLibrary::fInstance !=0) {
            delete PDGLibrary::fInstance;
            PDGLibrary::fInstance = 0;
         }
      }
  };
  friend struct Cleaner;
};

}      // genie namespace

#endif // _PDG_LIBRARY_H_
