//____________________________________________________________________________
/*!

\class    genie::PythiaHadronization

\brief    Provides access to the PYTHIA hadronization models.

          Is a concrete implementation of the HadronizationModelI interface.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  August 17, 2004

*/
//____________________________________________________________________________

#ifndef _PYTHIA_HADRONIZATION_H_
#define _PYTHIA_HADRONIZATION_H_

#include <TPythia6.h>

#include "Fragmentation/HadronizationModelI.h"

namespace genie {

class PythiaHadronization : public HadronizationModelI {

public:

  PythiaHadronization();
  PythiaHadronization(const char * param_set);
  virtual ~PythiaHadronization();

  //-- define PythiaHadronization interface

  void           Initialize   (void)                 const;
  TClonesArray * Hadronize    (const Interaction * ) const;


  //-- tmp - std interface violoating method for PYTHIA config

  TPythia6 * PYTHIA(void) const { return fPythia; }

private:

  TPythia6 * fPythia;
  
};

}         // genie namespace

#endif    // _PYTHIA_HADRONIZATION__H_

