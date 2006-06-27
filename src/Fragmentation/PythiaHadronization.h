//____________________________________________________________________________
/*!

\class    genie::PythiaHadronization

\brief    Provides access to the PYTHIA hadronization models. \n
          Is a concrete implementation of the HadronizationModelI interface.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  August 17, 2004

\cpright  Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE.
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
  PythiaHadronization(string config);
  virtual ~PythiaHadronization();

  //! define PythiaHadronization interface
  void           Initialize   (void)                 const;
  TClonesArray * Hadronize    (const Interaction * ) const;
  double         Weight       (void)                 const;

  //! overload the Algorithm::Configure() methods to load private data
  //! members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:

  void LoadConfig (void);

  TPythia6 * fPythia; ///< PYTHIA6 wrapper class

  //! configuration parameters
  double fSSBarSuppression;   ///< ssbar suppression
  double fGaussianPt2;        ///< gaussian pt2 distribution width
  double fNonGaussianPt2Tail; ///< non gaussian pt2 tail parameterization
  double fRemainingECutoff;   ///< remaining E cutoff for stopping fragmentation
};

}         // genie namespace

#endif    // _PYTHIA_HADRONIZATION__H_

