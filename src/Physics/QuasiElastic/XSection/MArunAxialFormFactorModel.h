//____________________________________________________________________________
/*!

\class    genie::MArunAxialFormFactorModel

\brief    Concrete implementation of the AxialFormFactorModelI interface.
          Computes the axial form factor using a running MA

\ref      1. I.D. Kakorin, K.S. Kuzmin, V.A. Naumov, "Running axial mass of the nucleon 
             as a phenomenological tool for calculating quasielastic neutrino–nucleus cross sections", 
             Eur.Phys.J.C 81 (2021) 1142 [arXiV: 2112.13745 [hep-ph]].
          2. I.D. Kakorin, K.S. Kuzmin, V.A. Naumov, "A unified empirical model for quasielastic 
             interactions of neutrino and antineutrino with nuclei", Phys.Part.Nucl.Lett. 17 (2020) 265-288.
          3. K.S. Kuzmin, V.A. Naumov, O.N. Petrova, "Quasielastic neutrino–nucleus interactions in the 
             empirical model of running axial mass of the nucleon", Phys.Part.Nucl. 48 (2017) 995-997.

\author   Hugh Gallagher <hugh.gallagher@tufts.edu>
          From code provided by:
          Igor Kakorin <idkakorin@gmail.com>
          Joint Institute for Nuclear Research, Dubna

\created  August 1, 2016
\updated  April 10, 2024

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org          
*/
//____________________________________________________________________________

#ifndef _MARUN_AXIAL_FORM_FACTOR_MODEL_H_
#define _MARUN_AXIAL_FORM_FACTOR_MODEL_H_

#include "Physics/QuasiElastic/XSection/AxialFormFactorModelI.h"

namespace genie {

class MArunAxialFormFactorModel : public AxialFormFactorModelI {

public:
  MArunAxialFormFactorModel();
  MArunAxialFormFactorModel(string config);
  virtual ~MArunAxialFormFactorModel();

  // implement the AxialFormFactorModelI interface
  double FA (const Interaction * interaction) const;

  // overload Algorithm's Configure()
  void   Configure  (const Registry & config);
  void   Configure  (string param_set);

private:

  void LoadConfig(void);

  double fMa;  ///< axial mass
  double fMa2;
  double fFA0; ///< FA(q2=0)
  double fE0; ///< E0 for calculating running axial mass: Ma*(1+E0/Enu)
};

}         // genie namespace

#endif    // _MARUN_AXIAL_FORM_FACTOR_MODEL_H_
