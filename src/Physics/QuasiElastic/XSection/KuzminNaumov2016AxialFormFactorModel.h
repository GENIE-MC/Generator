//____________________________________________________________________________
/*!

\class    genie::KuzminNaumov2016AxialFormFactorModel

\brief    Concrete implementation of the AxialFormFactorModelI interface.
          Computes the axial form factor using a running MA 

\ref      Konstantin S. Kuzmin and Vadim A. Naumov. 
          Running axial-vector mass of the nucleon for a precise evaluation of the 
          quasielastic (anti)neutrino–nucleus cross sectionsent. 2016 (in preparation).

\author   Hugh Gallagher <hugh.gallagher \at tufts.edu>
          From code provided by: 
          Igor Kakorin <idkakorin \at gmail.com>
          Joint Institute for Nuclear Research, Dubna

\created  August 1, 2016

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _KUZMIN_NAUMOV_2016_AXIAL_FORM_FACTOR_MODEL_H_
#define _KUZMIN_NAUMOV_2016_AXIAL_FORM_FACTOR_MODEL_H_

#include "Physics/QuasiElastic/XSection/AxialFormFactorModelI.h"

namespace genie {

class KuzminNaumov2016AxialFormFactorModel : public AxialFormFactorModelI {

public:
  KuzminNaumov2016AxialFormFactorModel();
  KuzminNaumov2016AxialFormFactorModel(string config);
  virtual ~KuzminNaumov2016AxialFormFactorModel();

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

#endif    // _KUZMIN_NAUMOV_2016_AXIAL_FORM_FACTOR_MODEL_H_
