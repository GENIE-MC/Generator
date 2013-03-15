//____________________________________________________________________________
/*!

\class    genie::LwlynSmithFF

\brief    Abstract Base Class:
          implements the QELFormFactorsModelI interface but can not be
          instantiated.

          Its sole purpose of existence is to transmit common implementation
          (related to the Llewellyn-Smith model for QEL vN scattering) to its
          concrete subclasses: LwlynSmithFFCC, LwlynSmithFFNC.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  May 03, 2004

\cpright  Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _LLEWELLYN_SMITH_FORM_FACTOR_MODEL_H_
#define _LLEWELLYN_SMITH_FORM_FACTOR_MODEL_H_

#include "Base/QELFormFactorsModelI.h"
#include "ElFF/ELFormFactors.h"

namespace genie {

class ELFormFactorsModelI;

class LwlynSmithFF : public QELFormFactorsModelI {

public:

  virtual ~LwlynSmithFF();

  // QELFormFactorModelI interface implementation
  virtual double F1V     (const Interaction * interaction) const;
  virtual double xiF2V   (const Interaction * interaction) const;
  virtual double FA      (const Interaction * interaction) const;
  virtual double Fp      (const Interaction * interaction) const;

  // Overload the Algorithm::Configure() methods to load private data
  // members from configuration options
  virtual void Configure(const Registry & config);
  virtual void Configure(string config);

protected:

  LwlynSmithFF();
  LwlynSmithFF(string name);
  LwlynSmithFF(string name, string config);

  virtual void LoadConfig (void);

  virtual double tau    (const Interaction * interaction) const;
  virtual double GVE    (const Interaction * interaction) const;
  virtual double GVM    (const Interaction * interaction) const;

  const ELFormFactorsModelI * fElFFModel;

  mutable ELFormFactors fELFF;

  double fMa;  ///< axial mass
  double fMa2;
  double fFA0; ///< Fa(q2=0)
  double fMuP;
  double fMuN;
  double fSin28w;
};

}       // genie namespace

#endif  

