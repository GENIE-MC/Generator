//____________________________________________________________________________
/*!

\class    genie::LlewellynSmithModel

\brief    Abstract Base Class:
          implements the QELFormFactorsModelI interface but can not be
          instantiated.

          Its sole purpose of existence is to transmit common implementation
          (related to the Llewellyn-Smith model for QEL vN scattering) to its
          concrete subclasses: LlewellynSmithModelCC, LlewellynSmithModelNC.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 03, 2004

*/
//____________________________________________________________________________

#ifndef _LLEWELLYN_SMITH_MODEL_H_
#define _LLEWELLYN_SMITH_MODEL_H_

#include "Base/QELFormFactorsModelI.h"

namespace genie {

class ELFormFactorsModelI;

class LlewellynSmithModel : public QELFormFactorsModelI {

public:

  virtual ~LlewellynSmithModel();

  //-- QELFormFactorModelI interface implementation
  virtual double F1V     (const Interaction * interaction) const;
  virtual double xiF2V   (const Interaction * interaction) const;
  virtual double FA      (const Interaction * interaction) const;
  virtual double Fp      (const Interaction * interaction) const;

  //-- overload the Algorithm::Configure() methods to load private data
  //   members from configuration options
  virtual void Configure(const Registry & config);
  virtual void Configure(string config);

protected:

  LlewellynSmithModel();
  LlewellynSmithModel(string name);
  LlewellynSmithModel(string name, string config);

  virtual void LoadSubAlg     (void);
  virtual void LoadConfigData (void);

  virtual double tau    (const Interaction * interaction) const;
  virtual double GVE    (const Interaction * interaction) const;
  virtual double GVM    (const Interaction * interaction) const;

  const ELFormFactorsModelI * fElFFModel;

  double fMa2;
  double fFA0;
};

}       // genie namespace

#endif  // _LLEWELLYN_SMITH_MODEL_H_

