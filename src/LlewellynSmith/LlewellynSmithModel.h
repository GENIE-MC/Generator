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

class LlewellynSmithModel : public QELFormFactorsModelI {

public:

  virtual ~LlewellynSmithModel();

  //-- QELFormFactorModelI interface implementation

  virtual double F1V     (const Interaction * interaction) const;
  virtual double xiF2V   (const Interaction * interaction) const;
  virtual double FA      (const Interaction * interaction) const;
  virtual double Fp      (const Interaction * interaction) const;

protected:

  virtual double q2_4Mnucl2 (const Interaction * interaction) const;
  virtual double F1N        (const Interaction * interaction) const;
  virtual double munF2N     (const Interaction * interaction) const;
  virtual double GNE        (const Interaction * interaction) const;
  virtual double GNE0       (const Interaction * interaction) const;
  virtual double GNM        (const Interaction * interaction) const;
  virtual double GNM0       (const Interaction * interaction) const;
  virtual double GVE        (double q2) const;
  virtual double GVM        (double q2) const;

  LlewellynSmithModel();
  LlewellynSmithModel(const char * param_set);

//ClassDef(LlewellynSmithModel, 0)
};

}       // genie namespace

#endif  // _LLEWELLYN_SMITH_MODEL_H_

