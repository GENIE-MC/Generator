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

  virtual double tau    (const Interaction * interaction) const;
  virtual double GVE    (const Interaction * interaction) const;
  virtual double GVM    (const Interaction * interaction) const;

  LlewellynSmithModel();
  LlewellynSmithModel(const char * param_set);
};

}       // genie namespace

#endif  // _LLEWELLYN_SMITH_MODEL_H_

