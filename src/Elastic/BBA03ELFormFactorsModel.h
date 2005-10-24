//____________________________________________________________________________
/*!

\class    genie::BBA03ELFormFactorsModel

\brief    Concrete implementation of the ELFormFactorsModelI interface.
          Computes elastic form factors using the BBA2003 parameterization.

\ref      H.Budd, NuINT-02 proceedings

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  October 19, 2005

*/
//____________________________________________________________________________

#ifndef _BBA2003_EL_FORM_FACTORS_MODEL_H_
#define _BBA2003_EL_FORM_FACTORS_MODEL_H_

#include "Base/ELFormFactorsModelI.h"

namespace genie {

typedef struct SBBA2003Fit
{
  double a2, a4, a6, a8, a10, a12;
}
BBA2003Fit_t;

class BBA03ELFormFactorsModel : public ELFormFactorsModelI {

public:
  BBA03ELFormFactorsModel();
  BBA03ELFormFactorsModel(string config);
  virtual ~BBA03ELFormFactorsModel();

  // implement the ELFormFactorsModelI interface
  double Gep (double q2) const;
  double Gmp (double q2) const;
  double Gen (double q2) const;
  double Gmn (double q2) const;

  // overload Algorithm's Configure() to load the BBA2003Fit_t
  // structs from the configuration Registry
  void   Configure  (const Registry & config);
  void   Configure  (string param_set);

private:

  // the actual BBA2003 inverse polynomial fit function
  double BBA03Fit(double q2, double g0, const BBA2003Fit_t & fp) const;

  // fill data members from the configuration Registry
  void LoadBBA2003Params(void);

  // BBA2003 model parameters.
  // These are filled in by the configuration Registry entries
  // or are set to defaults.
  BBA2003Fit_t fGep;   ///< BBA2003 fit coefficients for Gep
  BBA2003Fit_t fGmp;   ///< BBA2003 fit coefficients for Gmp
  BBA2003Fit_t fGmn;   ///< BBA2003 fit coefficients for Gmn
  double       fGenA;  ///< Krutov parameterization for Gen
  double       fGenB;  ///< Krutov parameterization for Gen
  double       fQ2Max; ///< Gep/Gmp assummed const for Q2 > Q2Max
  double       fMv2;   ///< Elactic vector mass ^ 2
};

}         // genie namespace

#endif    // _BBA2003_EL_FORM_FACTORS_MODEL_H_
