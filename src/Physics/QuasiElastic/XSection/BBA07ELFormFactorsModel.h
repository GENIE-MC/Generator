//____________________________________________________________________________
/*!

\class    genie::BBA07ELFormFactorsModel

\brief    Computes elastic form factors using the BBA2007 parameterization.
          Concrete implementation of the ELFormFactorsModelI interface.

\ref      A.Bodek, R.Bradford, H.Budd and S.Avvakumov,
          Euro.Phys.J.C53 (2008);[arXiv:0708.1946 [hep-ex]]


\author   Igor Kakorin <kakorin@jinr.ru>
          Joint Institute for Nuclear Research \n

          adapted from  fortran code provided by:
          Konstantin Kuzmin <kkuzmin@theor.jinr.ru>, \n
          Joint Institute for Nuclear Research
          Institute for Theoretical and Experimental Physics \n

          Vladimir Lyubushkin, \n
          Joint Institute for Nuclear Research \n

          Vadim Naumov <vnaumov@theor.jinr.ru>, \n
          Joint Institute for Nuclear Research  \n

          based on code of:
          Costas Andreopoulos <c.andreopoulos \at cern.ch> \n
          University of Liverpool

\created  Dec 01, 2017

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org

*/
//____________________________________________________________________________

#ifndef _BBA2007_EL_FORM_FACTORS_MODEL_H_
#define _BBA2007_EL_FORM_FACTORS_MODEL_H_

#include "Physics/QuasiElastic/XSection/ELFormFactorsModelI.h"

namespace genie {

typedef struct SBBA2007Fit
{
  double a1, b1, b2, b3, p1, p2, p3, p4, p5, p6, p7;
}
BBA2007Fit_t;

class BBA07ELFormFactorsModel : public ELFormFactorsModelI {

public:
  BBA07ELFormFactorsModel();
  BBA07ELFormFactorsModel(string config);
  virtual ~BBA07ELFormFactorsModel();

  // implement the ELFormFactorsModelI interface
  double Gep (const Interaction * interaction) const;
  double Gmp (const Interaction * interaction) const;
  double Gen (const Interaction * interaction) const;
  double Gmn (const Interaction * interaction) const;

  // overload Algorithm's Configure() to load the BBA2007Fit_t
  // structs from the configuration Registry
  void   Configure  (const Registry & config);
  void   Configure  (string param_set);

private:

  // fill data members from the configuration Registry
  void LoadConfig(void);


  double AN (double x,double c1, double c2, double c3,double c4,double c5, double c6, double c7) const;



  // model parameters.
  BBA2007Fit_t fGep;   ///< BBA2007 fit coefficients for Gep
  BBA2007Fit_t fGen;   ///< BBA2007 fit coefficients for Gen
  BBA2007Fit_t fGmp;   ///< BBA2007 fit coefficients for Gmp
  BBA2007Fit_t fGmn;   ///< BBA2007 fit coefficients for Gmn
  double       fMuP;   ///< Anomalous proton magnetic moment
  double       fMuN;   ///< Anomalous neutron magnetic moment
};

}         // genie namespace

#endif    // _BBA2007_EL_FORM_FACTORS_MODEL_H_
