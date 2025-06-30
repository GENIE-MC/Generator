//____________________________________________________________________________
/*!

  \class    genie::EngelFormFactor

  \brief    Form Factor for BertuzzoDNuCOHXSec...

  \ref      J. Engel
            Phys.Lett. B264, 114 (1991)

  \author   Iker de Icaza <i.de-icaza-astiz \at sussex.ac.uk>
            University of Sussex

            Costas Andreopoulos <c.andreopoulos \at cern.ch>
            University of Liverpool

  \created  June 12, 2020

  \cpright  Copyright (c) 2003-2025, The GENIE Collaboration
            For the full text of the license visit http://copyright.genie-mc.org

*/
//____________________________________________________________________________

#ifndef _ENGEL_FORM_FACTOR_H_
#define _ENGEL_FORM_FACTOR_H_

#include "Framework/Algorithm/Algorithm.h"
#include "Framework/Interaction/Target.h"

namespace genie {

class EngelFormFactor : public Algorithm {

public:
  EngelFormFactor();
  EngelFormFactor(string config);
  virtual ~EngelFormFactor();

  double FormFactor(const double Q, const Target & target) const ;
  // The Q has to be in GeV
  // The returned FF is in natural units

  //-- override the Algorithm::Configure methods to load configuration
  //   data to private data members
  void Configure (const Registry & config);
  void Configure (string param_set);


private:

  void LoadConfig(void);

};

} // genie namespace
#endif  // _ENGEL_FORM_FACTOR_H_
