/*!___________________________________________________________________________

\class    genie::RSHelicityAmplModelI

\brief    Pure abstract base class. Defines the RSHelicityAmplModelI interface.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  July 10, 2004

____________________________________________________________________________*/

#ifndef _REIN_SEGHAL_HELICITY_AMPL_MODEL_I_H_
#define _REIN_SEGHAL_HELICITY_AMPL_MODEL_I_H_

#include "Algorithm/Algorithm.h"
#include "Interaction/Interaction.h"
#include "ReinSeghal/FKR.h"

namespace genie {

class RSHelicityAmplModelI : public Algorithm
{
public:

  virtual ~RSHelicityAmplModelI();

  //-- define the RSHelicityAmplModelI interface

  virtual double AmpMinus1 (const Interaction * interaction, const FKR & fkr) const = 0;
  virtual double AmpPlus1  (const Interaction * interaction, const FKR & fkr) const = 0;
  virtual double AmpMinus3 (const Interaction * interaction, const FKR & fkr) const = 0;
  virtual double AmpPlus3  (const Interaction * interaction, const FKR & fkr) const = 0;
  virtual double Amp0Minus (const Interaction * interaction, const FKR & fkr) const = 0;
  virtual double Amp0Plus  (const Interaction * interaction, const FKR & fkr) const = 0;

protected:

  RSHelicityAmplModelI();
  RSHelicityAmplModelI(const char * param_set);
};

}        // namespace

#endif   // _REIN_SEGHAL_HELICITY_AMPL_MODEL_I_H_



