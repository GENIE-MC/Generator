//____________________________________________________________________________
/*!

\class   genie::WeightCalculator

\brief

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created October 22, 2005

*/
//____________________________________________________________________________

#ifndef _WEIGHT_CALCULATOR_H_
#define _WEIGHT_CALCULATOR_H_

#include "ReWeight/MCModel.h"

namespace genie {

class EventRecord;

class WeightCalculator {

public :

  WeightCalculator();
  ~WeightCalculator();

  void OldCrossSectionModel(const MCModel & model);
  void NewCrossSectionModel(const MCModel & model);

  double ReWeight(const EventRecord & event);

private :

  MCModel fOldMCModel;
  MCModel fNewMCModel;
};

}      // genie namespace

#endif // _WEIGHT_CALCULATOR_H_
