//____________________________________________________________________________
/*!

\class    genie::KNODistribution

\brief    Describes a KNO Distribution.
          Its data are loaded from its XML configuration file.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  August 21, 2004

*/
//____________________________________________________________________________

#ifndef _KNO_DISTRIBUTION_H_
#define _KNO_DISTRIBUTION_H_

#include <TSpline.h>

#include "Algorithm/Algorithm.h"

namespace genie {

class KNODistribution : public Algorithm {

public:

  KNODistribution();
  KNODistribution(const char * param_set);
  virtual ~KNODistribution();

  double Value      (double n_avn) const;

  //-- Overload Algorithm's Configure() to make sure that the configuration
  //   registry (which is read by its config XML file) is translated to the
  //   KNO spline, when the algorithm is served by the AlgFactory
  
  void   Configure  (const Registry & config);
  void   Configure  (string param_set);
  
public:

  void KNOFromXmlConfig2Spline(void);

  TSpline * fKNOSpline;
  double    fMaxScaledMultiplicity; // n/<n>
};

}         // genie namespace

#endif    // _KNO_DISTRIBUTION_H_

