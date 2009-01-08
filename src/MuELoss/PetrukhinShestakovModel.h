//____________________________________________________________________________
/*!

\class    genie::mueloss::PetrukhinShestakovModel

\brief    Bethe-Heitler, Petrukhin-Shestakov model for the energy loss of muons
          due to bremsstrahlung.
          Concrete implementation of the MuELossI interface.

\ref      W.Lohmann, R.Kopp and R.Voss,
          Energy Loss of Muons in the Energy Range 1-10000 GeV, CERN 85-03

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  December 10, 2003

*/
//____________________________________________________________________________

#ifndef _PETRUKHIN_SHESTAKOV_MODEL_H_
#define _PETRUKHIN_SHESTAKOV_MODEL_H_

#include "MuELoss/MuELossI.h"
#include "Numerical/GSFunc.h"

namespace genie {

class IntegratorI;

namespace mueloss {

class PetrukhinShestakovModel : public MuELossI
{
public:
  PetrukhinShestakovModel();
  PetrukhinShestakovModel(string config);
  virtual ~PetrukhinShestakovModel();

  //! implement the MuELossI interface
  double        dE_dx   (double E, MuELMaterial_t material) const;
  MuELProcess_t Process (void) const { return eMupBremsstrahlung; }

  //! overload the Algorithm::Configure() methods to load private data
  //!  members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:
  void LoadConfig (void);
  const IntegratorI * fIntegrator;
};

} // mueloss namespace
} // genie   namespace

//____________________________________________________________________________
/*!
\class    genie::mueloss::PetrukhinShestakovIntegrand

\brief    Auxiliary scalar function for the internal integration in Petrukhin
          Shestakov model

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  December 10, 2003
*/
//____________________________________________________________________________

namespace genie   {
namespace mueloss {

class PetrukhinShestakovIntegrand : public GSFunc
{
public:
  PetrukhinShestakovIntegrand(double E, double Z);
  ~PetrukhinShestakovIntegrand();

  double operator () (const vector<double> & x);

private:
  double fE;
  double fZ;
};

} // mueloss namespace
} // genie   namespace

#endif // _PETRUKHIN_SHESTAKOV_MODEL_H_
