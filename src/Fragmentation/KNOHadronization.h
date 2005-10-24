//____________________________________________________________________________
/*!

\class    genie::KNOHadronization

\brief    The KNO hadronization model.

          This hadronization scheme is similar to the one originally used
          in NeuGEN by G.Barr, G.F.Pearce, H.Gallagher. \n

          Is a concrete implementation of the HadronizationModelI interface.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  August 17, 2004

*/
//____________________________________________________________________________

#ifndef _KNO_HADRONIZATION_H_
#define _KNO_HADRONIZATION_H_

#include <vector>

#include "Decay/DecayModelI.h"
#include "Fragmentation/HadronizationModelI.h"
#include "Fragmentation/MultiplicityProbModelI.h"

using std::vector;

namespace genie {

class KNOHadronization : public HadronizationModelI {

public:

  KNOHadronization();
  KNOHadronization(string config);
  virtual ~KNOHadronization();

  //-- implement the HadronizationModelI interface

  void           Initialize   (void)                 const;
  TClonesArray * Hadronize    (const Interaction * ) const;

private:

  const MultiplicityProbModelI * MultiplicityProbabilityModel (void) const;
  const DecayModelI *            DecayModel                   (void) const;

  vector<int> * GenerateFSHadronCodes (int mult, int maxQ, double W) const;
  int           GenerateBaryonPdgCode (int mult, int maxQ)           const;
  int           HadronShowerCharge    (const Interaction * proc)     const;
  void          HandleDecays          (TClonesArray * particle_list) const;
};

}         // genie namespace

#endif    // _KNO_HADRONIZATION_H_

