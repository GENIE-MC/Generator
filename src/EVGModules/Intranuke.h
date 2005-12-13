//____________________________________________________________________________
/*!

\class   genie::Intranuke

\brief   The INTRANUKE cascading MC for intranuclear rescattering.

         Is a concrete implementation of the EventRecordVisitorI interface.

\ref     R.Merenyi et al., Phys.Rev.D45 (1992)
         R.D.Ransome, Nucl.Phys.B 139 (2005)

         The original INTRANUKE cascade MC was developed (in fortran) for the
         NeuGEN MC by G.F.Pearce, R.Edgecock, W.A.Mann and H.Gallagher.

\author  Hugh Gallagher <gallag@minos.phy.tufts.edu>, Tufts University
         Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk> CCLRC, Rutherford Lab

\created September 20, 2005

*/
//____________________________________________________________________________

#ifndef _INTRANUKE_H_
#define _INTRANUKE_H_

#include "EVGCore/EventRecordVisitorI.h"
#include "EVGModules/IntranukeProc.h"

class TLorentzVector;
class TVector3;

namespace genie {

class GHepParticle;

class Intranuke : public EventRecordVisitorI {

public :
  Intranuke();
  Intranuke(string config);
  ~Intranuke();

  //-- implement the EventRecordVisitorI interface
  void ProcessEventRecord(GHepRecord * event_rec) const;

private:

  bool        CanRescatter   (const GHepParticle * p) const;
  TVector3    Hadronic3P     (GHepRecord * event) const;
  void        StepParticle   (GHepParticle * p, double step) const;
  void        StepParticleMF (GHepParticle * p, double step, int Z) const;
  double      MeanFreePath   (double K) const;
  bool        IsInNucleus    (const GHepParticle * p, double R0) const;
  INukeProc_t ParticleFate   (const GHepParticle * p) const;
};

}      // genie namespace

#endif // _INTRANUKE_H_
