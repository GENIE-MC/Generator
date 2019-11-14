//____________________________________________________________________________
/*!

\class    genie::QELEventGeneratorSM

\brief    Generates values for the kinematic variables describing QEL neutrino
          interaction events for Smith-Moniz model.
          Is a concrete implementation of the EventRecordVisitorI interface.

\ref      [1] R.A.Smith and E.J.Moniz, Nuclear Physics  B43, (1972) 605-622 \n
          [2] K.S. Kuzmin, V.V. Lyubushkin, V.A.Naumov Eur. Phys. J. C54, (2008) 517-538

\author   Igor Kakorin <kakorin@jinr.ru>, Joint Institute for Nuclear Research \n
          adopted from  fortran code provided by
          Konstantin Kuzmin <kkuzmin@theor.jinr.ru>, \n
          Joint Institute for Nuclear Research,  Institute for Theoretical and Experimental Physics \n
          Vladimir Lyubushkin, \n
          Joint Institute for Nuclear Research \n
          Vadim Naumov <vnaumov@theor.jinr.ru>, \n
          Joint Institute for Nuclear Research  \n
          based on code of Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk> \n
          University of Liverpool & STFC Rutherford Appleton Lab


\created  May 05, 2017

\cpright  Copyright (c) 2003-2017, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________


#ifndef _QEL_EVENT_GENERATORSM_H_
#define _QEL_EVENT_GENERATORSM_H_

#include "Physics/Common/KineGeneratorWithCache.h"
#include "Framework/Conventions/Controls.h"
#include "Physics/QuasiElastic/XSection/SmithMonizUtils.h"


namespace genie {

class QELEventGeneratorSM: public KineGeneratorWithCache {

public :
  QELEventGeneratorSM();
  QELEventGeneratorSM(string config);
 ~QELEventGeneratorSM();

  // implement the EventRecordVisitorI interface
  void ProcessEventRecord(GHepRecord * event_rec) const;

  // overload the Algorithm::Configure() methods to load private data
  // members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:

  mutable SmithMonizUtils * sm_utils;

  void   LoadConfig     (void);
  double ComputeMaxXSec(const Interaction * in) const;
  void AddTargetNucleusRemnant (GHepRecord * evrec) const; ///< add a recoiled nucleus remnant

  double ComputeMaxXSec2 (const Interaction * in) const;
  double MaxXSec2        (GHepRecord * evrec) const;
  double FindMaxXSec2    (const Interaction * in) const;
  void   CacheMaxXSec2   (const Interaction * in, double xsec) const;
  CacheBranchFx * AccessCacheBranch2 (const Interaction * in) const;

  double ComputeMaxDiffv (const Interaction * in) const;
  double MaxDiffv        (GHepRecord * evrec) const;
  double FindMaxDiffv    (const Interaction * in) const;
  void   CacheMaxDiffv   (const Interaction * in, double xsec) const;
  CacheBranchFx * AccessCacheBranchDiffv (const Interaction * in) const;

  mutable KinePhaseSpace_t fkps;
  

  bool fGenerateNucleonInNucleus;           ///< generate struck nucleon in nucleus
  double fQ2Min;                            ///< Q2-threshold for seeking the second maximum

  double fSafetyFacor_nu;


}; // class definition

} // genie namespace

#endif // _QEL_EVENT_GENERATORSM_H_
