//____________________________________________________________________________
/*!

\class    genie::HadronizationModelBase

\brief    An abstract class. It avoids implementing the HadronizationModelI
          interface, leaving it for the concrete subclasses (KNO, Pythia,...).
          It propagates some common methods to the concrete subclasses.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  August 17, 2004

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _HADRONIZATION_MODEL_BASE_H_
#define _HADRONIZATION_MODEL_BASE_H_

#include "Physics/Hadronization/HadronizationModelI.h"

namespace genie {

class HadronizationModelBase : public HadronizationModelI {

public:

  //! Don't implement the HadronizationModelI interface
  //! Leave it for the concrete implementations (KNO, Pythia,...)

  virtual void           Initialize       (void)                                 const = 0;
  virtual TClonesArray * Hadronize        (const Interaction * )                 const = 0;
  virtual double         Weight           (void)                                 const = 0;
  virtual PDGCodeList *  SelectParticles  (const Interaction*)                   const = 0;
  virtual TH1D *         MultiplicityProb (const Interaction*, Option_t* opt="") const = 0;

protected:
  HadronizationModelBase();
  HadronizationModelBase(string name);
  HadronizationModelBase(string name, string config);
  virtual ~HadronizationModelBase();

  //! Various utility methods common to hadronization models

  double Wmin               (void) const;
  double MaxMult            (const Interaction * i) const;
  void   ApplyRijk          (const Interaction * i, bool norm, TH1D * mp) const;
  TH1D * CreateMultProbHist (double maxmult) const;

  //! configuration data common to all hadronizers
  double   fWcut;        ///< neugen's Rijk applied for W<Wcut (see DIS/RES join scheme)
  double   fRvpCCm2;     ///< neugen's Rijk: vp,  CC, multiplicity = 2
  double   fRvpCCm3;     ///< neugen's Rijk: vp,  CC, multiplicity = 3
  double   fRvpNCm2;     ///< neugen's Rijk: vp,  NC, multiplicity = 2
  double   fRvpNCm3;     ///< neugen's Rijk: vp,  NC, multiplicity = 3
  double   fRvnCCm2;     ///< neugen's Rijk: vn,  CC, multiplicity = 2
  double   fRvnCCm3;     ///< neugen's Rijk: vn,  CC, multiplicity = 3
  double   fRvnNCm2;     ///< neugen's Rijk: vn,  NC, multiplicity = 2
  double   fRvnNCm3;     ///< neugen's Rijk: vn,  NC, multiplicity = 3
  double   fRvbpCCm2;    ///< neugen's Rijk: vbp, CC, multiplicity = 2
  double   fRvbpCCm3;    ///< neugen's Rijk: vbp, CC, multiplicity = 3
  double   fRvbpNCm2;    ///< neugen's Rijk: vbp, NC, multiplicity = 2
  double   fRvbpNCm3;    ///< neugen's Rijk: vbp, NC, multiplicity = 3
  double   fRvbnCCm2;    ///< neugen's Rijk: vbn, CC, multiplicity = 2
  double   fRvbnCCm3;    ///< neugen's Rijk: vbn, CC, multiplicity = 3
  double   fRvbnNCm2;    ///< neugen's Rijk: vbn, NC, multiplicity = 2
  double   fRvbnNCm3;    ///< neugen's Rijk: vbn, NC, multiplicity = 3
};

}         // genie namespace

#endif    // _HADRONIZATION_MODEL_BASE_H_

