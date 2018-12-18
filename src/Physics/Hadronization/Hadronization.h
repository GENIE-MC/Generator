//____________________________________________________________________________
/*!

\class    genie::Hadronizer

\brief    Base class for Hadronization classes
          Implements common configuration, allowing users to perform
          hadronization on a quark/diquark pair.
          Is a concerete implementation of the EventRecordVisitorI interface.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab
          Shivesh Mandalia <s.p.mandalia \at qmul.ac.uk>
          Queen Mary University of London

\created  November 19, 2019

\cpright  Copyright (c) 2003-2018, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _HADRONIZATION_H_
#define _HADRONIZATION_H_

class TClonesArray;
class TH1D;

#include "Framework/EventGen/EventRecordVisitorI.h"
#include "Framework/ParticleData/PDGCodeList.h"
#include "Framework/GHEP/GHepStatus.h"

namespace genie {

class GHepParticle;
class Interaction;

class Hadronization : public EventRecordVisitorI {

public:
    virtual ~Hadronization();

    // Overload the Algorithm::Configure() methods to load private data
    // members from configuration options
    virtual void Configure(const Registry & config);
    virtual void Configure(string config);

protected:
    Hadronization();
    Hadronization(string name);
    Hadronization(string name, string config);

    void  LoadConfigRemoveMe (void);

    virtual double         Wmin               (void) const;
    virtual double         Weight             (void) const;
    virtual TH1D *         CreateMultProbHist (double maxmult) const;
    virtual double         MaxMult            (const Interaction * i) const;
    virtual bool           AssertValidity     (const Interaction * i) const;
    virtual PDGCodeList *  SelectParticles    (const Interaction * i) const;
    virtual TClonesArray * Hadronize          (const Interaction * i) const;
    virtual void ApplyRijk (const Interaction * i, bool norm, TH1D * mp) const;
    virtual TH1D * MultiplicityProb (
            const Interaction * i, Option_t* opt="") const;

    bool fGenerateWeighted; ///< generate weighted or unweighted?

    const EventRecordVisitorI * fDecayer;

    // nuegen parameters
    double   fWcut;      ///< Rijk applied for W<Wcut (see DIS/RES join scheme)
    double   fRvpCCm2;   ///< Rijk: vp,  CC, multiplicity = 2
    double   fRvpCCm3;   ///< Rijk: vp,  CC, multiplicity = 3
    double   fRvpNCm2;   ///< Rijk: vp,  NC, multiplicity = 2
    double   fRvpNCm3;   ///< Rijk: vp,  NC, multiplicity = 3
    double   fRvnCCm2;   ///< Rijk: vn,  CC, multiplicity = 2
    double   fRvnCCm3;   ///< Rijk: vn,  CC, multiplicity = 3
    double   fRvnNCm2;   ///< Rijk: vn,  NC, multiplicity = 2
    double   fRvnNCm3;   ///< Rijk: vn,  NC, multiplicity = 3
    double   fRvbpCCm2;  ///< Rijk: vbp, CC, multiplicity = 2
    double   fRvbpCCm3;  ///< Rijk: vbp, CC, multiplicity = 3
    double   fRvbpNCm2;  ///< Rijk: vbp, NC, multiplicity = 2
    double   fRvbpNCm3;  ///< Rijk: vbp, NC, multiplicity = 3
    double   fRvbnCCm2;  ///< Rijk: vbn, CC, multiplicity = 2
    double   fRvbnCCm3;  ///< Rijk: vbn, CC, multiplicity = 3
    double   fRvbnNCm2;  ///< Rijk: vbn, NC, multiplicity = 2
    double   fRvbnNCm3;  ///< Rijk: vbn, NC, multiplicity = 3
};

}      // genie namespace
#endif // _HADRONIZATION_H_
