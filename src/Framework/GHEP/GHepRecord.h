//____________________________________________________________________________
/*!

\class    genie::GHepRecord

\brief    GENIE's GHEP MC event record.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  October 1, 2004

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _GHEP_RECORD_H_
#define _GHEP_RECORD_H_

#include <ostream>
#include <vector>

#include <TClonesArray.h>
#include <TBits.h>

#include "Framework/Conventions/GMode.h"
#include "Framework/Conventions/KinePhaseSpace.h"
#include "Framework/Interaction/Interaction.h" 
#include "Framework/GHEP/GHepStatus.h"

class TRootIOCtor;
class TLorentzVector;

using std::ostream;
using std::vector;

namespace genie {

class GHepRecord;
class GHepParticle;

ostream & operator << (ostream & stream, const GHepRecord & event);

class GHepRecord : public TClonesArray {

public :
  using TClonesArray::Print; // suppress clang 'hides overloaded virtual function [-Woverloaded-virtual]' warnings
  using TClonesArray::Copy;

  GHepRecord();
  GHepRecord(int size);
  GHepRecord(const GHepRecord & record);
  GHepRecord(TRootIOCtor*);
  virtual ~GHepRecord();

  // Methods to attach / get summary information

  virtual Interaction * Summary       (void) const;
  virtual void          AttachSummary (Interaction * interaction);

  // Provide a simplified wrapper of the 'new with placement'
  // TClonesArray object insertion method
  // ALWAYS use these methods to insert new particles as they check
  // for the compactness of the daughter lists.
  // Note that the record might be automatically re-arranged as the
  // result of your GHepParticle insertion

  virtual void AddParticle (const GHepParticle & p);
  virtual void AddParticle (int pdg, GHepStatus_t ist,
                     int mom1, int mom2, int dau1, int dau2,
                        const TLorentzVector & p, const TLorentzVector & v);
  virtual void AddParticle (int pdg, GHepStatus_t ist,
                     int mom1, int mom2, int dau1, int dau2,
                           double px, double py, double pz, double E,
                                    double x, double y, double z, double t);

  // Methods to search the GHEP record

  virtual GHepParticle * Particle     (int position) const;
  virtual GHepParticle * FindParticle (int pdg, GHepStatus_t ist, int start) const;

  virtual int ParticlePosition (int pdg, GHepStatus_t i, int start=0) const;
  virtual int ParticlePosition (GHepParticle * particle, int start=0) const;

  virtual vector<int> * GetStableDescendants(int position) const;

  // Return the mode (lepton+nucleon/nucleus, hadron+nucleon/nucleus, nucleon
  // decay etc...) by looking at the event entries

  GEvGenMode_t EventGenerationMode(void) const;

  // Easy access methods for the most frequently used GHEP entries

  virtual GHepParticle * Probe                            (void) const;
  virtual GHepParticle * TargetNucleus                    (void) const;
  virtual GHepParticle * RemnantNucleus                   (void) const;
  virtual GHepParticle * HitNucleon                       (void) const;
  virtual GHepParticle * HitElectron                      (void) const;
  virtual GHepParticle * FinalStatePrimaryLepton          (void) const;
  virtual GHepParticle * FinalStateHadronicSystem         (void) const;
  virtual int            ProbePosition                    (void) const;
  virtual int            TargetNucleusPosition            (void) const;
  virtual int            RemnantNucleusPosition           (void) const;
  virtual int            HitNucleonPosition               (void) const;
  virtual int            HitElectronPosition              (void) const;
  virtual int            FinalStatePrimaryLeptonPosition  (void) const;
  virtual int            FinalStateHadronicSystemPosition (void) const; 

  // Number of GHepParticle occurences in GHEP

  virtual unsigned int NEntries (int pdg, GHepStatus_t ist, int start=0) const;
  virtual unsigned int NEntries (int pdg, int start=0) const;

  // Methods to switch on/off and ask for event record flags

  virtual TBits * EventFlags   (void) const { return fEventFlags; }
  virtual TBits * EventMask    (void) const { return fEventMask; }
  virtual bool    IsUnphysical (void) const { return (fEventFlags->CountBits()>0); }
  virtual bool    Accept       (void) const;

  // Methods to set/get the event weight and cross sections

  virtual double Weight         (void) const  { return fWeight;   }
  virtual double Probability    (void) const  { return fProb;     }
  virtual double XSec           (void) const  { return fXSec;     }
  virtual double DiffXSec       (void) const  { return fDiffXSec; }
  virtual KinePhaseSpace_t DiffXSecVars  (void) const  { return fDiffXSecPhSp; }

  virtual void   SetWeight      (double wght) { fWeight   = (wght>0) ? wght : 0.; }
  virtual void   SetProbability (double prob) { fProb     = (prob>0) ? prob : 0.; }
  virtual void   SetXSec        (double xsec) { fXSec     = (xsec>0) ? xsec : 0.; }
  virtual void   SetDiffXSec    (double xsec, KinePhaseSpace_t ps) 
  { fDiffXSecPhSp = ps; 
    fDiffXSec = (xsec>0) ? xsec : 0.; 
  }

  // Set/get event vertex in detector coordinate system

  virtual TLorentzVector * Vertex (void) const { return fVtx; }

  virtual void SetVertex (double x, double y, double z, double t);
  virtual void SetVertex (const TLorentzVector & vtx);

  // Common event record operations

  virtual void Copy        (const GHepRecord & record);
  virtual void Clear       (Option_t * opt="");
  virtual void ResetRecord (void);
  virtual void CompactifyDaughterLists     (void);
  virtual void RemoveIntermediateParticles (void);

  // Set mask
  void SetUnphysEventMask(const TBits & mask);

  // Set/Get print level
  static void SetPrintLevel(int print_level);
  static int  GetPrintLevel();

  // Methods & operators to print the record

  void Print (ostream & stream) const;
  friend ostream & operator << (ostream & stream, const GHepRecord & event);

protected:

  // Attached interaction
  Interaction * fInteraction; ///< attached summary information

  // Vertex position
  TLorentzVector * fVtx;  ///< vertex in the detector coordinate system

  // Flags (and user-specified mask) for the generated event
  TBits * fEventFlags;    ///< event flags indicating various pathologies or an unphysical event
  TBits * fEventMask;     ///< an input bit-field mask allowing one to ignore bits set in fEventFlags

  // Event weight, probability and cross-sections
  double           fWeight;         ///< event weight
  double           fProb;           ///< event probability (for given flux neutrino and density-weighted path-length for target element)
  double           fXSec;           ///< cross section for selected event
  double           fDiffXSec;       ///< differential cross section for selected event kinematics
  KinePhaseSpace_t fDiffXSecPhSp;   ///< specifies which differential cross-section (dsig/dQ2, dsig/dQ2dW, dsig/dxdy,...)

  // Utility methods
  void InitRecord  (void);
  void CleanRecord (void);

  // Methods used by the daughter list compactifier
  virtual void UpdateDaughterLists    (void);
  virtual bool HasCompactDaughterList (int pos);
  virtual void SwapParticles          (int i, int j);
  virtual void FinalizeDaughterLists  (void);
  virtual int  FirstNonInitStateEntry (void);

  //
  static int fPrintLevel; //! print-level flag, see GHepRecord::Print()

private:

ClassDef(GHepRecord, 2)

};

}      // genie namespace

#endif // _GHEP_RECORD_H_
