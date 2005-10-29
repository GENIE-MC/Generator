//____________________________________________________________________________
/*!

\class   genie::GHepRecord

\brief   Generated Event Record: STDHEP-like record and Summary Information.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created October 1, 2004

*/
//____________________________________________________________________________

#ifndef _GHEP_RECORD_H_
#define _GHEP_RECORD_H_

#include <ostream>

#include <TClonesArray.h>

#include "Interaction/Interaction.h"
#include "GHEP/GHepStatus.h"

class TLorentzVector;

using std::ostream;

namespace genie {

class GHepParticle;

class GHepRecord : public TClonesArray {

public :

  GHepRecord();
  GHepRecord(int size);
  GHepRecord(const GHepRecord & record);
  ~GHepRecord();

  //-- methods to attach / get summary information

  virtual Interaction * GetInteraction    (void) const;
  virtual void          AttachInteraction (Interaction * interaction);

  //-- operations on the record

  virtual void Copy                    (const GHepRecord & record);
  virtual void ShiftVertex             (const TLorentzVector & vec4);
  virtual void ResetRecord             (void);
  virtual void CompactifyDaughterLists (void);

  //-- provide a simplified wrapper of the 'new with placement'
  //   TClonesArray object insertion method

  //   ALWAYS use these methods to insert new particles as they check
  //   for the compactness of the daughter lists.
  //   Note that the record might be automatically re-arranged as the
  //   result of your GHepParticle insertion

  virtual void AddParticle (const GHepParticle & p);
  virtual void AddParticle (int pdg, GHepStatus_t ist,
                     int mom1, int mom2, int dau1, int dau2,
                        const TLorentzVector & p, const TLorentzVector & v);
  virtual void AddParticle (int pdg, GHepStatus_t ist,
                     int mom1, int mom2, int dau1, int dau2,
                           double px, double py, double pz, double E,
                                    double x, double y, double z, double t);

  //-- methods to search the GHEP (STDHEP-like) record

  virtual GHepParticle * GetParticle    (int position) const;
  virtual GHepParticle * FindParticle   (int pdg, GHepStatus_t ist, int start) const;

  virtual int ParticlePosition(int pdg, GHepStatus_t ist,  int start=0) const;
  virtual int ParticlePosition(GHepParticle * particle, int start=0) const;

  virtual unsigned int NEntries (int pdg, GHepStatus_t ist, int start=0) const;
  virtual unsigned int NEntries (int pdg, int start=0) const;

  //-- methods to switch on/off and ask for event record flags

  virtual void SwitchIsPauliBlocked (bool on_off);
  virtual void SwitchIsBelowThrNRF  (bool on_off);
  virtual void SwitchGenericErrFlag (bool on_off);
  virtual bool IsPauliBlocked       (void) const { return fIsPauliBlocked; }
  virtual bool IsBelowThrNRF        (void) const { return fIsBelowThrNRF;  }
  virtual bool GenericErrFlag       (void) const { return fGenericErrFlag; }
  virtual bool IsUnphysical         (void) const;

  //-- methods to set / get the event weight and cross sections

  virtual double GetWeight   (void) const  { return fWeight;   }
  virtual double GetXSec     (void) const  { return fXSec;     }
  virtual double GetDiffXSec (void) const  { return fDiffXSec; }
  virtual void   SetWeight   (double wght) { fWeight   = (wght>0) ? wght : 0.; }
  virtual void   SetXSec     (double xsec) { fXSec     = (xsec>0) ? xsec : 0.; }
  virtual void   SetDiffXSec (double xsec) { fDiffXSec = (xsec>0) ? xsec : 0.; }

  //-- methods & operators to print the record

  void Print (ostream & stream) const;

  friend ostream & operator << (ostream & stream, const GHepRecord & event);

protected:

  // Summary information for the Initial State, Process Type & Kinematics
  Interaction * fInteraction;

  // Flags for the generated event
  bool fIsPauliBlocked;   ///< true for Pauli-blocked event
  bool fIsBelowThrNRF;    ///< true if it is below threshold in the nucleon rest frame
  bool fGenericErrFlag;   ///< true for etc problems

  // Misc info associated with the generated event
  double fWeight;         ///< event weight
  double fXSec;           ///< cross section for selected event
  double fDiffXSec;       ///< differential cross section for selected event kinematics

  // Utility methods
  void InitRecord  (void);
  void CleanRecord (void);

  // Methods used by the daughter list compactifier
  virtual void UpdateDaughterLists    (void);
  virtual bool HasCompactDaughterList (int pos);
  virtual void SwapParticles          (int i, int j);
  virtual void FinalizeDaughterLists  (void);
  virtual int  FirstNonInitStateEntry (void);

private:

ClassDef(GHepRecord, 1)

};

}      // genie namespace

#endif // _GHEP_RECORD_H_
