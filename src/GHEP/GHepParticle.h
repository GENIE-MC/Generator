//____________________________________________________________________________
/*!

\class   genie::GHepParticle

\brief   STDHEP-like event record entry that can fit a particle or a nucleus.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created November 18, 2004

*/
//____________________________________________________________________________

#ifndef _GHEP_PARTICLE_H_
#define _GHEP_PARTICLE_H_

#include <string>
#include <iostream>

#include <TObject.h>
#include <TLorentzVector.h>

#include "GHEP/GHepStatus.h"

using std::string;
using std::ostream;

namespace genie {

class GHepParticle : public TObject {

public :

  GHepParticle();
  GHepParticle(const GHepParticle & particle);

  // TParticle-like constructors for compatibility
  GHepParticle(
      int pdg, GHepStatus_t status,
            int mother1, int mother2, int daughter1, int daughter2,
                        const TLorentzVector & p, const TLorentzVector & v);
  GHepParticle(
       int pdg, GHepStatus_t status,
         int mother1, int mother2, int daughter1, int daughter2,
                           double px, double py, double pz, double E,
                                    double x, double y, double z, double t);
  ~GHepParticle();

  //-- Basic properties
  bool          IsNucleus      (void) const { return  fIsNucleus;          }
  bool          IsParticle     (void) const { return !fIsNucleus;          }
  bool          IsFake         (void) const { return  fIsFake;             }
  int           PdgCode        (void) const { return  fPdgCode;            }
  GHepStatus_t  Status         (void) const { return  fStatus;             }
  int           FirstMother    (void) const { return  fFirstMother;        }
  int           LastMother     (void) const { return  fLastMother;         }
  int           FirstDaughter  (void) const { return  fFirstDaughter;      }
  int           LastDaughter   (void) const { return  fLastDaughter;       }
  bool          HasDaughters   (void) const { return (fFirstDaughter!=-1); }

  //-- Returns the momentum & vertex 4-vectors
  TLorentzVector * P4 (void) const { return  fP4; }
  TLorentzVector * X4 (void) const { return  fX4; }

  //-- Hand over clones of the momentum & vertex 4-vectors (+ their ownership)
  TLorentzVector * GetP4 (void) const;
  TLorentzVector * GetX4 (void) const;

  //-- Returns the momentum & vertex 4-vectors components
  double Px     (void) const { return (fP4) ? fP4->Px()     : 0; }
  double Py     (void) const { return (fP4) ? fP4->Py()     : 0; }
  double Pz     (void) const { return (fP4) ? fP4->Pz()     : 0; }
  double E      (void) const { return (fP4) ? fP4->Energy() : 0; }
  double Energy (void) const { return this->E();                 }
  double Vx     (void) const { return (fX4) ? fX4->X()      : 0; }
  double Vy     (void) const { return (fX4) ? fX4->Y()      : 0; }
  double Vz     (void) const { return (fX4) ? fX4->Z()      : 0; }
  double Vt     (void) const { return (fX4) ? fX4->T()      : 0; }

  string Name   (void) const; // that corresponds to the PDG code
  double Mass   (void) const; // that corresponds to the PDG code
  double Charge (void) const; // that corresponds to the PDG code
  double KinE   (bool mass_from_pdg = false) const; // kinetic energy (E-m)

  bool Compare(const GHepParticle * p) const;

  //-- On/Off "shellness" if mass from PDG != mass from 4-P
  bool IsOnMassShell  (void) const;
  bool IsOffMassShell (void) const;

  //-- Relevant if entry is nucleus / Else -1 / Decoded from PDG code
  int  Z (void) const;
  int  A (void) const;

  //-- Set pdg code / status / parent-children links
  void SetPdgCode        (int c);
  void SetStatus         (GHepStatus_t s) { fStatus        = s; }
  void SetFirstMother    (int m)          { fFirstMother   = m; }
  void SetLastMother     (int m)          { fLastMother    = m; }
  void SetFirstDaughter  (int d)          { fFirstDaughter = d; }
  void SetLastDaughter   (int d)          { fLastDaughter  = d; }

  //-- Set the momentum & position 4-vectors
  void SetMomentum (const TLorentzVector & p4);
  void SetPosition (const TLorentzVector & v4);
  void SetMomentum (double px, double py, double pz, double E);
  void SetPosition (double x,  double y,  double z,  double t);
  void SetEnergy   (double E );
  void SetPx       (double px);
  void SetPy       (double py);
  void SetPz       (double pz);

  void CleanUp (void);
  void Reset   (void);
  void Clear   (Option_t * option);
  void Copy    (const GHepParticle & particle);
  void Print   (ostream & stream) const;
  void Print   (Option_t * opt)   const;

  bool             operator == (const GHepParticle & p) const;
  GHepParticle &   operator =  (const GHepParticle & p);
  friend ostream & operator << (ostream & stream, const GHepParticle & p);

private:

  void Init(void);
  void AssertIsKnownParticle(void) const;

  bool CompareFamily   (const GHepParticle * p) const;
  bool CompareMomentum (const GHepParticle * p) const;

  bool             fIsNucleus;      ///< nucleus flag
  bool             fIsFake;         ///< fake particle flag (rootino etc)
  int              fPdgCode;        ///< particle PDG code
  GHepStatus_t     fStatus;         ///< particle status
  int              fFirstMother;    ///< first mother idx
  int              fLastMother;     ///< last mother idx
  int              fFirstDaughter;  ///< first daughter idx
  int              fLastDaughter;   ///< last daughter idx
  TLorentzVector * fP4;             ///< momentum 4-vector
  TLorentzVector * fX4;             ///< position 4-vector 

ClassDef(GHepParticle, 1)

};

}      // genie namespace

#endif // _GHEP_PARTICLE_H_

