//____________________________________________________________________________
/*!

\class   genie::GHepParticle

\brief   STDHEP-like event record entry that can fit a particle or a nucleus.

\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

\created November 18, 2004

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _GHEP_PARTICLE_H_
#define _GHEP_PARTICLE_H_

#include <string>
#include <iostream>

#include <TObject.h>
#include <TLorentzVector.h>

#include "Framework/GHEP/GHepStatus.h"

class TRootIOCtor;

using std::string;
using std::ostream;

namespace genie {

class GHepParticle;
ostream & operator << (ostream & stream, const GHepParticle & p);

class GHepParticle : public TObject {

public :
  using TObject::Copy; // suppress clang 'hides overloaded virtual function [-Woverloaded-virtual]' warnings
  using TObject::Compare;

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

  GHepParticle(TRootIOCtor*);
 ~GHepParticle();

  // Basic properties
  int           Pdg            (void) const { return  fPdgCode;            }
  GHepStatus_t  Status         (void) const { return  fStatus;             }
  int           RescatterCode  (void) const { return  fRescatterCode;      }
  int           FirstMother    (void) const { return  fFirstMother;        }
  int           LastMother     (void) const { return  fLastMother;         }
  int           FirstDaughter  (void) const { return  fFirstDaughter;      }
  int           LastDaughter   (void) const { return  fLastDaughter;       }
  bool          HasDaughters   (void) const { return (fFirstDaughter!=-1); }
  bool          IsBound        (void) const { return  fIsBound;            }

  string Name   (void) const; ///< Name that corresponds to the PDG code
  double Mass   (void) const; ///< Mass that corresponds to the PDG code
  double Charge (void) const; ///< Chrg that corresponds to the PDG code

  // Returns the momentum & position 4-vectors
  const TLorentzVector * P4 (void) const { return fP4; }
  const TLorentzVector * X4 (void) const { return fX4; }
  TLorentzVector * P4 (void) { return fP4; }
  TLorentzVector * X4 (void) { return fX4; }

  // Hand over clones of the momentum & position 4-vectors (+ their ownership)
  TLorentzVector * GetP4 (void) const;
  TLorentzVector * GetX4 (void) const;

  // Returns the momentum & position 4-vectors components
  double Px     (void) const { return (fP4) ? fP4->Px()     : 0; } ///< Get Px
  double Py     (void) const { return (fP4) ? fP4->Py()     : 0; } ///< Get Py
  double Pz     (void) const { return (fP4) ? fP4->Pz()     : 0; } ///< Get Pz 
  double E      (void) const { return (fP4) ? fP4->Energy() : 0; } ///< Get energy
  double Energy (void) const { return this->E();                 } ///< Get energy
  double KinE   (bool mass_from_pdg = false) const;                ///< Get kinetic energy
  double Vx     (void) const { return (fX4) ? fX4->X()      : 0; } ///< Get production x
  double Vy     (void) const { return (fX4) ? fX4->Y()      : 0; } ///< Get production y
  double Vz     (void) const { return (fX4) ? fX4->Z()      : 0; } ///< Get production z
  double Vt     (void) const { return (fX4) ? fX4->T()      : 0; } ///< Get production time

  // Return removal energy /set only for bound nucleons/
  double RemovalEnergy (void) const { return fRemovalEnergy; } ///< Get removal energy 

  // Compare with another particle
  bool Compare            (const GHepParticle * p) const;
  bool ComparePdgCodes    (const GHepParticle * p) const;
  bool CompareStatusCodes (const GHepParticle * p) const;
  bool CompareFamily      (const GHepParticle * p) const;
  bool CompareMomentum    (const GHepParticle * p) const;

  // On/Off "shellness" if mass from PDG != mass from 4-P
  bool IsOnMassShell  (void) const;
  bool IsOffMassShell (void) const;

  // Relevant if GHEP entry is a nucleus, else=-1 / Decoded from PDG code
  int  Z (void) const;
  int  A (void) const;

  // Get the polarization. Most likely it is only the f/s primary lepton
  // for which this is usefull and might be set during event generation
  double PolzPolarAngle   (void) const { return fPolzTheta; }
  double PolzAzimuthAngle (void) const { return fPolzPhi; }
  bool   PolzIsSet        (void) const;
  void   GetPolarization  (TVector3 & polz);

  // Set pdg code and status codes
  void SetPdgCode  (int c);
  void SetStatus   (GHepStatus_t s) { fStatus = s; }

  // Set the rescattering code
  void SetRescatterCode(int code) { fRescatterCode = code; }

  // Set the mother/daughter links
  void SetFirstMother    (int m)          { fFirstMother   = m; }
  void SetLastMother     (int m)          { fLastMother    = m; }
  void SetFirstDaughter  (int d)          { fFirstDaughter = d; }
  void SetLastDaughter   (int d)          { fLastDaughter  = d; }

  // Set the momentum & position 4-vectors
  void SetMomentum (const TLorentzVector & p4);
  void SetPosition (const TLorentzVector & v4);
  void SetMomentum (double px, double py, double pz, double E);
  void SetPosition (double x,  double y,  double z,  double t);
  void SetEnergy   (double E );
  void SetPx       (double px);
  void SetPy       (double py);
  void SetPz       (double pz);

  // Set the polarization angles
  void SetPolarization(double theta, double phi);
  void SetPolarization(const TVector3 & polz);

  // Set the bould flag & removal energy (bound flag set automatically 
  // if a positive removal energy is set)
  void SetBound         (bool bound);
  void SetRemovalEnergy (double Erm);

  // Clean-up, reset, copy, print,...
  void CleanUp (void);
  void Reset   (void);
  void Clear   (Option_t * option);
  void Copy    (const GHepParticle & particle);
  void Print   (ostream & stream) const;
  void Print   (Option_t * opt)   const;

  // Overloaded operators
  bool             operator == (const GHepParticle & p) const;
  GHepParticle &   operator =  (const GHepParticle & p);
  friend ostream & operator << (ostream & stream, const GHepParticle & p);

private:

  void Init(void);
  void AssertIsKnownParticle(void) const;

  int              fPdgCode;        ///< particle PDG code
  GHepStatus_t     fStatus;         ///< particle status
  int              fRescatterCode;  ///< rescattering code
  int              fFirstMother;    ///< first mother idx
  int              fLastMother;     ///< last mother idx
  int              fFirstDaughter;  ///< first daughter idx
  int              fLastDaughter;   ///< last daughter idx
  TLorentzVector * fP4;             ///< momentum 4-vector (GeV)
  TLorentzVector * fX4;             ///< position 4-vector (in the target nucleus coordinate system / x,y,z in fm / t=0)
  double           fPolzTheta;      ///< polar polarization angle (rad)
  double           fPolzPhi;        ///< azimuthal polarization angle (rad)
  double           fRemovalEnergy;  ///< removal energy for bound nucleons (GeV)
  bool             fIsBound;        ///< 'is it a bound particle?' flag

ClassDef(GHepParticle, 2)

};

}      // genie namespace

#endif // _GHEP_PARTICLE_H_

