//____________________________________________________________________________
/*!

\class   genie::NtpMCGHepEntry

\brief   MINOS-style Ntuple class corresponding to a GHepParticle information.
         The Ntuple class is treated as a C-struct with public member data of
         basic-only types so that the ntuple can be easily analyzed in bare
         ROOT sessions (without loading the GENIE libraries).

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created October 1, 2004

\cpright  Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _NTP_GHEP_PARTICLE_H_
#define _NTP_GHEP_PARTICLE_H_

#include <ostream>

#include <TObject.h>

using std::ostream;

namespace genie {

class GHepParticle;

class NtpMCGHepEntry : public TObject {

public :
  NtpMCGHepEntry();
  NtpMCGHepEntry(const NtpMCGHepEntry & particle);
  virtual ~NtpMCGHepEntry();

  void Init (void);
  void Copy (const GHepParticle & particle);
  void Copy (const NtpMCGHepEntry & particle);

  void PrintToStream(ostream & stream) const;

  friend ostream & operator<< (ostream& stream, const NtpMCGHepEntry & stdhep);

  // Ntuple is treated like a C-struct with public data members and
  // rule-breaking field data members not prefaced by "f" and mostly lowercase.

  int     idx;      ///< particle index in event record
  int     ist;      ///< GHepStatus_t
  int     pdgc;     ///< PDG Code
  int     mom[2];   ///< First/Last mother
  int     dtr[2];   ///< First/Last daughter
  float   mass;     ///< Mass (GeV)
  float   p4[4];    ///< particle 4-momentum px,py,pz,E
  float   v4[4];    ///< particle 4-position vx,vy,vz,t

  ClassDef(NtpMCGHepEntry, 1)
};

}      // genie namespace

#endif // _NTP_GHEP_PARTICLE_H_
