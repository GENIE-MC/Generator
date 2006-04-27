//____________________________________________________________________________
/*!

\class   genie::NtpMCSummary

\brief   MINOS-style Ntuple class to hold an MC event summary information.
         The Ntuple class is treated as a C-struct with public member data of
         basic-only types so that the ntuple can be easily analyzed in bare
         ROOT sessions (without loading the GENIE libraries).

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created October 1, 2004

\cpright  Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _NTP_MC_SUMMARY_H_
#define _NTP_MC_SUMMARY_H_

#include <ostream>

#include <TObject.h>

using std::ostream;

namespace genie {

class EventRecord;

class NtpMCSummary : public TObject {

public :

  NtpMCSummary();
  NtpMCSummary(const NtpMCSummary & mcs);
  virtual ~NtpMCSummary();

  void Init(void);

  void Copy (const NtpMCSummary & mcs);
  void Copy (const EventRecord & evrec);

  void PrintToStream(ostream & stream) const;

  friend ostream & operator<< (ostream& stream, const NtpMCSummary & summary);

  // Ntuple is treated like a C-struct with public data members and
  // rule-breaking field data members not prefaced by "f" and mostly lowercase
  // with in -some cases- absurdely abbreviated & cryptic names for fast typing
  // from the ROOT command line
  int          probe;     ///< Probe (v,e...) PDG code
  int          fsl;       ///< Final state primary lepton PDG code
  int          tgt;       ///< Nuclear target PDG code
  int          nucl;      ///< Struck nucleon PDG code
  int          iqrk;      ///< Struck quark PDG code
  int          fqrk;      ///< Outgoing quark PDG code
  int          res;       ///< Resonance PDG code
  int          ch;        ///< Charm hadron PDG code
  int          Z;         ///< Nuclear target: Z
  int          A;         ///< Nuclear target: A
  int          N;         ///< Nuclear target: N
  int          scat;      ///< Scattering type (QEL, DIS, RES, COH, ...)
  int          proc;      ///< Process type (E/M, CC, NC, ...)
  double       x;         ///< Bjorken x
  double       y;         ///< Inelasticity
  double       z;         ///< Fragmentation z
  double       Q2;        ///< Momentum transfer (>0 = -q2)
  double       W;         ///< Hadronic invariant mass
  double       emfrac;    ///< E/M fraction
  double       xsec;      ///< Cross section at this energy
  double       dxsec;     ///< Cross section for the selected kinematics
  double       wght;      ///< Event weight (if weighted events are produced)
  double       v[4];      ///< Interaction vertex (x,y,z,t)
  double       p4p[4];    ///< Probe 4-P (px,py,pz,E) in LAB
  double       p4nucl[4]; ///< Struck nucleon 4-P (px,py,py,E) in LAB
  double       p4fsl[4];  ///< Final state primary lepton 4-P (px,py,py,E) in LAB
  double       p4fsl2[4]; ///<
  double       p4fsh[4];  ///< Hadronic system 4-P (px,py,py,E) in LAB
  unsigned int np;        ///< Number of generated protons
  unsigned int nn;        ///< Number of generated neutrons
  unsigned int npi0;      ///< Number of generated pi^0
  unsigned int npip;      ///< Number of generated pi^+
  unsigned int npim;      ///< Number of generated pi^-
  unsigned int nK0;       ///< Number of generated K^0
  unsigned int nKp;       ///< Number of generated K^+
  unsigned int nKm;       ///< Number of generated K^-

private:

ClassDef(NtpMCSummary, 1)
};

}      // genie namespace

#endif // _NTP_MC_SUMMARY_H_
