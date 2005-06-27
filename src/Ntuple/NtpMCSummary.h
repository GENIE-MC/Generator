//____________________________________________________________________________
/*!

\class   genie::NtpMCSummary

\brief   MINOS-style Ntuple Class to hold an MC event summary information.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created October 1, 2004

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

  int    probe;     ///< Probe (v,e...) PDG Code
  int    fsl;       ///< Final state primary lepton
  int    tgt;       ///< Nuclear target PDG Code
  int    nucl;      ///< Struck nucleon PDG Code
  int    iqrk;      ///< Struck quark PDG Code
  int    fqrk;      ///< Outgoing quark PDG Code
  int    Z;         ///< Nuclear target: Z
  int    A;         ///< Nuclear target: A
  int    N;         ///< Nuclear target: N    
  int    scat;      ///< Scattering type (QEL, DIS, RES, COH, ...)
  int    proc;      ///< Process type (E/M, CC, NC, ...)
  float  x;         ///< Bjorken x
  float  y;         ///< Inelasticity
  float  z;         ///< 
  float  Q2;        ///< Momentum transfer (>0 = -q2)
  float  W;         ///< Hadronic invariant mass  
  float  xsec;      ///< Cross section for this interaction
  float  v[4];      ///< Interaction vertex (x,y,z,t)
  float  p4p[4];    ///< Probe 4-momentum (px,py,pz,E) in LAB
  float  p4nucl[4]; ///< Struck nucleon 4-momentum (px,py,py,E) in LAB
  float  p4fsl[4];  ///< Final state primary lepton 4-momentum (px,py,py,E) in LAB
  float  p4fsl2[4]; ///< 
  float  p4fsh[4];  ///< Hadronic system 4-momentum (px,py,py,E) in LAB

  ClassDef(NtpMCSummary, 1)    
};

}      // genie namespace

#endif // _NTP_MC_SUMMARY_H_
