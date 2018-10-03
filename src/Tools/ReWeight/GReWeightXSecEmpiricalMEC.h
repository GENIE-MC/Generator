//____________________________________________________________________________
/*!

\class    genie::rew::GReWeightXSecEmpiricalMEC

\brief    Reweighting CCQE GENIE neutrino cross sections

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

          Jim Dobson <J.Dobson07 \at imperial.ac.uk>
          Imperial College London

\created  Aug 1, 2009

\cpright  Copyright (c) 2003-2017, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _G_REWEIGHT_EmpMEC_H_
#define _G_REWEIGHT_EmpMEC_H_

//#define _G_REWEIGHT_EmpMEC_DEBUG_

#include <map>
#include <string>

#include "Tools/ReWeight/GReWeightI.h"

using std::map;
using std::string;

class TFile;
class TNtupleD;

namespace genie {

class XSecAlgorithmI;
class Registry;

namespace rew {

class GReWeightXSecEmpiricalMEC : public GReWeightI {
public:
  GReWeightXSecEmpiricalMEC();
  ~GReWeightXSecEmpiricalMEC();

  // implement the GReWeightI interface
  bool IsHandled(GSyst_t syst);
  void SetSystematic(GSyst_t syst, double val);
  void Reset(void);
  void Reconfigure(void);
  double CalcWeight(const EventRecord &event);

private:
  void Init(void);

  XSecAlgorithmI *fXSecModelDef; ///< default model
  XSecAlgorithmI *fXSecModel;    ///< tweaked model
  Registry *fXSecModelConfig;    ///< config in tweaked model

  double fMq2d_TwkDial;
  double fMass_TwkDial;
  double fWidth_TwkDial;
  double fFracPN_NC_TwkDial;
  double fFracPN_CC_TwkDial;
  double fFracCCQE_TwkDial;
  double fFracNCQE_TwkDial;
  double fFracPN_EM_TwkDial;
  double fFracEMQE_TwkDial;

  double fMq2d_Def;
  double fMass_Def;
  double fWidth_Def;
  double fFracPN_NC_Def;
  double fFracPN_CC_Def;
  double fFracCCQE_Def;
  double fFracNCQE_Def;
  double fFracPN_EM_Def;
  double fFracEMQE_Def;

  double fMq2d_Curr;
  double fMass_Curr;
  double fWidth_Curr;
  double fFracPN_NC_Curr;
  double fFracPN_CC_Curr;
  double fFracCCQE_Curr;
  double fFracNCQE_Curr;
  double fFracPN_EM_Curr;
  double fFracEMQE_Curr;

  bool fAnyTwk;

#ifdef _G_REWEIGHT_EmpMEC_DEBUG_
  TFile *fTestFile;
  TNtupleD *fTestNtp;
#endif
};

} // namespace rew
} // namespace genie

#endif

