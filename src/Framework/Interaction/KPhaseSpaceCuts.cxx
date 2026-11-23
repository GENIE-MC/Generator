//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#include <cmath>
#include <cstdlib>

#include <TMath.h>

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Interaction/KPhaseSpaceCuts.h"
#include "Framework/Interaction/ProcessInfo.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Registry/Registry.h"

using namespace genie;

//____________________________________________________________________________
KPhaseSpaceCuts * KPhaseSpaceCuts::fInstance = 0;
//____________________________________________________________________________
KPhaseSpaceCuts::KPhaseSpaceCuts() :
fLoaded(false),
fHasEMQ2Min(false),
fHasWeakQ2Min(false),
fEMQ2Min(-1.),
fWeakQ2Min(-1.)
{
  fInstance = 0;
}
//____________________________________________________________________________
KPhaseSpaceCuts::~KPhaseSpaceCuts()
{
  fInstance = 0;
}
//____________________________________________________________________________
KPhaseSpaceCuts * KPhaseSpaceCuts::Instance(void)
{
  if(fInstance == 0) {
    static KPhaseSpaceCuts::Cleaner cleaner;
    cleaner.DummyMethodAndSilentCompiler();
    fInstance = new KPhaseSpaceCuts;
  }
  return fInstance;
}
//____________________________________________________________________________
bool KPhaseSpaceCuts::HasQ2MinCut(const Interaction * interaction) const
{
  this->LoadConfig();

  if(!interaction) return fHasEMQ2Min || fHasWeakQ2Min;

  const ProcessInfo & proc = interaction->ProcInfo();
  if(proc.IsEM()) {
    if(!fHasEMQ2Min) this->FailMissingEMQ2Min();
    return true;
  }

  if(proc.IsWeak()) {
    return fHasWeakQ2Min;
  }

  return false;
}
//____________________________________________________________________________
double KPhaseSpaceCuts::Q2MinCut(
  const Interaction * interaction, double default_Q2min) const
{
  this->LoadConfig();

  if(!interaction) return default_Q2min;

  const ProcessInfo & proc = interaction->ProcInfo();
  if(proc.IsEM()) {
    if(!fHasEMQ2Min) this->FailMissingEMQ2Min();
    return this->RaiseDefault(default_Q2min, fEMQ2Min);
  }

  if(proc.IsWeak()) {
    if(fHasWeakQ2Min) return this->RaiseDefault(default_Q2min, fWeakQ2Min);
  }

  return default_Q2min;
}
//____________________________________________________________________________
bool KPhaseSpaceCuts::HasSplineQ2MinCut(void) const
{
  this->LoadConfig();

  return fHasWeakQ2Min || fHasEMQ2Min;
}
//____________________________________________________________________________
double KPhaseSpaceCuts::SplineQ2MinCut(void) const
{
  this->LoadConfig();

  if(fHasWeakQ2Min)     return fWeakQ2Min;
  if(fHasEMQ2Min)       return fEMQ2Min;
  return -1.;
}
//____________________________________________________________________________
string KPhaseSpaceCuts::SplineQ2MinCutSource(void) const
{
  this->LoadConfig();

  if(fHasWeakQ2Min || fHasEMQ2Min) return "CommonPhaseSpaceCuts.xml";
  return "";
}
//____________________________________________________________________________
bool KPhaseSpaceCuts::HasSplineQ2MinCutForProbe(int probe_pdg) const
{
  this->LoadConfig();

  if(pdg::IsNeutrino(probe_pdg) || pdg::IsAntiNeutrino(probe_pdg)) {
    return fHasWeakQ2Min;
  }
  if(pdg::IsChargedLepton(probe_pdg)) {
    return fHasEMQ2Min;
  }
  return fHasWeakQ2Min || fHasEMQ2Min;
}
//____________________________________________________________________________
double KPhaseSpaceCuts::SplineQ2MinCutForProbe(int probe_pdg) const
{
  this->LoadConfig();

  if(pdg::IsNeutrino(probe_pdg) || pdg::IsAntiNeutrino(probe_pdg)) {
    if(fHasWeakQ2Min) return fWeakQ2Min;
    return -1.;
  }
  if(pdg::IsChargedLepton(probe_pdg)) {
    if(fHasEMQ2Min) return fEMQ2Min;
    return -1.;
  }
  if(fHasWeakQ2Min) return fWeakQ2Min;
  if(fHasEMQ2Min)   return fEMQ2Min;
  return -1.;
}
//____________________________________________________________________________
string KPhaseSpaceCuts::SplineQ2MinCutSourceForProbe(int probe_pdg) const
{
  this->LoadConfig();

  if(!this->HasSplineQ2MinCutForProbe(probe_pdg)) return "";
  return "CommonPhaseSpaceCuts.xml";
}
//____________________________________________________________________________
void KPhaseSpaceCuts::LoadConfig(void) const
{
  if(fLoaded) return;

  fLoaded = true;

  AlgConfigPool * confp = AlgConfigPool::Instance();
  const Registry * r = confp->CommonList("PhaseSpaceCuts", "Default");

  if(!r) {
    LOG("KPhaseSpaceCuts", pWARN)
      << "CommonPhaseSpaceCuts.xml [Default] was not found";
    return;
  }

  if(r->Exists("EM-Q2-min")) {
    fEMQ2Min = this->ValidateQ2Min("EM-Q2-min", r->GetDouble("EM-Q2-min"));
    fHasEMQ2Min = true;
  }

  if(r->Exists("Weak-Q2-min")) {
    fWeakQ2Min = this->ValidateQ2Min("Weak-Q2-min", r->GetDouble("Weak-Q2-min"));
    fHasWeakQ2Min = true;
  }
}
//____________________________________________________________________________
double KPhaseSpaceCuts::ValidateQ2Min(const char * name, double Q2min) const
{
  if(!std::isfinite(Q2min) || Q2min <= 0.) {
    LOG("KPhaseSpaceCuts", pFATAL)
      << "Invalid " << name << " phase-space cut: " << Q2min
      << " GeV^2. Q2 minimum cuts must be finite and > 0.";
    exit(78);
  }
  return Q2min;
}
//____________________________________________________________________________
double KPhaseSpaceCuts::RaiseDefault(double default_Q2min, double cut_Q2min) const
{
  return TMath::Max(default_Q2min, cut_Q2min);
}
//____________________________________________________________________________
void KPhaseSpaceCuts::FailMissingEMQ2Min(void) const
{
  LOG("KPhaseSpaceCuts", pFATAL)
    << "Electromagnetic interactions require EM-Q2-min in "
    << "CommonPhaseSpaceCuts.xml [Default].";
  exit(78);
}
//____________________________________________________________________________
