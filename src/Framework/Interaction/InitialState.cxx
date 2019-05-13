//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab 

         Changes required to implement the GENIE Boosted Dark Matter module 
         were installed by Josh Berger (Univ. of Wisconsin)

         CMEnergy() method added by Andy Furmanski (Univ. of Manchester) 
         and Joe Johnston (Univ of Pittsburgh)
*/
//____________________________________________________________________________

#include <cassert>
#include <sstream>
#include <iomanip>

#include <TRootIOCtor.h>

#include "Framework/Interaction/InitialState.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGUtils.h"

using namespace genie;

using std::endl;
using std::setprecision;
using std::setw;
using std::ostringstream;

ClassImp(InitialState)

//____________________________________________________________________________
namespace genie {
 ostream & operator << (ostream & stream, const InitialState & init_state)
 {
   init_state.Print(stream);
   return stream;
 }
}
//___________________________________________________________________________
InitialState::InitialState() :
TObject()
{
  this->Init();
}
//___________________________________________________________________________
InitialState::InitialState(int target_pdgc, int probe_pdgc) :
TObject()
{
  this->Init(target_pdgc, probe_pdgc);
}
//___________________________________________________________________________
InitialState::InitialState(int Z, int A, int probe_pdgc) :
TObject()
{
  int target_pdgc = pdg::IonPdgCode(A,Z);
  this->Init(target_pdgc, probe_pdgc);
}
//___________________________________________________________________________
InitialState::InitialState(const Target & tgt, int probe_pdgc) :
TObject()
{
  int target_pdgc = tgt.Pdg();
  this->Init(target_pdgc, probe_pdgc);
}
//___________________________________________________________________________
InitialState::InitialState(const InitialState & init_state) :
TObject()
{
  this->Init();
  this->Copy(init_state);
}
//___________________________________________________________________________
InitialState::InitialState(TRootIOCtor*) :
TObject(),
fProbePdg(0),
fTgt(0), 
fProbeP4(0), 
fTgtP4(0)
{

}
//___________________________________________________________________________
InitialState::~InitialState()
{
  this->CleanUp();
}
//___________________________________________________________________________
void InitialState::Init(void)
{
  fProbePdg  = 0;
  fTgt       = new Target();
  fProbeP4   = new TLorentzVector(0, 0, 0, 0);
  fTgtP4     = new TLorentzVector(0, 0, 0, 0);
}
//___________________________________________________________________________
void InitialState::Init(int target_pdgc, int probe_pdgc)
{
  TParticlePDG * t = PDGLibrary::Instance()->Find(target_pdgc);
  TParticlePDG * p = PDGLibrary::Instance()->Find(probe_pdgc );

  assert(t && p);

  double m = p->Mass();
  double M = t->Mass();

  fProbePdg  = probe_pdgc;
  fTgt       = new Target(target_pdgc);
  fProbeP4   = new TLorentzVector(0, 0, 0, m);
  fTgtP4     = new TLorentzVector(0, 0, 0, M);
}
//___________________________________________________________________________
void InitialState::CleanUp(void)
{
  delete fTgt;
  delete fProbeP4;
  delete fTgtP4;
}
//___________________________________________________________________________
void InitialState::Reset(void)
{
  this->CleanUp();
  this->Init();
}
//___________________________________________________________________________
void InitialState::Copy(const InitialState & init_state)
{
  fProbePdg = init_state.fProbePdg;

  fTgt->Copy(*init_state.fTgt);

  this -> SetProbeP4 ( *init_state.fProbeP4 );
  this -> SetTgtP4   ( *init_state.fTgtP4   );
}
//___________________________________________________________________________
int InitialState::TgtPdg(void) const
{
  assert(fTgt);
  return fTgt->Pdg();
}
//___________________________________________________________________________
TParticlePDG * InitialState::Probe(void) const
{
  TParticlePDG * p = PDGLibrary::Instance()->Find(fProbePdg);
  return p;
}
//___________________________________________________________________________
void InitialState::SetPdgs(int tgt_pdgc, int probe_pdgc)
{
  this->CleanUp();
  this->Init(tgt_pdgc, probe_pdgc);
}
//___________________________________________________________________________
void InitialState::SetTgtPdg(int tgt_pdgc)
{
  int probe_pdgc = this->ProbePdg();

  this->CleanUp();
  this->Init(tgt_pdgc, probe_pdgc);
}
//___________________________________________________________________________
void InitialState::SetProbePdg(int probe_pdgc)
{
  int tgt_pdgc = this->TgtPdg();

  this->CleanUp();
  this->Init(tgt_pdgc, probe_pdgc);
}
//___________________________________________________________________________
void InitialState::SetProbeE(double E)
{
  fProbeP4 -> SetE  ( E );
  fProbeP4 -> SetPx ( 0.);
  fProbeP4 -> SetPy ( 0.);
  fProbeP4 -> SetPz ( E );
}
//___________________________________________________________________________
void InitialState::SetProbeP4(const TLorentzVector & P4)
{
  fProbeP4 -> SetE  ( P4.E()  );
  fProbeP4 -> SetPx ( P4.Px() );
  fProbeP4 -> SetPy ( P4.Py() );
  fProbeP4 -> SetPz ( P4.Pz() );
}
//___________________________________________________________________________
void InitialState::SetTgtP4(const TLorentzVector & P4)
{
  fTgtP4 -> SetE  ( P4.E()  );
  fTgtP4 -> SetPx ( P4.Px() );
  fTgtP4 -> SetPy ( P4.Py() );
  fTgtP4 -> SetPz ( P4.Pz() );
}
//___________________________________________________________________________
bool InitialState::IsNuP(void) const
{
  int  prob = fProbePdg;
  int  nucl = fTgt->HitNucPdg();
  bool isvp = pdg::IsNeutrino(prob) && pdg::IsProton(nucl);

  return isvp;
}
//___________________________________________________________________________
bool InitialState::IsNuN(void) const
{
  int  prob = fProbePdg;
  int  nucl = fTgt->HitNucPdg();
  bool isvn = pdg::IsNeutrino(prob) && pdg::IsNeutron(nucl);

  return isvn;
}
//___________________________________________________________________________
bool InitialState::IsNuBarP(void) const
{
  int  prob  = fProbePdg;
  int  nucl  = fTgt->HitNucPdg();
  bool isvbp = pdg::IsAntiNeutrino(prob) && pdg::IsProton(nucl);

  return isvbp;
}
//___________________________________________________________________________
bool InitialState::IsNuBarN(void) const
{
  int  prob  = fProbePdg;
  int  nucl  = fTgt->HitNucPdg();
  bool isvbn = pdg::IsAntiNeutrino(prob) && pdg::IsNeutron(nucl);

  return isvbn;
}
//___________________________________________________________________________
bool InitialState::IsDMP(void) const
{
// Check if DM - proton interaction
  int  prob = fProbePdg;
  int  nucl = fTgt->HitNucPdg();
  bool isdp = pdg::IsDarkMatter(prob) && pdg::IsProton(nucl);

  return isdp;
}
//___________________________________________________________________________
bool InitialState::IsDMN(void) const
{
// Check if DM - neutron interaction
  int  prob = fProbePdg;
  int  nucl = fTgt->HitNucPdg();
  bool isdn = pdg::IsDarkMatter(prob) && pdg::IsNeutron(nucl);

  return isdn;
}
//___________________________________________________________________________
TLorentzVector * InitialState::GetTgtP4(RefFrame_t ref_frame) const
{
// Return the target 4-momentum in the specified reference frame
// Note: the caller adopts the TLorentzVector object

  switch (ref_frame) {

       //------------------ NUCLEAR TARGET REST FRAME:
       case (kRfTgtRest) :
       {
             // for now make sure that the target nucleus is always at
             // rest and it is only the struck nucleons that can move:
             // so the [target rest frame] = [LAB frame]

             return this->GetTgtP4(kRfLab);
       }

       //------------------ STRUCK NUCLEON REST FRAME:
       case (kRfHitNucRest) :
       {
             // make sure that 'struck nucleon' properties were set in
             // the nuclear target object
             assert(fTgt->HitNucIsSet());
             TLorentzVector * pnuc4 = fTgt->HitNucP4Ptr();

             // compute velocity vector (px/E, py/E, pz/E)
             double bx = pnuc4->Px() / pnuc4->Energy();
             double by = pnuc4->Py() / pnuc4->Energy();
             double bz = pnuc4->Pz() / pnuc4->Energy();

             // BOOST
             TLorentzVector * p4 = new TLorentzVector(*fTgtP4);
             p4->Boost(-bx,-by,-bz);

             return p4;
             break;
       }
       //------------------ LAB:
       case (kRfLab) :
       {
             TLorentzVector * p4 = new TLorentzVector(*fTgtP4);
             return p4;
             break;
       }
       default:
             LOG("Interaction", pERROR) << "Uknown reference frame";
  }
  return 0;
}
//___________________________________________________________________________
TLorentzVector * InitialState::GetProbeP4(RefFrame_t ref_frame) const
{
// Return the probe 4-momentum in the specified reference frame
// Note: the caller adopts the TLorentzVector object

  switch (ref_frame) {

       //------------------ NUCLEAR TARGET REST FRAME:
       case (kRfTgtRest) :
       {
             // for now assure that the target nucleus is always at rest
             // and it is only the struck nucleons that can move:
             // so the [target rest frame] = [LAB frame]

             TLorentzVector * p4 = new TLorentzVector(*fProbeP4);
             return p4;
       }

       //------------------ STRUCK NUCLEON REST FRAME:
       case (kRfHitNucRest) :
       {
             // make sure that 'struck nucleon' properties were set in
             // the nuclear target object

             assert( fTgt->HitNucP4Ptr() != 0 );

             TLorentzVector * pnuc4 = fTgt->HitNucP4Ptr();

             // compute velocity vector (px/E, py/E, pz/E)

             double bx = pnuc4->Px() / pnuc4->Energy();
             double by = pnuc4->Py() / pnuc4->Energy();
             double bz = pnuc4->Pz() / pnuc4->Energy();

             // BOOST

             TLorentzVector * p4 = new TLorentzVector(*fProbeP4);

             p4->Boost(-bx,-by,-bz);

             return p4;

             break;
       }
       //------------------ LAB:
       case (kRfLab) :
       {
             TLorentzVector * p4 = new TLorentzVector(*fProbeP4);
             return p4;

             break;
       }
       default:

             LOG("Interaction", pERROR) << "Uknown reference frame";
  }
  return 0;
}
//___________________________________________________________________________
double InitialState::ProbeE(RefFrame_t ref_frame) const
{
  TLorentzVector * p4 = this->GetProbeP4(ref_frame);
  double E = p4->Energy();

  delete p4;
  return E;
}

//___________________________________________________________________________
double InitialState::CMEnergy() const
{
  TLorentzVector * k4 = this->GetProbeP4(kRfLab);
  TLorentzVector * p4 = fTgt->HitNucP4Ptr();
  
  *k4 += *p4; // now k4 represents centre-of-mass 4-momentum
  double s = k4->Dot(*k4); // dot-product with itself
  double E = TMath::Sqrt(s);

  delete k4;
  
  return E;
}
//___________________________________________________________________________
string InitialState::AsString(void) const
{
// Code-ify the interaction in a string to be used as (part of a) keys
// Template:
//     nu_pdg:code;tgt-pdg:code;

  ostringstream init_state;

  if (this->Probe()->Mass() > 0) {
    init_state << "dm_mass:" << this->Probe()->Mass() << ";";
  }
  else {
    init_state << "nu-pdg:"  << this->ProbePdg()  << ";";
  }
  init_state << "tgt-pdg:" << this->Tgt().Pdg() << ";";

  return init_state.str();
}
//___________________________________________________________________________
void InitialState::Print(ostream & stream) const
{
  stream << "[-] [Init-State] " << endl;

  stream << " |--> probe        : "
         << "PDG-code = " << fProbePdg
         << " (" << this->Probe()->GetName() << ")" << endl;

  stream << " |--> nucl. target : "
         << "Z = "          << fTgt->Z()
         << ", A = "        << fTgt->A()
         << ", PDG-Code = " << fTgt->Pdg();

  TParticlePDG * tgt = PDGLibrary::Instance()->Find( fTgt->Pdg() );
  if(tgt) {
    stream << " (" << tgt->GetName() << ")";
  }
  stream << endl;

  stream << " |--> hit nucleon  : ";
  int nuc_pdgc = fTgt->HitNucPdg();

  if ( pdg::IsNeutronOrProton(nuc_pdgc) ) {
    TParticlePDG * p = PDGLibrary::Instance()->Find(nuc_pdgc);
    stream << "PDC-Code = " << nuc_pdgc << " (" << p->GetName() << ")";
  } else {
    stream << "no set";
  }
  stream << endl;

  stream << " |--> hit quark    : ";
  int qrk_pdgc = fTgt->HitQrkPdg();

  if ( pdg::IsQuark(qrk_pdgc) || pdg::IsAntiQuark(qrk_pdgc)) {
    TParticlePDG * p = PDGLibrary::Instance()->Find(qrk_pdgc);
    stream << "PDC-Code = " << qrk_pdgc << " (" << p->GetName() << ") ";
    stream << (fTgt->HitSeaQrk() ? "[sea]" : "[valence]");
  } else {
    stream << "no set";
  }
  stream << endl;

  stream << " |--> probe 4P     : "
         << "(E = "   << setw(12) << setprecision(6) << fProbeP4->E()
         << ", Px = " << setw(12) << setprecision(6) << fProbeP4->Px()
         << ", Py = " << setw(12) << setprecision(6) << fProbeP4->Py()
         << ", Pz = " << setw(12) << setprecision(6) << fProbeP4->Pz()
         << ")"
         << endl;
  stream << " |--> target 4P    : "
         << "(E = "   << setw(12) << setprecision(6) << fTgtP4->E()
         << ", Px = " << setw(12) << setprecision(6) << fTgtP4->Px()
         << ", Py = " << setw(12) << setprecision(6) << fTgtP4->Py()
         << ", Pz = " << setw(12) << setprecision(6) << fTgtP4->Pz()
         << ")"
         << endl;

  if ( pdg::IsNeutronOrProton(nuc_pdgc) ) {

    TLorentzVector * nuc_p4 = fTgt->HitNucP4Ptr();

    stream << " |--> nucleon 4P   : "
           << "(E = "   << setw(12) << setprecision(6) << nuc_p4->E()
           << ", Px = " << setw(12) << setprecision(6) << nuc_p4->Px()
           << ", Py = " << setw(12) << setprecision(6) << nuc_p4->Py()
           << ", Pz = " << setw(12) << setprecision(6) << nuc_p4->Pz()
           << ")";
  }
}
//___________________________________________________________________________
bool InitialState::Compare(const InitialState & init_state) const
{
  int            probe  = init_state.ProbePdg();
  const Target & target = init_state.Tgt();

  bool equal = (fProbePdg == probe) && (*fTgt == target);

  return equal;
}
//___________________________________________________________________________
bool InitialState::operator == (const InitialState & init_state) const
{
  return this->Compare(init_state);
}
//___________________________________________________________________________
InitialState & InitialState::operator = (const InitialState & init_state)
{
  this->Copy(init_state);
  return (*this);
}
//___________________________________________________________________________


