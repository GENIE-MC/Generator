//____________________________________________________________________________
/*
 Copyright (c) 2003-2023, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
 University of Liverpool & STFC Rutherford Appleton Laboratory

 Changes required to implement the Electron Velocity module
 were installed by Brinden Carlson (Univ. of Florida)
*/
//____________________________________________________________________________

#include <sstream>

#include <TParticlePDG.h>
#include <TRootIOCtor.h>

#include "Framework/Conventions/Constants.h"
#include "Framework/Interaction/Target.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Utils/PrintUtils.h"

using std::endl;
using std::ostringstream;

using namespace genie;
using namespace genie::constants;

ClassImp(Target)

//____________________________________________________________________________
namespace genie {
 ostream & operator << (ostream & stream, const Target & target)
 {
   target.Print(stream);
   return stream;
 }
}
//___________________________________________________________________________
Target::Target() :
TObject()
{
  this->Init();
}
//___________________________________________________________________________
Target::Target(int pdgc) :
TObject()
{
  this->Init();
  this->SetId(pdgc);
}
//___________________________________________________________________________
Target::Target(int ZZ, int AA) :
TObject()
{
  this->Init();
  this->SetId(ZZ,AA);
}
//___________________________________________________________________________
Target::Target(int ZZ, int AA, int hit_part_pdgc) :
TObject()
{
  this->Init();
  this->SetId(ZZ,AA);
  this->SetHitPartPdg(hit_part_pdgc);
}
//___________________________________________________________________________
Target::Target(const Target & tgt) :
TObject()
{
  this->Init();
  this->Copy(tgt);
}
//___________________________________________________________________________
Target::Target(TRootIOCtor*) :
TObject(),
fZ(0),
fA(0),
fTgtPDG(0),
fHitPartPDG(0),
fHitSeaQrk(false),
<<<<<<< HEAD
fHitNucP4(0),
fHitEleP4(0)
=======
fHitPartP4(nullptr)
>>>>>>> c3e0f096f4c4d645988286b20dd77cde5d717adb
{

}
//___________________________________________________________________________
Target::~Target()
{
  this->CleanUp();
}
//___________________________________________________________________________
void Target::Reset(void)
{
  this->CleanUp();
  this->Init();
}
//___________________________________________________________________________
void Target::Init(void)
{
<<<<<<< HEAD
  fZ         = 0;
  fA         = 0;
  fTgtPDG    = 0;
  fHitNucPDG = 0;
  fHitQrkPDG = 0;
  fHitSeaQrk = false;
  fHitNucP4  = new TLorentzVector(0,0,0,kNucleonMass);
  fHitEleP4  = new TLorentzVector(0,0,0,kElectronMass);
  fHitNucRad = 0.;
=======
  fZ          = 0;
  fA          = 0;
  fTgtPDG     = 0;
  fHitPartPDG = 0;
  fHitQrkPDG  = 0;
  fHitSeaQrk  = false;
  fHitPartP4  = new TLorentzVector(0,0,0,kNucleonMass);
  fHitPartRad = 0.;
>>>>>>> c3e0f096f4c4d645988286b20dd77cde5d717adb
}
//___________________________________________________________________________
void Target::CleanUp(void)
{
<<<<<<< HEAD
  delete fHitNucP4;
  delete fHitEleP4;
=======
  delete fHitPartP4;
>>>>>>> c3e0f096f4c4d645988286b20dd77cde5d717adb
}
//___________________________________________________________________________
void Target::Copy(const Target & tgt)
{
  fTgtPDG = tgt.fTgtPDG;

  if( pdg::IsIon(fTgtPDG) ) {

     fZ          = tgt.fZ; // copy A,Z
     fA          = tgt.fA;
     fHitPartPDG = tgt.fHitPartPDG; // struck nucleon PDG
     fHitQrkPDG  = tgt.fHitQrkPDG; // struck quark PDG
     fHitSeaQrk  = tgt.fHitSeaQrk; // struck quark is from sea?

     //// valgrind warns about this ... try something else
     // (*fHitPartP4) = (*tgt.fHitPartP4);
     const TLorentzVector& p4 = *(tgt.fHitPartP4);
     //  *fHitPartP4 = p4; // nope
     //// this works for valgrind
     fHitPartP4->SetX(p4.X());
     fHitPartP4->SetY(p4.Y());
     fHitPartP4->SetZ(p4.Z());
     fHitPartP4->SetT(p4.T());

     fHitPartRad = tgt.fHitPartRad;

     // look-up the nucleus in the isotopes chart
     this->ForceNucleusValidity();

     // make sure the hit nucleus constituent object is a valid 
     // particle usable in the simulation
     this->ForceHitPartValidity();
  }
  if (tgt.fHitNucPDG == 0 && tgt.fHitQrkPDG == 0 && tgt.fHitSeaQrk == 0){ //No interaction with nucleus -> interaction with electron
    const TLorentzVector& p4 = *(tgt.fHitEleP4);
    *fHitEleP4 = *tgt.fHitEleP4 ;
  }
}
//___________________________________________________________________________
void Target::SetId(int pdgc)
{
  fTgtPDG = pdgc;
  if( pdg::IsIon(pdgc) ) {
     fZ = pdg::IonPdgCodeToZ(pdgc);
     fA = pdg::IonPdgCodeToA(pdgc);
  }

  this->ForceNucleusValidity(); // search at the isotopes chart
  //this->AutoSetHitPart();      // struck nuc := tgt for free nucleon tgt
}
//___________________________________________________________________________
void Target::SetId(int ZZ, int AA)
{
  fTgtPDG = pdg::IonPdgCode(AA,ZZ);
  fZ = ZZ;
  fA = AA;

  this->ForceNucleusValidity(); // search at the isotopes chart
  //this->AutoSetHitPart();      // struck nuc := tgt for free nucleon tgt
}
//___________________________________________________________________________
void Target::SetHitPartPdg(int pdgc)
{
  fHitPartPDG = pdgc;
  bool is_valid = this->ForceHitPartValidity();  // valid particle

  // If it is a valid struck nucleon pdg code, initialize its 4P:
  // at-rest + on-mass-shell
  if(is_valid) {
    double M = PDGLibrary::Instance()->Find(pdgc)->Mass();
    fHitPartP4->SetPxPyPzE(0,0,0,M);
  }
}
//___________________________________________________________________________
void Target::SetHitQrkPdg(int pdgc)
{
  if(pdg::IsQuark(pdgc) || pdg::IsAntiQuark(pdgc)) fHitQrkPDG = pdgc;
}
//___________________________________________________________________________
void Target::SetHitPartP4(const TLorentzVector & p4)
{
  if(fHitPartP4) delete fHitPartP4;
  fHitPartP4 = new TLorentzVector(p4);
}
//___________________________________________________________________________
void Target::SetHitEleP4(const TLorentzVector & p4)
{
  if(fHitEleP4) delete fHitEleP4;
  fHitEleP4 = new TLorentzVector(p4);
}
//___________________________________________________________________________
void Target::SetHitSeaQrk(bool tf)
{
  fHitSeaQrk = tf;
}
//___________________________________________________________________________
void Target::ForceHitPartOnMassShell(void)
{
  if(this->HitPartIsSet()) {
     double m = this->HitPartMass();
     double p = this->HitPartP4Ptr()->P();
     double e = TMath::Sqrt(p*p+m*m);
     this->HitPartP4Ptr()->SetE(e);
  }
}
//___________________________________________________________________________
void Target::SetHitPartPosition(double r)
{
  fHitPartRad = r;
}
//___________________________________________________________________________
double Target::Charge(void) const
{
// Shortcut for commonly used code for extracting the nucleus charge from PDG
//
  TParticlePDG * p = PDGLibrary::Instance()->Find(fTgtPDG);
  if(p) return p->Charge() / 3.; // in +e
  return 0;
}
//___________________________________________________________________________
double Target::Mass(void) const
{
// Shortcut for commonly used code for extracting the nucleus mass from PDG
//
  TParticlePDG * p = PDGLibrary::Instance()->Find(fTgtPDG);
  if(p) return p->Mass(); // in GeV
  return 0.;
}
//___________________________________________________________________________
double Target::HitPartMass(void) const
{
  if(!fHitPartPDG) {
    LOG("Target", pWARN) << "Returning struck particle mass = 0";
    return 0;
  }
  return PDGLibrary::Instance()->Find(fHitPartPDG)->Mass();
}
//___________________________________________________________________________
int Target::HitQrkPdg(void) const
{
  return fHitQrkPDG;
}
//___________________________________________________________________________
TLorentzVector * Target::HitPartP4Ptr(void) const
{
  if(!fHitPartP4) {
    LOG("Target", pWARN) << "Returning NULL struck particle 4-momentum";
    return 0;
  }

  return fHitPartP4;
}
//___________________________________________________________________________
TLorentzVector * Target::HitEleP4Ptr(void) const
{
  if(!fHitEleP4) {
    LOG("Target", pWARN) << "Returning NULL struck electron 4-momentum";
    return 0;
  }
  return fHitEleP4;
}
//___________________________________________________________________________
bool Target::IsFreeNucleon(void) const
{
  return (fA == 1 && (fZ == 0 || fZ == 1));
}
//___________________________________________________________________________
bool Target::IsProton(void) const
{
  return (fA == 1 && fZ == 1);
}
//___________________________________________________________________________
bool Target::IsElectron(void) const
{
  return (fA == 0 && fZ == 0); //No nucleons
}
//___________________________________________________________________________
bool Target::IsNeutron(void) const
{
  return (fA == 1 && fZ == 0);
}
//___________________________________________________________________________
bool Target::IsNucleus(void) const
{
  return (fA > 1); // IsValidNucleus() was ensured when A,Z were set
}
//___________________________________________________________________________
bool Target::IsParticle(void) const
{
  TParticlePDG * p = PDGLibrary::Instance()->Find(fTgtPDG);
  return (p && fA==0 && fZ==0);
}
//___________________________________________________________________________
bool Target::HitPartIsSet(void) const
{
  bool ok =
     pdg::IsNucleon(fHitPartPDG)          ||
     pdg::Is2NucleonCluster(fHitPartPDG)  || 
     pdg::IsElectron(fHitPartPDG);

  return ok;
}
//___________________________________________________________________________
bool Target::HitEleIsSet(void) const
{
  bool ok = fHitNucPDG == 0 &&
            fHitQrkPDG == 0 &&
            fHitSeaQrk == 0; //Hit no quarks

  return ok;
}
//___________________________________________________________________________
bool Target::HitQrkIsSet(void) const
{
  return (
     pdg::IsQuark(fHitQrkPDG) || pdg::IsAntiQuark(fHitQrkPDG)
  );
}
//___________________________________________________________________________
bool Target::HitSeaQrk(void) const
{
  return fHitSeaQrk;
}
//___________________________________________________________________________
int Target::HitPartPdg(void) const
{
  return fHitPartPDG;
}
//___________________________________________________________________________
bool Target::IsValidNucleus(void) const
{
  //-- it is valid if it is a free nucleon...
  if(this->IsFreeNucleon()) return true;

  //-- ... or a nucleus that can be found in the MINOS ion PDG extensions
  int pdg_code = pdg::IonPdgCode(fA, fZ);
  TParticlePDG * p = PDGLibrary::Instance()->Find(pdg_code);
  if(p) return true;

  return false;
}
//___________________________________________________________________________
bool Target::IsEvenEven(void) const
{
  if( this->IsNucleus() ) {
    int NN = this->N();
    int ZZ = this->Z();
    if( NN % 2 == 0 && ZZ % 2 == 0 ) return true;
  }
  return false;
}
//___________________________________________________________________________
bool Target::IsEvenOdd(void) const
{
  if( this->IsNucleus() ) {
      if(! this->IsEvenEven() && ! this->IsOddOdd() ) return true;
  }
  return false;
}
//___________________________________________________________________________
bool Target::IsOddOdd(void) const
{
  if( this->IsNucleus() ) {
    int NN = this->N();
    int ZZ = this->Z();
    if( NN % 2 == 1 && ZZ % 2 == 1 ) return true;
  }
  return false;
}
//___________________________________________________________________________
bool Target::ForceHitPartValidity(void)
{
// resets the struck part pdg-code if it is found not to be a valid one

  bool valid =
      pdg::IsNucleon(fHitPartPDG)          ||
      pdg::Is2NucleonCluster (fHitPartPDG) ||
      pdg::IsElectron(fHitPartPDG)         ||
      (fHitPartPDG==0); /* not set */

  return valid;
}
//___________________________________________________________________________
void Target::ForceNucleusValidity(void)
{
// resets the target pdg-code if it is found not to be a valid one

  if( ! this->IsValidNucleus() ) {
    LOG("Target", pWARN) << "Invalid target -- Reseting to Z = 0, A = 0";
    fZ = 0;
    fA = 0;
  }
}
//___________________________________________________________________________
void Target::AutoSetHitPart(void)
{
// for free nucleon targets -> (auto)set struck nucleon = target

  if( this->IsFreeNucleon() ) {
    if( this->IsProton() ) this->SetHitPartPdg(kPdgProton);
    else                   this->SetHitPartPdg(kPdgNeutron);
  }
}
//___________________________________________________________________________
string Target::AsString(void) const
{
  ostringstream s;

  s << this->Pdg();
  if(this->HitPartIsSet())
     s << "[Part=" << this->HitPartPdg() << "]";
  if(this->HitQrkIsSet()) {
     s << "[q=" << this->HitQrkPdg();
     s << (this->HitSeaQrk() ? "(s)" : "(v)");
     s << "]";
  }

  return s.str();
}
//___________________________________________________________________________
void Target::Print(ostream & stream) const
{
  stream << " target PDG code = " << fTgtPDG << endl;

  if( this->IsNucleus() || this->IsFreeNucleon() ) {
      stream << " Z = " << fZ << ", A = " << fA << endl;
  }

  if( this->HitPartIsSet() ) {
    TParticlePDG * p = PDGLibrary::Instance()->Find(fHitPartPDG);
    stream << " struck Part = " << p->GetName()
           << ", P4 = " << utils::print::P4AsString(fHitPartP4) << endl;
  }

  if( this->HitQrkIsSet() ) {
    TParticlePDG * q = PDGLibrary::Instance()->Find(fHitQrkPDG);
    stream << " struck quark = " << q->GetName()
           << " (from sea: "
           << utils::print::BoolAsYNString(this->HitSeaQrk())
           << ")";
  }
  if( this->HitEleIsSet() ) {
    TParticlePDG * q = PDGLibrary::Instance()->Find(fHitQrkPDG);
    stream << " struck electron = ";
  }
}
//___________________________________________________________________________
bool Target::Compare(const Target & target) const
{
  int  tgt_pdg         = target.Pdg();
  int  struck_part_pdg = target.HitPartPdg();
  int  struck_qrk_pdg  = target.HitQrkPdg();
  bool struck_sea_qrk  = target.HitSeaQrk();

  bool equal = ( fTgtPDG     == tgt_pdg         ) &&
               ( fHitPartPDG == struck_part_pdg ) &&
               ( fHitQrkPDG  == struck_qrk_pdg  ) &&
               ( fHitSeaQrk  == struck_sea_qrk  );
  return equal;
}
//___________________________________________________________________________
bool Target::operator == (const Target & target) const
{
  return this->Compare(target);
}
//___________________________________________________________________________
Target & Target::operator = (const Target & target)
{
  this->Copy(target);
  return (*this);
}
//___________________________________________________________________________
