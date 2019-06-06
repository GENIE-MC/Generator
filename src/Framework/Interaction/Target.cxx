//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Oct 09, 2007 - CA
   Hit nucleon not auto-set for hit nucleon targets
 @ May 05, 2010 - CR
   Adding special ctor for ROOT I/O purposes so as to avoid memory leak due to
   memory allocated in the default ctor when objects of this class are read by
   the ROOT Streamer.
 @ Nov 28, 2011 - CA
   Now a nucleon-cluster ID is an accepted option for SetHitNucPdg().
 @ Mar 18, 2016 - JJ (SD)
   Added methods to store and retrieve the struck nucleon position. Position
   is stored as a double, indicating the distance from the center of the
   nucleus.
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
Target::Target(int ZZ, int AA, int hit_nucleon_pdgc) :
TObject()
{
  this->Init();
  this->SetId(ZZ,AA);
  this->SetHitNucPdg(hit_nucleon_pdgc);
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
fHitNucPDG(0),
fHitSeaQrk(false),
fHitNucP4(0)
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
  fZ         = 0;
  fA         = 0;
  fTgtPDG    = 0;
  fHitNucPDG = 0;
  fHitQrkPDG = 0;
  fHitSeaQrk = false;
  fHitNucP4  = new TLorentzVector(0,0,0,kNucleonMass);
  fHitNucRad = 0.;
}
//___________________________________________________________________________
void Target::CleanUp(void)
{
  delete fHitNucP4;
}
//___________________________________________________________________________
void Target::Copy(const Target & tgt)
{
  fTgtPDG = tgt.fTgtPDG;

  if( pdg::IsIon(fTgtPDG) ) {

     fZ         = tgt.fZ; // copy A,Z
     fA         = tgt.fA;
     fHitNucPDG = tgt.fHitNucPDG; // struck nucleon PDG
     fHitQrkPDG = tgt.fHitQrkPDG; // struck quark PDG
     fHitSeaQrk = tgt.fHitSeaQrk; // struck quark is from sea?

     //// valgrind warns about this ... try something else
     // (*fHitNucP4) = (*tgt.fHitNucP4);
     const TLorentzVector& p4 = *(tgt.fHitNucP4);
     //  *fHitNucP4 = p4; // nope
     //// this works for valgrind
     fHitNucP4->SetX(p4.X());
     fHitNucP4->SetY(p4.Y());
     fHitNucP4->SetZ(p4.Z());
     fHitNucP4->SetT(p4.T());

     fHitNucRad = tgt.fHitNucRad;

     // look-up the nucleus in the isotopes chart
     this->ForceNucleusValidity();

     // make sure the hit nucleus constituent object is either
     // a nucleon (p or n) or a di-nucleon cluster (p+p, p+n, n+n)
     this->ForceHitNucValidity();
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
  //this->AutoSetHitNuc();      // struck nuc := tgt for free nucleon tgt
}
//___________________________________________________________________________
void Target::SetId(int ZZ, int AA)
{
  fTgtPDG = pdg::IonPdgCode(AA,ZZ);
  fZ = ZZ;
  fA = AA;

  this->ForceNucleusValidity(); // search at the isotopes chart
  //this->AutoSetHitNuc();      // struck nuc := tgt for free nucleon tgt
}
//___________________________________________________________________________
void Target::SetHitNucPdg(int nucl_pdgc)
{
  fHitNucPDG = nucl_pdgc;
  bool is_valid = this->ForceHitNucValidity();  // p, n or a di-nucleon

  // If it is a valid struck nucleon pdg code, initialize its 4P:
  // at-rest + on-mass-shell
  if(is_valid) {
    double M = PDGLibrary::Instance()->Find(nucl_pdgc)->Mass();
    fHitNucP4->SetPxPyPzE(0,0,0,M);
  }
}
//___________________________________________________________________________
void Target::SetHitQrkPdg(int pdgc)
{
  if(pdg::IsQuark(pdgc) || pdg::IsAntiQuark(pdgc)) fHitQrkPDG = pdgc;
}
//___________________________________________________________________________
void Target::SetHitNucP4(const TLorentzVector & p4)
{
  if(fHitNucP4) delete fHitNucP4;
  fHitNucP4 = new TLorentzVector(p4);
}
//___________________________________________________________________________
void Target::SetHitSeaQrk(bool tf)
{
  fHitSeaQrk = tf;
}
//___________________________________________________________________________
void Target::ForceHitNucOnMassShell(void)
{
  if(this->HitNucIsSet()) {
     double m = this->HitNucMass();
     double p = this->HitNucP4Ptr()->P();
     double e = TMath::Sqrt(p*p+m*m);
     this->HitNucP4Ptr()->SetE(e);
  }
}
//___________________________________________________________________________
void Target::SetHitNucPosition(double r)
{
  fHitNucRad = r;
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
double Target::HitNucMass(void) const
{
  if(!fHitNucPDG) {
    LOG("Target", pWARN) << "Returning struck nucleon mass = 0";
    return 0;
  }
  return PDGLibrary::Instance()->Find(fHitNucPDG)->Mass();
}
//___________________________________________________________________________
int Target::HitQrkPdg(void) const
{
  return fHitQrkPDG;
}
//___________________________________________________________________________
TLorentzVector * Target::HitNucP4Ptr(void) const
{
  if(!fHitNucP4) {
    LOG("Target", pWARN) << "Returning NULL struck nucleon 4-momentum";
    return 0;
  }

  return fHitNucP4;
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
bool Target::HitNucIsSet(void) const
{
  bool ok =
     pdg::IsNucleon(fHitNucPDG)          ||
     pdg::Is2NucleonCluster (fHitNucPDG);

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
int Target::HitNucPdg(void) const
{
  return fHitNucPDG;
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
bool Target::ForceHitNucValidity(void)
{
// resets the struck nucleon pdg-code if it is found not to be a valid one

  bool valid =
      pdg::IsNucleon(fHitNucPDG)          ||
      pdg::Is2NucleonCluster (fHitNucPDG) ||
      (fHitNucPDG==0); /* not set */

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
void Target::AutoSetHitNuc(void)
{
// for free nucleon targets -> (auto)set struck nucleon = target

  if( this->IsFreeNucleon() ) {
    if( this->IsProton() ) this->SetHitNucPdg(kPdgProton);
    else                   this->SetHitNucPdg(kPdgNeutron);
  }
}
//___________________________________________________________________________
string Target::AsString(void) const
{
  ostringstream s;

  s << this->Pdg();
  if(this->HitNucIsSet())
     s << "[N=" << this->HitNucPdg() << "]";
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

  if( this->HitNucIsSet() ) {
    TParticlePDG * p = PDGLibrary::Instance()->Find(fHitNucPDG);
    stream << " struck nucleon = " << p->GetName()
              << ", P4 = " << utils::print::P4AsString(fHitNucP4) << endl;
  }

  if( this->HitQrkIsSet() ) {
    TParticlePDG * q = PDGLibrary::Instance()->Find(fHitQrkPDG);
    stream << " struck quark = " << q->GetName()
           << " (from sea: "
           << utils::print::BoolAsYNString(this->HitSeaQrk())
           << ")";
  }
}
//___________________________________________________________________________
bool Target::Compare(const Target & target) const
{
  int  tgt_pdg        = target.Pdg();
  int  struck_nuc_pdg = target.HitNucPdg();
  int  struck_qrk_pdg = target.HitQrkPdg();
  bool struck_sea_qrk = target.HitSeaQrk();

  bool equal = ( fTgtPDG    == tgt_pdg        ) &&
               ( fHitNucPDG == struck_nuc_pdg ) &&
               ( fHitQrkPDG == struck_qrk_pdg ) &&
               ( fHitSeaQrk == struck_sea_qrk );
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
