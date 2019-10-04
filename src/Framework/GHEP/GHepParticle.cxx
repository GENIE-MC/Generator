//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Sep 15, 2009 - CA   
   Remove redundant IsFake(), IsNucleus() and IsParticle() methods. Removed 
   the corresponding priv data members. Use pdg::IsPseudoParticle(), 
   pdg::IsIon() and pdg::IsParticle instead.
 @ Sep 15, 2009 - CA   
   Added 'rescattering code' to allow intranuclear hadron transport MCs to
   store a hadronic reaction code which can not always be easily recreated
   from the particle list.
 @ May 05, 2010 - CR
   Adding special ctor for ROOT I/O purposes so as to avoid memory leak due to
   memory allocated in the default ctor when objects of this class are read by 
   the ROOT Streamer. 

*/
//____________________________________________________________________________

#include <cstdlib>
#include <cassert>
#include <iomanip>

#include <TMath.h>
#include <TRootIOCtor.h>

#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Utils/PrintUtils.h"

using namespace genie;
using namespace genie::constants;

using std::setw;
using std::setprecision;
using std::setfill;
using std::ios;
using std::cout;

const double kPCutOff    = 1e-15;
const double kOffShellDm = 0.002; // 2 MeV

ClassImp(GHepParticle)

//____________________________________________________________________________
namespace genie {
 ostream & operator<< (ostream& stream, const GHepParticle & particle)
 {
   particle.Print(stream);
   return stream;
 }
}
//___________________________________________________________________________
GHepParticle::GHepParticle() :
TObject()
{
  this->Init();
}
//___________________________________________________________________________
// TParticle-like constructor
GHepParticle::GHepParticle(int pdg, GHepStatus_t status,
        int mother1, int mother2, int daughter1, int daughter2,
        const TLorentzVector & p, const TLorentzVector & v) :
TObject(),
fStatus(status),
fFirstMother(mother1),
fLastMother(mother2),
fFirstDaughter(daughter1),
fLastDaughter(daughter2)
{
  this->SetPdgCode(pdg);

  fP4 = new TLorentzVector(p);
  fX4 = new TLorentzVector(v);

  fRescatterCode  = -1;
  fPolzTheta      = -999; 
  fPolzPhi        = -999;    
  fIsBound        = false;
  fRemovalEnergy  = 0.;
}
//___________________________________________________________________________
// TParticle-like constructor
GHepParticle::GHepParticle(int pdg, GHepStatus_t status,
        int mother1, int mother2, int daughter1, int daughter2,
        double px, double py, double pz, double En,
        double x, double y, double z, double t) :
TObject(),
fStatus(status),
fFirstMother(mother1),
fLastMother(mother2),
fFirstDaughter(daughter1),
fLastDaughter(daughter2)
{
  this->SetPdgCode(pdg);

  fP4 = new TLorentzVector(px,py,pz,En);
  fX4 = new TLorentzVector(x,y,z,t);

  fRescatterCode  = -1;
  fPolzTheta      = -999; 
  fPolzPhi        = -999;   
  fIsBound        = false;
  fRemovalEnergy  = 0.; 
}
//___________________________________________________________________________
// Copy constructor
GHepParticle::GHepParticle(const GHepParticle & particle) :
TObject()
{
  this->Init();
  this->Copy(particle);
}
//___________________________________________________________________________
GHepParticle::GHepParticle(TRootIOCtor*) :
TObject(),
fPdgCode(0),
fStatus(kIStUndefined),
fRescatterCode(-1),
fFirstMother(-1),
fLastMother(-1),
fFirstDaughter(-1),
fLastDaughter(-1),
fP4(0), 
fX4(0),
fPolzTheta(-999.),
fPolzPhi(-999.),
fRemovalEnergy(0),
fIsBound(false)
{

}
//___________________________________________________________________________
GHepParticle::~GHepParticle()
{
  this->CleanUp();
}
//___________________________________________________________________________
string GHepParticle::Name(void) const
{
  this->AssertIsKnownParticle();

  TParticlePDG * p = PDGLibrary::Instance()->Find(fPdgCode);
  return p->GetName();
}
//___________________________________________________________________________
double GHepParticle::Mass(void) const
{
  this->AssertIsKnownParticle();

  TParticlePDG * p = PDGLibrary::Instance()->Find(fPdgCode);
  return p->Mass();
}
//___________________________________________________________________________
double GHepParticle::Charge(void) const
{
  this->AssertIsKnownParticle();

  TParticlePDG * p = PDGLibrary::Instance()->Find(fPdgCode);
  return p->Charge();
}
//___________________________________________________________________________
double GHepParticle::KinE(bool mass_from_pdg) const
{
  if(!fP4) {
    LOG("GHepParticle", pWARN) << "4-momentum not yet set!";
    return 0;
  }

  double En = fP4->Energy();
  double M = ( (mass_from_pdg) ? this->Mass() : fP4->M() );
  double K = En - M;

  K = TMath::Max(K,0.);
  return K;
}
//___________________________________________________________________________
int GHepParticle::Z(void) const
{
// Decoding Z from the PDG code

  if(!pdg::IsIon(fPdgCode)) 
    return -1;
  else
    return pdg::IonPdgCodeToZ(fPdgCode);
}
//___________________________________________________________________________
int GHepParticle::A(void) const
{
// Decoding A from the PDG code

  if(!pdg::IsIon(fPdgCode)) 
    return -1;
  else
    return pdg::IonPdgCodeToA(fPdgCode);
}
//___________________________________________________________________________
TLorentzVector * GHepParticle::GetP4(void) const 
{ 
// see GHepParticle::P4() for a method that does not create a new object and
// transfers its ownership 

  if(fP4) {
     TLorentzVector * p4 = new TLorentzVector(*fP4); 
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
     LOG("GHepParticle", pDEBUG) 
          << "Return vp = " << utils::print::P4AsShortString(p4);
#endif
     return p4;
  } else {
    LOG("GHepParticle", pWARN) << "NULL 4-momentum TLorentzVector";
    return 0;
  }
}
//___________________________________________________________________________
TLorentzVector * GHepParticle::GetX4(void) const 
{ 
// see GHepParticle::X4() for a method that does not create a new object and
// transfers its ownership

  if(fX4) {
     TLorentzVector * x4 = new TLorentzVector(*fX4); 
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
     LOG("GHepParticle", pDEBUG) 
         << "Return x4 = " << utils::print::X4AsString(x4);
#endif
     return x4;
  } else {
    LOG("GHepParticle", pWARN) << "NULL 4-position TLorentzVector";
    return 0;
  }
}
//___________________________________________________________________________
void GHepParticle::SetPdgCode(int code)
{
  fPdgCode = code;
  this->AssertIsKnownParticle();
}
//___________________________________________________________________________
void GHepParticle::SetMomentum(const TLorentzVector & p4)
{
  if(fP4)
      fP4->SetPxPyPzE( p4.Px(), p4.Py(), p4.Pz(), p4.Energy() );
  else
      fP4 = new TLorentzVector(p4);
}
//___________________________________________________________________________
void GHepParticle::SetMomentum(double px, double py, double pz, double En)
{
  if(fP4)
      fP4->SetPxPyPzE(px, py, pz, En);
  else
      fP4 = new TLorentzVector(px, py, pz, En);
}
//___________________________________________________________________________
void GHepParticle::SetPosition(const TLorentzVector & v4)
{
  this->SetPosition(v4.X(), v4.Y(), v4.Z(), v4.T());
}
//___________________________________________________________________________
void GHepParticle::SetPosition(double x, double y, double z, double t)
{
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("GHepParticle", pDEBUG) 
            << "Setting position to (x = " << x << ", y = " 
                               << y << ", z = " << z << ", t = " << t << ")";
#endif

  if(fX4) fX4->SetXYZT(x,y,z,t);
  else    fX4 = new TLorentzVector(x,y,z,t);
}
//___________________________________________________________________________
void GHepParticle::SetEnergy(double En)
{
  this->SetMomentum(this->Px(), this->Py(), this->Pz(), En);
}
//___________________________________________________________________________
void GHepParticle::SetPx(double px)
{
  this->SetMomentum(px, this->Py(), this->Pz(), this->E());
}
//___________________________________________________________________________
void GHepParticle::SetPy(double py)
{
  this->SetMomentum(this->Px(), py, this->Pz(), this->E());
}
//___________________________________________________________________________
void GHepParticle::SetPz(double pz)
{
  this->SetMomentum(this->Px(), this->Py(), pz, this->E());
}
//___________________________________________________________________________
bool GHepParticle::IsOnMassShell(void) const
{
  this->AssertIsKnownParticle();

  TParticlePDG * p = PDGLibrary::Instance()->Find(fPdgCode);

  double Mpdg = p->Mass();
  double M4p  = (fP4) ? fP4->M() : 0.;

//  return utils::math::AreEqual(Mpdg, M4p);

  return (TMath::Abs(M4p-Mpdg) < kOffShellDm);
}
//___________________________________________________________________________
bool GHepParticle::IsOffMassShell(void) const
{
  return (! this->IsOnMassShell());
}
//___________________________________________________________________________
bool GHepParticle::PolzIsSet(void) const
{
// checks whether the polarization angles have been set

  return (fPolzTheta > -999 && fPolzPhi > -999);
}
//___________________________________________________________________________
void GHepParticle::GetPolarization(TVector3 & polz)
{
// gets the polarization vector

  if(! this->PolzIsSet() ) {
      polz.SetXYZ(0.,0.,0.);
      return;
  }
  polz.SetX( TMath::Sin(fPolzTheta) * TMath::Cos(fPolzPhi) );
  polz.SetY( TMath::Sin(fPolzTheta) * TMath::Sin(fPolzPhi) );
  polz.SetZ( TMath::Cos(fPolzTheta) );
}
//___________________________________________________________________________
void GHepParticle::SetPolarization(double theta, double phi)
{
// sets the polarization angles

  if(theta>=0 && theta<=kPi && phi>=0 && phi<2*kPi) 
  {
    fPolzTheta = theta; 
    fPolzPhi   = phi;    

  } else {
    LOG("GHepParticle", pERROR) 
      << "Invalid polarization angles (polar = " << theta 
      << ", azimuthal = " << phi << ")";
  }
}
//___________________________________________________________________________
void GHepParticle::SetPolarization(const TVector3 & polz)
{
// sets the polarization angles

  double p = polz.Mag();
  if(! (p>0) ) {
    LOG("GHepParticle", pERROR) 
           << "Input polarization vector has non-positive norm! Ignoring it";
    return;
  }

  double theta = TMath::ACos(polz.z()/p);
  double phi   = kPi + TMath::ATan2(-polz.y(), -polz.x());

  this->SetPolarization(theta,phi);
}
//___________________________________________________________________________
void GHepParticle::SetBound(bool bound)
{
  // only set it for p or n
  bool is_nucleon = pdg::IsNeutronOrProton(fPdgCode);
  if(!is_nucleon && bound) {
    LOG("GHepParticle", pERROR) 
       << "Refusing to set the bound flag for particles other than nucleons";
    LOG("GHepParticle", pERROR) 
       << "(Requested for pdg = " << fPdgCode << ")";
    return;
  }
  // if the particles isn't bound then make sure that its removal energy = 0
  if(!bound) {
   fRemovalEnergy = 0;
  }
  // set the flag
  fIsBound = bound;
}
//___________________________________________________________________________
void GHepParticle::SetRemovalEnergy(double Erm)
{
  fRemovalEnergy = TMath::Max(Erm, 0.); // non-negative

  // if a value was set, make sure that the IsBound flag is turned on
  if(fRemovalEnergy>0) this->SetBound(true);
}
//___________________________________________________________________________
void GHepParticle::Init(void)
{
  fPdgCode       = 0;
  fStatus        = kIStUndefined;
  fRescatterCode = -1;
  fFirstMother   = -1;
  fLastMother    = -1;
  fFirstDaughter = -1;
  fLastDaughter  = -1;
  fPolzTheta     = -999; 
  fPolzPhi       = -999;    
  fIsBound       = false;
  fRemovalEnergy = 0.;
  fP4            = new TLorentzVector(0,0,0,0);
  fX4            = new TLorentzVector(0,0,0,0);
}
//___________________________________________________________________________
void GHepParticle::CleanUp(void)
{
// deallocate memory

  if(fP4) delete fP4;
  if(fX4) delete fX4;
  fP4 = 0;
  fX4 = 0;
}
//___________________________________________________________________________
void GHepParticle::Reset(void)
{
// deallocate memory + initialize

  this->CleanUp();
  this->Init();
}
//___________________________________________________________________________
void GHepParticle::Clear(Option_t * /*option*/)
{
// implement the Clear(Option_t *) method so that the GHepParticle when is a
// member of a GHepRecord, gets deleted properly when calling TClonesArray's
// Clear("C")

  this->CleanUp();
}
//___________________________________________________________________________
void GHepParticle::Print(ostream & stream) const
{
  stream << "\n |";
  stream << setfill(' ') << setw(14) << this->Name()           << " | ";
  stream << setfill(' ') << setw(14) << this->Pdg()            << " | ";
  stream << setfill(' ') << setw(6)  << this->Status()         << " | ";
  stream << setfill(' ') << setw(3)  << this->FirstMother()    << " | ";
  stream << setfill(' ') << setw(3)  << this->LastMother()     << " | ";
  stream << setfill(' ') << setw(3)  << this->FirstDaughter()  << " | ";
  stream << setfill(' ') << setw(3)  << this->LastDaughter()   << " | ";
  stream << std::fixed  << setprecision(3);
  stream << setfill(' ') << setw(6)  << this->Px()             << " | ";
  stream << setfill(' ') << setw(6)  << this->Py()             << " | ";
  stream << setfill(' ') << setw(6)  << this->Pz()             << " | ";
  stream << setfill(' ') << setw(6)  << this->E()              << " | ";
  stream << setfill(' ') << setw(6)  << this->Mass()           << " | ";
  
  int rescat_code = this->RescatterCode();
  if( rescat_code != -1 ) 
  {
     stream << setfill(' ') << setw(5)  << rescat_code << " | ";
  }
}
//___________________________________________________________________________
void GHepParticle::Print(Option_t * /*opt*/) const
{
// implement the TObject's Print(Option_t *) method

  this->Print(cout);
}
//___________________________________________________________________________
bool GHepParticle::Compare(const GHepParticle * p) const
{
// Do the comparisons in steps & put the ones that cost most
// in the inner-most {}

  bool same_particle = (this->fPdgCode == p->fPdgCode);
  bool same_status   = (this->fStatus  == p->fStatus );

  if( !same_particle || !same_status )  return false;
  else {
     if ( ! this->CompareFamily(p) )    return false;
     else {
       if( ! this->CompareMomentum(p) ) return false;
       else                             return true;
     }
  }
}
//___________________________________________________________________________
bool GHepParticle::ComparePdgCodes(const GHepParticle * p) const
{
  return (this->fPdgCode == p->fPdgCode);
}
//___________________________________________________________________________
bool GHepParticle::CompareStatusCodes(const GHepParticle * p) const
{
  return (this->fStatus == p->fStatus);
}
//___________________________________________________________________________
bool GHepParticle::CompareFamily(const GHepParticle * p) const
{
  bool same_family  = (
          this->fFirstMother == p->fFirstMother &&
                  this->fLastMother  == p->fLastMother &&
                          this->fFirstDaughter == p->fFirstDaughter &&
                                   this->fLastDaughter  == p->fLastDaughter
  );
  return same_family;
}
//___________________________________________________________________________
bool GHepParticle::CompareMomentum(const GHepParticle * p) const
{
  double dE  = TMath::Abs( this->E()  - p->E()  );
  double dPx = TMath::Abs( this->Px() - p->Px() );
  double dPy = TMath::Abs( this->Py() - p->Py() );
  double dPz = TMath::Abs( this->Pz() - p->Pz() );

  bool same_momentum =
       (dE < kPCutOff && dPx < kPCutOff && dPy < kPCutOff && dPz < kPCutOff);

  return same_momentum;
}
//___________________________________________________________________________
void GHepParticle::Copy(const GHepParticle & particle)
{
  this->SetStatus           (particle.Status()          );
  this->SetPdgCode          (particle.Pdg()             );
  this->SetRescatterCode    (particle.RescatterCode()   );
  this->SetFirstMother      (particle.FirstMother()     );
  this->SetLastMother       (particle.LastMother()      );
  this->SetFirstDaughter    (particle.FirstDaughter()   );
  this->SetLastDaughter     (particle.LastDaughter()    );

  this->SetMomentum (*particle.P4());
  this->SetPosition (*particle.X4());

  this->fPolzTheta = particle.fPolzTheta;
  this->fPolzPhi   = particle.fPolzPhi;

  this->fIsBound       = particle.fIsBound;
  this->fRemovalEnergy = particle.fRemovalEnergy;
}
//___________________________________________________________________________
void GHepParticle::AssertIsKnownParticle(void) const
{
  TParticlePDG * p = PDGLibrary::Instance()->Find(fPdgCode);
  if(!p) {
    LOG("GHepParticle", pFATAL)
      << "\n** You are attempting to insert particle with PDG code = " 
      << fPdgCode << " into the event record."
      << "\n** This particle can not be found in "
      << "$GENIE/data/evgen/catalogues/pdg/genie_pdg_table.txt";
    gAbortingInErr = true;
    exit(1);
  }
}
//___________________________________________________________________________
bool GHepParticle::operator == (const GHepParticle & p) const
{
  return (this->Compare(&p));
}
//___________________________________________________________________________
GHepParticle & GHepParticle::operator = (const GHepParticle & p)
{
  this->Copy(p);
  return (*this);
}
//___________________________________________________________________________

