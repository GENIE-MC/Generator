//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ May 05, 2010 - CR
   Adding special ctor for ROOT I/O purposes so as to avoid memory leak due to
   memory allocated in the default ctor when objects of this class are read by
   the ROOT Streamer.

*/
//____________________________________________________________________________

#include <TMath.h>
#include <TLorentzVector.h>
#include <TRootIOCtor.h>

#include "Framework/Interaction/Kinematics.h"
#include "Framework/Messenger/Messenger.h"

using std::endl;

using namespace genie;

ClassImp(Kinematics)

//____________________________________________________________________________
namespace genie {
 ostream & operator << (ostream & stream, const Kinematics & kinematics)
 {
   kinematics.Print(stream);
   return stream;
 }
}
//___________________________________________________________________________
Kinematics::Kinematics() :
TObject()
{
  this->Init();
}
//____________________________________________________________________________
Kinematics::Kinematics(const Kinematics & kinematics) :
TObject()
{
  this->Init();
  this->Copy(kinematics);
}
//____________________________________________________________________________
Kinematics::Kinematics(TRootIOCtor*) :
TObject(),
fP4Fsl(0),
fP4HadSyst(0)
{

}
//____________________________________________________________________________
Kinematics::~Kinematics()
{
  this->CleanUp();
}
//____________________________________________________________________________
void Kinematics::Init(void)
{
  fKV.clear();

  fP4Fsl     = new TLorentzVector;
  fP4HadSyst = new TLorentzVector;
}
//____________________________________________________________________________
void Kinematics::CleanUp(void)
{
  fKV.clear();

  delete fP4Fsl;
  delete fP4HadSyst;
}
//____________________________________________________________________________
void Kinematics::Reset(void)
{
  fKV.clear();

  this->SetFSLeptonP4 (0,0,0,0);
  this->SetHadSystP4  (0,0,0,0);
}
//____________________________________________________________________________
void Kinematics::Copy(const Kinematics & kinematics)
{
  this->Reset();

  map<KineVar_t, double>::const_iterator iter;

  for(iter = kinematics.fKV.begin(); iter != kinematics.fKV.end(); ++iter) {
    KineVar_t kv  = iter->first;
    double    val = iter->second;
    this->SetKV(kv,val);
  }

  this->SetFSLeptonP4 (*kinematics.fP4Fsl);
  this->SetHadSystP4  (*kinematics.fP4HadSyst);
}
//____________________________________________________________________________
double Kinematics::x(bool selected) const
{
// returns the running or selected value of Bjorken scaling variable x

  KineVar_t kvar = (selected) ? kKVSelx : kKVx;

  if(this->KVSet(kvar)) { return this->GetKV(kvar); }
  else {
    LOG("Interaction", pWARN) << "Kinematic variable x was not set";
  }
  return -99999;
}
//____________________________________________________________________________
double Kinematics::y(bool selected) const
{
// returns the running or selected value of inelasticity y

  KineVar_t kvar = (selected) ? kKVSely : kKVy;

  if(this->KVSet(kvar)) { return this->GetKV(kvar); }
  else {
    LOG("Interaction", pWARN) << "Kinematic variable y was not set";
  }
  return -99999;
}
//____________________________________________________________________________
double Kinematics::Q2(bool selected) const
{
// returns the running or selected value of momentum transfer Q2 (>0)

  if(selected) {
    if      (this->KVSet(kKVSelQ2) ) { return     this->GetKV(kKVSelQ2); }
    else if (this->KVSet(kKVSelq2) ) { return -1* this->GetKV(kKVSelq2); }
  } else {
    if      (this->KVSet(kKVQ2) )    { return     this->GetKV(kKVQ2); }
    else if (this->KVSet(kKVq2) )    { return -1* this->GetKV(kKVq2); }
  }

  LOG("Interaction", pWARN) << "Kinematic variable Q2 was not set";
  return -99999;
}
//____________________________________________________________________________
double Kinematics::q2(bool selected) const
{
// returns the running or selected value of momentum transfer q2 (<0)

  if(selected) {
    if      (this->KVSet(kKVSelQ2) ) { return -1* this->GetKV(kKVSelQ2); }
    else if (this->KVSet(kKVSelq2) ) { return     this->GetKV(kKVSelq2); }
  } else {
    if      (this->KVSet(kKVQ2) )    { return -1* this->GetKV(kKVQ2); }
    else if (this->KVSet(kKVq2) )    { return     this->GetKV(kKVq2); }
  }

  LOG("Interaction", pWARN) << "Kinematic variable q2 was not set";
  return -99999;
}
//____________________________________________________________________________
double Kinematics::W(bool selected) const
{
// returns the running or selected value of invariant hadronic mass W

  KineVar_t kvar = (selected) ? kKVSelW : kKVW;

  if(this->KVSet(kvar)) { return this->GetKV(kvar); }
  else {
    LOG("Interaction", pWARN) << "Kinematic variable W was not set";
  }
  return -99999;
}
//____________________________________________________________________________
double Kinematics::t(bool selected) const
{
// returns the running or selected value of invariant hadronic mass W

  KineVar_t kvar = (selected) ? kKVSelt : kKVt;

  if(this->KVSet(kvar)) { return this->GetKV(kvar); }
  else {
    LOG("Interaction", pWARN) << "Kinematic variable t was not set";
  }
  return -99999;
}
//____________________________________________________________________________
double Kinematics::Logx(bool selected) const
{
  double xs = this->x(selected);
  return (xs>0) ? TMath::Log(xs) : -99999;
}
//____________________________________________________________________________
double Kinematics::Logy(bool selected) const
{
  double ys = this->y(selected);
  return (ys>0) ? TMath::Log(ys) : -99999;
}
//____________________________________________________________________________
double Kinematics::LogQ2(bool selected) const
{
  double Q2s = this->Q2(selected);
  return (Q2s>0) ? TMath::Log(Q2s) : -99999;
}
//____________________________________________________________________________
double Kinematics::LogW(bool selected) const
{
  double Ws = this->W(selected);
  return (Ws>0) ? TMath::Log(Ws) : -99999;
}
//____________________________________________________________________________
double Kinematics::Log10x(bool selected) const
{
  double xs = this->x(selected);
  return (xs>0) ? TMath::Log10(xs) : -99999;
}
//____________________________________________________________________________
double Kinematics::Log10y(bool selected) const
{
  double ys = this->y(selected);
  return (ys>0) ? TMath::Log10(ys) : -99999;
}
//____________________________________________________________________________
double Kinematics::Log10Q2(bool selected) const
{
  double Q2s = this->Q2(selected);
  return (Q2s>0) ? TMath::Log10(Q2s) : -99999;
}
//____________________________________________________________________________
double Kinematics::Log10W(bool selected) const
{
  double Ws = this->W(selected);
  return (Ws>0) ? TMath::Log10(Ws) : -99999;
}
//____________________________________________________________________________
void Kinematics::Setx(double xbj, bool selected)
{
// sets the running or selected value of Bjorken scaling variable x

  if(xbj<0 || xbj>1) {
     LOG("Interaction", pWARN)
                      << "Setting unphysical value for x (x = " << xbj << ")";
  }
  KineVar_t kvar = (selected) ? kKVSelx : kKVx;
  this->SetKV(kvar, xbj);
}
//____________________________________________________________________________
void Kinematics::Sety(double inel_y, bool selected)
{
// sets the running or selected value of inelasticity y

  if(inel_y<0 || inel_y>1) {
     LOG("Interaction", pWARN)
                     << "Setting unphysical value for y (y = " << inel_y << ")";
  }
  KineVar_t kvar = (selected) ? kKVSely : kKVy;
  this->SetKV(kvar, inel_y);
}
//____________________________________________________________________________
void Kinematics::SetQ2(double Qsqrd, bool selected)
{
// sets the running or selected value of momentum transfer Q2 (>0)

  if(Qsqrd<0) {
     LOG("Interaction", pWARN)
                 << "Setting unphysical value for Q2 (Q2 = " << Qsqrd << ")";
  }
  KineVar_t kvar = (selected) ? kKVSelQ2 : kKVQ2;
  this->SetKV(kvar, Qsqrd);
}
//____________________________________________________________________________
void Kinematics::Setq2(double qsqrd, bool selected)
{
// sets the running or selected value of momentum transfer q2 (<0)

  if(qsqrd>0) {
     LOG("Interaction", pWARN)
                 << "Setting unphysical value for q2 (q2 = " << qsqrd << ")";
  }
  KineVar_t kvar = (selected) ? kKVSelq2 : kKVq2;
  this->SetKV(kvar, qsqrd);
}
//____________________________________________________________________________
void Kinematics::SetW(double hadr_mass_W, bool selected)
{
// sets the running or selected value of invariant hadronic mass W

  if(hadr_mass_W<0) {
     LOG("Interaction", pWARN)
                   << "Setting unphysical value for W (W = " << hadr_mass_W << ")";
  }
  KineVar_t kvar = (selected) ? kKVSelW : kKVW;
  this->SetKV(kvar, hadr_mass_W);
}
//____________________________________________________________________________
void Kinematics::Sett(double tval, bool selected)
{
  KineVar_t kvar = (selected) ? kKVSelt : kKVt;
  this->SetKV(kvar, tval);
}
//____________________________________________________________________________
void Kinematics::SetFSLeptonP4(const TLorentzVector & p4)
{
  fP4Fsl->SetPxPyPzE(p4.Px(), p4.Py(), p4.Pz(), p4.E());
}
//____________________________________________________________________________
void Kinematics::SetFSLeptonP4(double px, double py, double pz, double E)
{
  fP4Fsl->SetPxPyPzE(px,py,pz,E);
}
//____________________________________________________________________________
void Kinematics::SetHadSystP4(const TLorentzVector & p4)
{
  fP4HadSyst->SetPxPyPzE(p4.Px(), p4.Py(), p4.Pz(), p4.E());
}
//____________________________________________________________________________
void Kinematics::SetHadSystP4(double px, double py, double pz, double E)
{
  fP4HadSyst->SetPxPyPzE(px,py,pz,E);
}
//____________________________________________________________________________
bool Kinematics::KVSet(KineVar_t kv) const
{
  if(fKV.count(kv) == 1) return true;
  else return false;
}
//____________________________________________________________________________
double Kinematics::GetKV(KineVar_t kv) const
{
  if(this->KVSet(kv)) {
     map<KineVar_t, double>::const_iterator iter = fKV.find(kv);
     return iter->second;
  } else {
    LOG("Interaction", pWARN)
        << "Kinematic variable: " << KineVar::AsString(kv) << " was not set";
  }
  return -99999;
}
//____________________________________________________________________________
void Kinematics::SetKV(KineVar_t kv, double value)
{
  LOG("Interaction", pDEBUG)
            << "Setting " << KineVar::AsString(kv) << " to " << value;

  if(this->KVSet(kv)) {
     fKV[kv] = value;
  } else {
     fKV.insert( map<KineVar_t, double>::value_type(kv,value) );
  }
}
//____________________________________________________________________________
void Kinematics::ClearRunningValues(void)
{
// clear the running values (leave the selected ones)
//
  fKV.erase( kKVx  );
  fKV.erase( kKVy  );
  fKV.erase( kKVQ2 );
  fKV.erase( kKVq2 );
  fKV.erase( kKVW  );
  fKV.erase( kKVt  );
}
//____________________________________________________________________________
void Kinematics::UseSelectedKinematics(void)
{
// copy the selected kinematics into the running ones
//
  map<KineVar_t, double>::const_iterator iter;
  iter = fKV.find(kKVSelx);
  if(iter != fKV.end()) this->Setx(iter->second);
  iter = fKV.find(kKVSely);
  if(iter != fKV.end()) this->Sety(iter->second);
  iter = fKV.find(kKVSelQ2);
  if(iter != fKV.end()) this->SetQ2(iter->second);
  iter = fKV.find(kKVSelq2);
  if(iter != fKV.end()) this->Setq2(iter->second);
  iter = fKV.find(kKVSelW);
  if(iter != fKV.end()) this->SetW(iter->second);
  iter = fKV.find(kKVSelt);
  if(iter != fKV.end()) this->Sett(iter->second);
}
//____________________________________________________________________________
void Kinematics::Print(ostream & stream) const
{
  stream << "[-] [Kinematics]" << endl;

  map<KineVar_t, double>::const_iterator iter;

  for(iter = fKV.begin(); iter != fKV.end(); ++iter) {
    KineVar_t kv  = iter->first;
    double    val = iter->second;
    stream << " |--> " << KineVar::AsString(kv) << " = " << val << endl;
  }
}
//____________________________________________________________________________
Kinematics & Kinematics::operator = (const Kinematics & kinematics)
{
  this->Copy(kinematics);
  return (*this);
}
//___________________________________________________________________________
