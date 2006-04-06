//____________________________________________________________________________
/*!

\class    genie::Kinematics

\brief    Kinematic variables for an event

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 08, 2004

*/
//____________________________________________________________________________

#include <TMath.h>

#include "Interaction/Kinematics.h"
#include "Messenger/Messenger.h"

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
Kinematics::Kinematics()
{
  this->Reset();
}
//____________________________________________________________________________
Kinematics::Kinematics(const Kinematics & kinematics)
{
  this->Reset();
  this->Copy(kinematics);
}
//____________________________________________________________________________
Kinematics::~Kinematics()
{
  this->Reset();
}
//____________________________________________________________________________
void Kinematics::Reset(void)
{
  fKV.clear();
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
}
//____________________________________________________________________________
double Kinematics::x(void) const
{
  if(this->KVSet(kKVx)) { return this->GetKV(kKVx); }
  else {
    LOG("Interaction", pWARN) << "Kinematic variable x was not set";
  }
  return -99999;
}
//____________________________________________________________________________
double Kinematics::y(void) const
{
  if(this->KVSet(kKVy)) { return this->GetKV(kKVy); }
  else {
    LOG("Interaction", pWARN) << "Kinematic variable y was not set";
  }
  return -99999;
}
//____________________________________________________________________________
double Kinematics::Q2(void) const
{
  if      (this->KVSet(kKVQ2) ) { return     this->GetKV(kKVQ2); }
  else if (this->KVSet(kKVq2) ) { return -1* this->GetKV(kKVq2); }
  else {
    LOG("Interaction", pWARN) << "Kinematic variable Q2 was not set";
  }
  return -99999;
}
//____________________________________________________________________________
double Kinematics::q2(void) const
{
  if      (this->KVSet(kKVQ2) ) { return -1* this->GetKV(kKVQ2); }
  else if (this->KVSet(kKVq2) ) { return     this->GetKV(kKVq2); }
  else {
    LOG("Interaction", pWARN) << "Kinematic variable q2 was not set";
  }
  return -99999;
}
//____________________________________________________________________________
double Kinematics::W(void) const
{
  if(this->KVSet(kKVW)) { return this->GetKV(kKVW); }
  else {
    LOG("Interaction", pWARN) << "Kinematic variable W was not set";
  }
  return -99999;
}
//____________________________________________________________________________
double Kinematics::Logx(void) const
{
  double x = this->x();
  return (x>0) ? TMath::Log(x) : -99999;
}
//____________________________________________________________________________
double Kinematics::Logy(void) const
{
  double y = this->y();
  return (y>0) ? TMath::Log(y) : -99999;
}
//____________________________________________________________________________
double Kinematics::LogQ2(void) const
{
  double Q2 = this->Q2();
  return (Q2>0) ? TMath::Log(Q2) : -99999;
}
//____________________________________________________________________________
double Kinematics::LogW(void) const
{
  double W = this->W();
  return (W>0) ? TMath::Log(W) : -99999;
}
//____________________________________________________________________________
double Kinematics::Log10x(void) const
{
  double x = this->x();
  return (x>0) ? TMath::Log10(x) : -99999;
}
//____________________________________________________________________________
double Kinematics::Log10y(void) const
{
  double y = this->y();
  return (y>0) ? TMath::Log10(y) : -99999;
}
//____________________________________________________________________________
double Kinematics::Log10Q2(void) const
{
  double Q2 = this->Q2();
  return (Q2>0) ? TMath::Log10(Q2) : -99999;
}
//____________________________________________________________________________
double Kinematics::Log10W(void) const
{
  double W = this->W();
  return (W>0) ? TMath::Log10(W) : -99999;
}
//____________________________________________________________________________
void Kinematics::Setx(double x)
{
  if(x<0 || x>1) {
     LOG("Interaction", pWARN)
                      << "Setting unphysical value for x (x = " << x << ")";
  }
  this->SetKV(kKVx, x);
}
//____________________________________________________________________________
void Kinematics::Sety(double y)
{
  if(y<0 || y>1) {
     LOG("Interaction", pWARN)
                     << "Setting unphysical value for y (y = " << y << ")";
  }
  this->SetKV(kKVy, y);
}
//____________________________________________________________________________
void Kinematics::SetQ2(double Q2)
{
  if(Q2<0) {
     LOG("Interaction", pWARN)
                 << "Setting unphysical value for Q2 (Q2 = " << Q2 << ")";
  }
  this->SetKV(kKVQ2, Q2);
}
//____________________________________________________________________________
void Kinematics::Setq2(double q2)
{
  if(q2>0) {
     LOG("Interaction", pWARN)
                 << "Setting unphysical value for q2 (q2 = " << q2 << ")";
  }
  this->SetKV(kKVq2, q2);
}
//____________________________________________________________________________
void Kinematics::SetW(double W)
{
  if(W<0) {
     LOG("Interaction", pWARN)
                   << "Setting unphysical value for W (W = " << W << ")";
  }
  this->SetKV(kKVW, W);
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

