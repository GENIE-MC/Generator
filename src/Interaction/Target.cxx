//____________________________________________________________________________
/*!

\class    genie::Target

\brief    A Neutrino Interaction Target. Is a transparent encapsulation of 
          quite different physical systems such as a nuclear target, a
          'spectator' nuclear target with a struck nucleon, a free nucleon or
          a free particle (eg a e- target in the inverse muon decay reaction)
          
\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 03, 2004

*/
//____________________________________________________________________________

#include <TParticlePDG.h>

#include "Conventions/Constants.h"
#include "Interaction/Target.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGLibrary.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "Utils/PrintUtils.h"

using namespace genie;
using namespace genie::constants;

using std::endl;

//____________________________________________________________________________
namespace genie {
 ostream & operator<< (ostream& stream, const Target & target)
 {
   target.Print(stream);

   return stream;
 }
}
//___________________________________________________________________________
Target::Target()
{
  Init();
}
//___________________________________________________________________________
Target::Target(int pdgc)
{
  Init();

  fTgtPDG = pdgc;

  if( pdg::IsIon(pdgc) ) {

     int Z = pdg::IonPdgCodeToZ(pdgc);
     int A = pdg::IonPdgCodeToA(pdgc);

     // set Z,A & fix struck nucleon PDG if tgt = free nucleon

     this->SetZA(Z,A);     
  }  
}
//___________________________________________________________________________
Target::Target(int Z, int A)
{
  Init();

  // set Z,A & fix struck nucleon PDG if tgt = free nucleon
  
  fTgtPDG = pdg::IonPdgCode(A,Z);

  this->SetZA(Z,A);  
}
//___________________________________________________________________________
Target::Target(int Z, int A, int struck_nucleon_pdgc)
{
  Init();

  fZ = Z;
  fA = A;

  fTgtPDG = pdg::IonPdgCode(A,Z);
  
  this->ForceNucleusValidity(); // search for this nucleus at the PDG Ions

  this->SetStruckNucleonPDGCode(struck_nucleon_pdgc);  
}
//___________________________________________________________________________
Target::Target(const Target & tgt)
{
  this->Copy(tgt);
}
//___________________________________________________________________________
Target::~Target()
{
  if(fStruckNucP4) delete fStruckNucP4;
}
//___________________________________________________________________________
void Target::Copy(const Target & tgt)
{
  Init();

  fZ = tgt.fZ; // copy A,Z
  fA = tgt.fA;

  fStruckNucPDG = tgt.fStruckNucPDG; // copy struck nucleon & tgt PDG
  fTgtPDG       = tgt.fTgtPDG;

  if(tgt.fStruckNucP4) {

     fStruckNucP4 = new TLorentzVector(*tgt.fStruckNucP4);
  }

  ForceNucleusValidity(); // search for this nucleus at the isotopes chart
  ForceStruckNucleonValidity(); // must be p or n
}
//___________________________________________________________________________
int Target::Z(void) const
{
  return fZ;
}
//___________________________________________________________________________
int Target::N(void) const
{
  return (fA-fZ);
}
//___________________________________________________________________________
int Target::A(void) const
{
  return fA;
}
//___________________________________________________________________________
int Target::PDGCode(void) const
{
  return fTgtPDG;
}
//___________________________________________________________________________
bool Target::IsFreeNucleon(void) const
{
  return ( fA == 1 && (fZ == 0 || fZ == 1) );
}
//___________________________________________________________________________
bool Target::IsProton(void) const
{
  return ( fA == 1 && fZ == 1 );
}
//___________________________________________________________________________
bool Target::IsNeutron(void) const
{
  return ( fA == 1 && fZ == 0 );
}
//___________________________________________________________________________
bool Target::IsNucleus(void) const
{
  return ( fA > 1 ); // IsValidNucleus() was ensured when A,Z were set
}
//___________________________________________________________________________
bool Target::IsParticle(void) const
{
  TParticlePDG * p = PDGLibrary::Instance()->Find(fTgtPDG);

  return (p && fA==0 && fZ==0);
}
//___________________________________________________________________________
bool Target::StruckNucleonIsSet(void) const
{
  return pdg::IsNeutronOrProton(fStruckNucPDG);
}
//___________________________________________________________________________
void Target::SetZA(int Z, int A)
{
  fZ = Z;
  fA = A;

  this->ForceNucleusValidity(); // search at the isotopes chart

  // if the target is a free nucleon, then the struck nucleon pdg code is
  // automaticaly set

  if( this->IsFreeNucleon() ) {

    if( this->IsProton() ) this->SetStruckNucleonPDGCode(kPdgProton);
    else                   this->SetStruckNucleonPDGCode(kPdgNeutron);
  }
}
//___________________________________________________________________________
void Target::SetStruckNucleonPDGCode(int nucl_pdgc)
{
  fStruckNucPDG = nucl_pdgc;

  bool is_valid = this->ForceStruckNucleonValidity();  // must be p or n

  // If it is a valid struck nucleon pdg code, initialize its 4P:
  // at-rest + on-mass-shell
  
  if( is_valid ) {

    double M = PDGLibrary::Instance()->Find(nucl_pdgc)->Mass();

    TLorentzVector p4(0,0,0,M);

    this->SetStruckNucleonP4(p4);
  }
}
//___________________________________________________________________________
void Target::SetStruckNucleonP4(const TLorentzVector & p4)
{
  if(fStruckNucP4) delete fStruckNucP4;

  fStruckNucP4 = new TLorentzVector(p4);
}
//___________________________________________________________________________
int Target::StruckNucleonPDGCode(void) const
{
  return fStruckNucPDG;
}
//___________________________________________________________________________
TLorentzVector * Target::StruckNucleonP4(void) const
{
  if(!fStruckNucP4) {
    LOG("Target", pWARN) << "Returning NULL struck nucleon 4-momentum";
    return 0;
  }

  return fStruckNucP4;
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
double Target::Charge(void) const
{
// Shortcut for commonly used code for extracting the nucleus charge from PDG
//
  TParticlePDG * p = PDGLibrary::Instance()->Find(fTgtPDG);

  if(p) return p->Charge() / 3.; // in +e

  return 0;
}
//___________________________________________________________________________
double Target::StruckNucleonMass(void) const
{
  if(!fStruckNucPDG) {
    LOG("Target", pWARN) << "Returning struck nucleon mass = 0";
    return 0;
  }

  return PDGLibrary::Instance()->Find(fStruckNucPDG)->Mass();
}
//___________________________________________________________________________
void Target::Init(void)
{
  fZ = 0;
  fA = 0;

  fTgtPDG       = 0;
  fStruckNucPDG = 0;
  
  fStruckNucP4  = new TLorentzVector(0,0,0,kNucleonMass);
}
//___________________________________________________________________________
bool Target::ForceStruckNucleonValidity(void)
{
// resets the struck nucleon pdg-code if it found not to be a valid one
  
  bool is_valid =
               pdg::IsProton(fStruckNucPDG) || pdg::IsNeutron(fStruckNucPDG);

  if( ! is_valid ) {

    LOG("Target", pWARN)<< "Reseting to struck nucleon to 'Rootino'";

    fStruckNucPDG = 0;
  }

  return is_valid;
}
//___________________________________________________________________________
void Target::ForceNucleusValidity(void)
{
  if( ! this->IsValidNucleus() ) {

    LOG("Target", pWARN) << "Invalid target -- Reseting to Z = 0, A = 0";

    fZ = 0;
    fA = 0;
  }
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
double Target::BindEnergy(void) const
{
  if( this->IsNucleus() ) {

    // Compute the average binding energy using the
    // semi-empirical, hardcoded, formula from Wapstra (Handbuch der
    // Physik, XXXVIII/1)

    // Eventually, this piece of code should not be here. The Target object
    // should either ask an external calculator or lookup an external table.
  
    double a = 15.835;
    double b =  18.33;
    double s =  23.20;
    double d =   0.714;

    double delta = 0;                        /*E-O*/
    if ( this->IsOddOdd()   ) delta =  11.2; /*O-O*/
    if ( this->IsEvenEven() ) delta = -11.2; /*E-E*/

    double N = (double) this->N();
    double Z = (double) this->Z();
    double A = (double) this->A();
  
    double BE = a*A - b*pow(A,0.667) - s*pow(N-Z,2.0)/A -
                           d*pow(Z,2.0)/pow(A,0.333) - delta/sqrt(A); // MeV

    return ( 1e-3 * BE ); // GeV

  } else return 0;  
}
//___________________________________________________________________________
double Target::BindEnergyPerNucleon(void) const
{
  if( this->IsNucleus() ) {
  
     return ( this->BindEnergy() / this->A() );

  } else return 0;      
}
//___________________________________________________________________________
double Target::BindEnergyLastNucleon(void) const
{
  if( this->IsNucleus() ) {

     //-- temporarily, return the binding energy per nucleon rather than the
     //   separation energy of the last nucleon
  
     return ( this->BindEnergy() / this->N() );

  } else return 0;   
}
//___________________________________________________________________________
bool Target::IsEvenEven(void) const
{
  if( this->IsNucleus() ) {

    int N = this->N();
    int Z = this->Z();

    if( N % 2 == 0 && Z % 2 == 0 ) return true;
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

    int N = this->N();
    int Z = this->Z();

    if( N % 2 == 1 && Z % 2 == 1 ) return true;
  }
  return false;
}
//___________________________________________________________________________
void Target::Print(ostream & stream) const
{
  stream << " target PDG code = " << fTgtPDG << endl;

  if( this->IsNucleus() || this->IsFreeNucleon() ) {
      stream << " Z = " << fZ << ", A = " << fA << endl;
  }
  
  if( this->StruckNucleonIsSet() ) {

    TParticlePDG * p = PDGLibrary::Instance()->Find(fStruckNucPDG);

    stream << " struck nucleon = " << p->GetName() 
               << ", P4 = " << print_utils::P4AsString(fStruckNucP4) << endl;
  }
}
//___________________________________________________________________________


