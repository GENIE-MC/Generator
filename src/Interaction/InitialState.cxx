//____________________________________________________________________________
/*!

\class    genie::InitialState

\brief    Initial State information

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 02, 2004
 
*/
//____________________________________________________________________________

#include <cassert>
#include <iomanip>

#include "Interaction/InitialState.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGLibrary.h"
#include "PDG/PDGUtils.h"

using namespace genie;

using std::endl;
using std::setprecision;
using std::setw;

//____________________________________________________________________________
namespace genie {
 ostream & operator<< (ostream& stream, const InitialState & init_state)
 {
   init_state.Print(stream);

   return stream;
 }
}
//___________________________________________________________________________
InitialState::InitialState()
{
  this->Initialize();
}
//___________________________________________________________________________
InitialState::InitialState(const Target & tgt, int probe_pdgc)
{
  this->Initialize();

  fTarget = new Target(tgt);     // set Target properties
  fProbePdgC = probe_pdgc;       // set Probe PDG code

  // set default, on-mass-shell 4-momenta

  double m = PDGLibrary::Instance() -> Find (  probe_pdgc   ) -> Mass();
  double M = PDGLibrary::Instance() -> Find ( tgt.PDGCode() ) -> Mass();
  
  fProbeP4    = new TLorentzVector(0, 0, 0, m);
  fTargetP4   = new TLorentzVector(0, 0, 0, M);
}
//___________________________________________________________________________
InitialState::InitialState(const InitialState & init_state)
{
  this->Initialize();
  
  this->Copy(init_state);
}
//___________________________________________________________________________
InitialState::~InitialState()
{
  delete fTarget;
  delete fProbeP4;
  delete fTargetP4;
}
//___________________________________________________________________________
void InitialState::Initialize(void)
{
  fProbePdgC  = 0;
  fTarget     = new Target();
  fProbeP4    = new TLorentzVector(0, 0, 0, 0);
  fTargetP4   = new TLorentzVector(0, 0, 0, 0);
}
//___________________________________________________________________________
void InitialState::Copy(const InitialState & init_state)
{
  fProbePdgC  = init_state.fProbePdgC;
  
  fTarget->Copy( *init_state.fTarget );

  double px = init_state.fProbeP4->Px();
  double py = init_state.fProbeP4->Py();
  double pz = init_state.fProbeP4->Pz();
  double E  = init_state.fProbeP4->Energy();
  
  fProbeP4->SetXYZT(px, py, pz, E);
  
  px = init_state.fTargetP4->Px();
  py = init_state.fTargetP4->Py();
  pz = init_state.fTargetP4->Pz();
  E  = init_state.fTargetP4->Energy();

  fTargetP4->SetXYZT(px, py, pz, E);
}
//___________________________________________________________________________
TParticlePDG * InitialState::GetProbe(void) const
{
  TParticlePDG * p = PDGLibrary::Instance()->Find(fProbePdgC);

  return p;
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
void InitialState::SetTargetP4(const TLorentzVector & P4)
{
  fTargetP4 -> SetE  ( P4.E()  );
  fTargetP4 -> SetPx ( P4.Px() );
  fTargetP4 -> SetPy ( P4.Py() );
  fTargetP4 -> SetPz ( P4.Pz() );
}
//___________________________________________________________________________
TLorentzVector * InitialState::GetTargetP4(RefFrame_t ref_frame) const
{
// Return the target 4-momentum in the specified reference frame
// Note: the caller adopts the TLorentzVector object

  switch (ref_frame) {

       //------------------ CENTER OF MOMENTUM FRAME:
       case (kRfCenterOfMass) :
             return 0;
             break;

       //------------------ NUCLEAR TARGET REST FRAME:
       case (kRfTargetAtRest) :
       {
             // for now make sure that the target nucleus is always at 
             // rest and it is only the struck nucleons that can move:
             // so the [target rest frame] = [LAB frame]

             return this->GetTargetP4(kRfLab);
       }
       
       //------------------ STRUCK NUCLEON REST FRAME:
       case (kRfStruckNucAtRest) :
       {
             // make sure that 'struck nucleon' properties were set in
             // the nuclear target object

             assert( fTarget->StruckNucleonP4() != 0 );

             TLorentzVector * pnuc4 = fTarget->StruckNucleonP4();

             // compute velocity vector (px/E, py/E, pz/E)

             double bx = pnuc4->Px() / pnuc4->Energy();
             double by = pnuc4->Py() / pnuc4->Energy();
             double bz = pnuc4->Pz() / pnuc4->Energy();

             // BOOST

             TLorentzVector * p4 = new TLorentzVector(*fTargetP4);

             p4->Boost(-bx,-by,-bz);

             return p4;

             break;
       }
       //------------------ LAB:
       case (kRfLab) :
       {
             TLorentzVector * p4 = new TLorentzVector(*fTargetP4);
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
    
       //------------------ CENTER OF MOMENTUM FRAME:
       case (kRfCenterOfMass) :
             return 0;
             break;
             
       //------------------ NUCLEAR TARGET REST FRAME:
       case (kRfTargetAtRest) :     
       {
             // for now assure that the target nucleus is always at rest
             // and it is only the struck nucleons that can move:
             // so the [target rest frame] = [LAB frame]
       
             TLorentzVector * p4 = new TLorentzVector(*fProbeP4);
             return p4;
       }
       
       //------------------ STRUCK NUCLEON REST FRAME:
       case (kRfStruckNucAtRest) :
       {
             // make sure that 'struck nucleon' properties were set in
             // the nuclear target object
             
             assert( fTarget->StruckNucleonP4() != 0 );

             TLorentzVector * pnuc4 = fTarget->StruckNucleonP4();

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
void InitialState::Print(ostream & stream) const
{
  stream << "[-] [Init-State] " << endl;

  stream << " |--> probe       : " 
         << "PDG-code = " << fProbePdgC 
         << " (" << this->GetProbe()->GetName() << ")" << endl;

  stream << " |--> target      : " 
         << "Z = "          << fTarget->Z() 
         << ", A = "        << fTarget->A() 
         << ", PDG-Code = " << fTarget->PDGCode(); 

  TParticlePDG * tgt = PDGLibrary::Instance()->Find( fTarget->PDGCode() );

  if(tgt) {
    stream << " (" << tgt->GetName() << ")";
  }
  stream << endl;

  stream << " |--> hit nucleon : ";

  int nuc_pdgc = fTarget->StruckNucleonPDGCode();

  if ( pdg::IsNeutronOrProton(nuc_pdgc) ) {

    TParticlePDG * p = PDGLibrary::Instance()->Find(nuc_pdgc);

    stream << "PDC-Code = " << nuc_pdgc << " (" << p->GetName() << ")";
  }
  stream << endl;

  stream << " |--> probe 4P    : " 
         << "(E = "   << setw(12) << setprecision(6) << fProbeP4->E()
         << ", Px = " << setw(12) << setprecision(6) << fProbeP4->Px()
         << ", Py = " << setw(12) << setprecision(6) << fProbeP4->Py()
         << ", Pz = " << setw(12) << setprecision(6) << fProbeP4->Pz() 
         << ")"
         << endl;
  stream << " |--> target 4P   : " 
         << "(E = "   << setw(12) << setprecision(6) << fTargetP4->E()
         << ", Px = " << setw(12) << setprecision(6) << fTargetP4->Px()
         << ", Py = " << setw(12) << setprecision(6) << fTargetP4->Py()
         << ", Pz = " << setw(12) << setprecision(6) << fTargetP4->Pz() 
         << ")"
         << endl;

  if ( pdg::IsNeutronOrProton(nuc_pdgc) ) {

    TLorentzVector * nuc_p4 = fTarget->StruckNucleonP4();

    stream << " |--> nucleon 4P  : " 
           << "(E = "   << setw(12) << setprecision(6) << nuc_p4->E()
           << ", Px = " << setw(12) << setprecision(6) << nuc_p4->Px()
           << ", Py = " << setw(12) << setprecision(6) << nuc_p4->Py()
           << ", Pz = " << setw(12) << setprecision(6) << nuc_p4->Pz() 
           << ")";
  } 
}  
//___________________________________________________________________________  


