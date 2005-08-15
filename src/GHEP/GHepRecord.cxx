//____________________________________________________________________________
/*!

\class   genie::GHepRecord

\brief   Generated Event Record: STDHEP-like record and Summary Information.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created October 1, 2004

*/
//____________________________________________________________________________

#include <cassert>
#include <string>
#include <algorithm>
#include <vector>

#include <TLorentzVector.h>

#include "GHEP/GHepParticle.h"
#include "GHEP/GHepRecord.h"
#include "GHEP/GHepStatus.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGUtils.h"

using std::string;
using std::vector;
using namespace genie;

ClassImp(GHepRecord)

//___________________________________________________________________________
GHepRecord::GHepRecord() :
TClonesArray("genie::GHepParticle")
{
  this->InitGHepRecord();
}
//___________________________________________________________________________
GHepRecord::GHepRecord(int size) :
TClonesArray("genie::GHepParticle", size)
{
  this->InitGHepRecord();
}
//___________________________________________________________________________
GHepRecord::GHepRecord(const GHepRecord & record) :
TClonesArray("genie::GHepParticle", record.GetEntries())
{
  this->InitGHepRecord();
  this->Copy(record);
}
//___________________________________________________________________________
GHepRecord::~GHepRecord()
{
  if(fInteraction) delete fInteraction;
}
//___________________________________________________________________________
Interaction * GHepRecord::GetInteraction(void) const
{
  if(!fInteraction)
  {
    LOG("GHEP", pWARN) << "Returning NULL interaction";
  }
  return fInteraction;
}
//___________________________________________________________________________
void GHepRecord::AttachInteraction(Interaction * interaction)
{
  fInteraction = interaction;
}
//___________________________________________________________________________
GHepParticle * GHepRecord::GetParticle(int position) const
{
// Returns the GHepParticle from the specified position of the event record.

  if( position >=0 && position < this->GetEntries() ) {

     GHepParticle * particle = (GHepParticle *) (*this)[position];
     if(particle) return particle;
  }
  LOG("GHEP", pWARN)
           << "\n **** Returning NULL GHepParticle: "
                                     << "None found at position" << position;
  return 0;
}
//___________________________________________________________________________
GHepParticle * GHepRecord::FindParticle(int pdg, int status, int start) const
{
// Returns the first GHepParticle with the input pdg-code and status
// starting from the specified position of the event record.

  for(int i = start; i < this->GetEntries(); i++) {

     GHepParticle * p = (GHepParticle *) (*this)[i];
     if(p->Status() == status && p->PdgCode() == pdg) return p;
  }
  LOG("GHEP", pWARN)
        << "\n **** Returning NULL GHepParticle: "
              << "None found with pdgc = " << pdg << " and status = "
                                << status << " with position >= " << start;

  LOG("GHEP", pWARN) << "Returning NULL GHepParticle";

  return 0;
}
//___________________________________________________________________________
int GHepRecord::ParticlePosition(int pdg, int status, int start) const
{
// Returns the position of the first GHepParticle with the input pdg-code
// and status starting from the specified position of the event record.

  for(int i = start; i < this->GetEntries(); i++) {

     GHepParticle * p = (GHepParticle *) (*this)[i];
     if(p->Status() == status && p->PdgCode() == pdg) return i;
  }
  LOG("GHEP", pWARN) << "Returning Invalid StdHep Record position";

  return -1;
}
//___________________________________________________________________________
int GHepRecord::ParticlePosition(GHepParticle * particle, int start) const
{
// Returns the position of the first match with the specified GHepParticle
// starting from the specified position of the event record.

  for(int i = start; i < this->GetEntries(); i++) {

     GHepParticle * p = (GHepParticle *) (*this)[i];
     if( p->Compare(particle) ) return i;
  }
  LOG("GHEP", pWARN) << "Returning Invalid StdHep Record position";

  return -1;
}
//___________________________________________________________________________
void GHepRecord::ShiftVertex(const TLorentzVector & vec4)
{
// Shifts the event record entries in space (used when events generated at
// the origin (0,0,0,0) are distributed within the volume described by an
// input ROOT geometry)

  double x0 = vec4.X();
  double y0 = vec4.Y();
  double z0 = vec4.Z();
  double t0 = vec4.T();

  for(int i = 0; i < this->GetEntries(); i++) {

     GHepParticle * p = (GHepParticle *) (*this)[i];

     double vx = x0 + p->Vx();
     double vy = y0 + p->Vy();
     double vz = z0 + p->Vz();
     double vt = t0 + p->Vt();

     p->SetVertex(vx, vy, vz, vt);
  }
}
//___________________________________________________________________________
void GHepRecord::AddParticle(const GHepParticle & p)
{
// Provides a simplified method for inserting entries in the TClonesArray

  unsigned int pos = this->GetEntries();
  LOG("GHEP", pINFO)
    << "Adding particle with pdgc = " << p.PdgCode() << " at slot = " << pos;

  new ((*this)[pos]) GHepParticle(p);

  // Update the mother's daughter list. If the newly inserted particle broke
  // compactification, then run CompactifyDaughterLists()
  this->UpdateDaughterLists();
}
//___________________________________________________________________________
void GHepRecord::AddParticle(
  int pdg, int status, int mom1, int mom2, int dau1, int dau2,
                         const TLorentzVector & p, const TLorentzVector & v)
{
// Provides a 'simplified' method for inserting entries in the TClonesArray

  unsigned int pos = this->GetEntries();
  LOG("GHEP", pINFO)
           << "Adding particle with pdgc = " << pdg << " at slot = " << pos;

  new ((*this)[pos]) GHepParticle(pdg,status, mom1,mom2,dau1,dau2, p, v);

  // Update the mother's daughter list. If the newly inserted particle broke
  // compactification, then run CompactifyDaughterLists()
  this->UpdateDaughterLists();
}
//___________________________________________________________________________
void GHepRecord::AddParticle(
  int pdg, int status, int mom1, int mom2, int dau1, int dau2,
                         double px, double py, double pz, double E,
                                     double x, double y, double z, double t)
{
// Provides a 'simplified' method for inserting entries in the TClonesArray

  unsigned int pos = this->GetEntries();
  LOG("GHEP", pINFO)
           << "Adding particle with pdgc = " << pdg << " at slot = " << pos;

  new ( (*this)[pos] ) GHepParticle (
            pdg, status, mom1, mom2, dau1, dau2, px, py, pz, E, x, y, z, t);

  // Update the mother's daughter list. If the newly inserted particle broke
  // compactification, then run CompactifyDaughterLists()
  this->UpdateDaughterLists();
}
//___________________________________________________________________________
void GHepRecord::UpdateDaughterLists(void)
{
  int pos = this->GetEntries() - 1; // position of last entry

  LOG("GHEP", pINFO)
     << "Updating the daughter-list for the mother of particle at: " << pos;

  GHepParticle * p = this->GetParticle(pos);
  assert(p);

  int mom_pos = p->FirstMother();
  LOG("GHEP", pINFO) << "Mother particle is at slot: " << mom_pos;
  GHepParticle * mom = this->GetParticle(mom_pos);
  if(!mom) return; // may not have mom (eg init state)

  int dau1 = mom->FirstDaughter();
  int dau2 = mom->LastDaughter();

  // handles the case where the daughter list was initially empty
  if(dau1 == -1) {
     mom->SetFirstDaughter(pos);
     mom->SetLastDaughter(pos);
     LOG("GHEP", pINFO)
        << "Done! Daughter-list is compact: [" << pos << ", " << pos << "]";
     return;
  }
  // handles the case where the new daughter is added at the slot just before
  // an already compact daughter list
  if(pos == dau1-1) {
     mom->SetFirstDaughter(pos);
     LOG("GHEP", pINFO)
       << "Done! Daughter-list is compact: [" << pos << ", " << dau2 << "]";
     return;
  }
  // handles the case where the new daughter is added at the slot just after
  // an already compact daughter list
  if(pos == dau2+1) {
     mom->SetLastDaughter(pos);
     LOG("GHEP", pINFO)
       << "Done! Daughter-list is compact: [" << dau1 << ", " << pos << "]";
     return;
  }

  // If you are here, then the last particle insertion broke the daughter
  // list compactification - Run the compactifier
  LOG("GHEP", pINFO)
                   << "Daughter-list is not compact - Running compactifier";
  this->CompactifyDaughterLists();
}
//___________________________________________________________________________
void GHepRecord::CompactifyDaughterLists(void)
{
  int n     = this->GetEntries();
  int start = this->FirstNonInitStateEntry();

  for(int i=0; i<n; i++) {
     bool compact = this->HasCompactDaughterList(i);
     if(!compact) {
        int ndau = 0;     // number of daughters
        int dau1 = start; // 1st daughter position
        int k=start+1;
        for(; k<n; k++) {
            GHepParticle * p = this->GetParticle(k);
            if(p->FirstMother() == i) {
                 ndau++;
                 this->SwapParticles(start,k);
                 start++;
            }
        }
        if(ndau>0) {
          this->GetParticle(i)->SetFirstDaughter(dau1);
          this->GetParticle(i)->SetLastDaughter(dau1+ndau);
        } else {
          this->GetParticle(i)->SetFirstDaughter(-1);
          this->GetParticle(i)->SetLastDaughter(-1);
        }
     } //!compact
     LOG("GHEP", pINFO)
       << "Compactifying daughter-list of particle at: " << i << " - Done!";
  }
  this->FinalizeDaughterLists();
}
//___________________________________________________________________________
bool GHepRecord::HasCompactDaughterList(int pos)
{
  LOG("GHEP", pDEBUG) << "Examining daughter-list of particle at: " << pos;

  vector<int> daughters;
  GHepParticle * p = 0;
  TIter iter(this);
  int i=0;
  while( (p = (GHepParticle *)iter.Next()) ) {
    if(p->FirstMother() == pos) {
       LOG("GHEP", pDEBUG) << "Particle at: " << i << " is a daughter";
       daughters.push_back(i);
    }
    i++;
  }

  bool is_compact = true;
  if(daughters.size()>1) {
    sort(daughters.begin(), daughters.end());
    vector<int>::iterator diter = daughters.begin();
    int prev = *diter;
    for(; diter != daughters.end(); ++diter) {
       int curr = *diter;
       is_compact = is_compact && (TMath::Abs(prev-curr)<=1);
       prev = curr;
    }
  }
  LOG("GHEP", pINFO)
         << "Daughter-list of particle at: " << pos << " is "
                                << (is_compact ? "" : "not") << " compact";
  return is_compact;
}
//___________________________________________________________________________
int GHepRecord::FirstNonInitStateEntry(void)
{
  GHepParticle * p = 0;
  TIter iter(this);
  int pos = 0;
  while( (p = (GHepParticle *)iter.Next()) ) {
    int ist = p->Status();
    if(ist != kIStInitialState && ist != kIstNucleonTarget) return pos;
    pos++;
  }
  return pos;
}
//___________________________________________________________________________
void GHepRecord::SwapParticles(int i, int j)
{
  LOG("GHEP", pINFO) << "Swapping GHepParticles : " << i << " <--> " << j;

  int n = this->GetEntries();
  assert(i>=0 && j>=0 && i<n && j<n);

  if(i==j) return;

  GHepParticle * pi  = this->GetParticle(i);
  GHepParticle * pj  = this->GetParticle(j);
  GHepParticle * tmp = new GHepParticle(*pi);

  pi->Copy(*pj);
  pj->Copy(*tmp);

  delete tmp;

  // tell their daughters
  if(pi->HasDaughters()) {
    for(int k=pi->FirstDaughter(); k<=pi->LastDaughter(); k++)
      this->GetParticle(k)->SetFirstMother(j);
  }
  if(pj->HasDaughters()) {
    for(int k=pj->FirstDaughter(); k<=pj->LastDaughter(); k++)
      this->GetParticle(k)->SetFirstMother(i);
  }
}
//___________________________________________________________________________
void GHepRecord::FinalizeDaughterLists(void)
{
// Update all daughter-lists based on particle 'first mother' field.
// To work correctly, the daughter-lists must have been compactified first.

  GHepParticle * p1 = 0;
  TIter iter1(this);
  int i1=0;
  while( (p1 = (GHepParticle *)iter1.Next()) ) {
    int dau1 = -1;
    int dau2 = -1;
    GHepParticle * p2 = 0;
    TIter iter2(this);
    int i2=0;
    while( (p2 = (GHepParticle *)iter2.Next()) ) {

       if(p2->FirstMother() == i1) {
          dau1 = (dau1<0) ? i2 : TMath::Min(dau1,i2);
          dau2 = (dau2<0) ? i2 : TMath::Max(dau2,i2);
       }
       i2++;
    }
    i1++;
    p1 -> SetFirstDaughter (dau1);
    p1 -> SetLastDaughter  (dau2);
  }
}
//___________________________________________________________________________
void GHepRecord::SwitchIsPauliBlocked(bool on_off)
{
  string status = (on_off) ? "ON" : "OFF";

  LOG("GHEP", pINFO) << "Switching Pauli Blocking flag: " << status;

  fEventIsPauliBlocked = on_off;
}
//___________________________________________________________________________
bool GHepRecord::IsForbidden(void) const
{
// Summarizes record flags (only one flag so far)

  return fEventIsPauliBlocked;
}
//___________________________________________________________________________
void GHepRecord::InitGHepRecord(void)
{
  fInteraction = 0;

  fEventIsPauliBlocked = false;
}
//___________________________________________________________________________
void GHepRecord::ResetGHepRecord(void)
{
  if (fInteraction) delete fInteraction;

  this->Clear("C");
  this->InitGHepRecord();
}
//___________________________________________________________________________
void GHepRecord::Copy(const GHepRecord & record)
{
  // clean up
  this->ResetGHepRecord();

  // copy event record entries
  unsigned int ientry = 0;
  GHepParticle * p = 0;
  TIter ghepiter(&record);
  while ( (p = (GHepParticle *) ghepiter.Next()) )
                         new ( (*this)[ientry++] ) GHepParticle(*p);

  // copy summary
  fInteraction = new Interaction( *record.fInteraction );

  // copy flags
  fEventIsPauliBlocked = record.fEventIsPauliBlocked;
}
//___________________________________________________________________________
