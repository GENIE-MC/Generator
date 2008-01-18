//____________________________________________________________________________
/*
 Copyright (c) 2003-2008, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         STFC, Rutherford Appleton Laboratory - October 1, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <cassert>
#include <algorithm>
#include <iomanip>

#include <TLorentzVector.h>
#include <TVector3.h>
#include <TSystem.h>

#include "Conventions/GBuild.h"
#include "Conventions/Units.h"
#include "GHEP/GHepParticle.h"
#include "GHEP/GHepRecord.h"
#include "GHEP/GHepStatus.h"
#include "GHEP/GHepFlags.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGUtils.h"
#include "PDG/PDGCodes.h"
#include "Utils/PrintUtils.h"

using std::endl;
using std::setw;
using std::setprecision;
using std::setfill;
using std::ios;

using namespace genie;

ClassImp(GHepRecord)

//___________________________________________________________________________
namespace genie {
 ostream & operator << (ostream & stream, const GHepRecord & rec)
 {
   rec.Print(stream);
   return stream;
 }
}
//___________________________________________________________________________
GHepRecord::GHepRecord() :
TClonesArray("genie::GHepParticle")
{
  this->InitRecord();
}
//___________________________________________________________________________
GHepRecord::GHepRecord(int size) :
TClonesArray("genie::GHepParticle", size)
{
  this->InitRecord();
}
//___________________________________________________________________________
GHepRecord::GHepRecord(const GHepRecord & record) :
TClonesArray("genie::GHepParticle", record.GetEntries())
{
  this->InitRecord();
  this->Copy(record);
}
//___________________________________________________________________________
GHepRecord::~GHepRecord()
{
  this->CleanRecord();
}
//___________________________________________________________________________
Interaction * GHepRecord::Summary(void) const
{
  if(!fInteraction) {
    LOG("GHEP", pWARN) << "Returning NULL interaction";
  }
  return fInteraction;
}
//___________________________________________________________________________
void GHepRecord::AttachSummary(Interaction * interaction)
{
  fInteraction = interaction;
}
//___________________________________________________________________________
GHepParticle * GHepRecord::Particle(int position) const
{
// Returns the GHepParticle from the specified position of the event record.

  if( position >=0 && position < this->GetEntries() ) {

     GHepParticle * particle = (GHepParticle *) (*this)[position];
     if(particle) return particle;
  }
  LOG("GHEP", pWARN)
        << "No GHepParticle found with: (pos = "
                              << position << ") - Returning NULL";
  return 0;
}
//___________________________________________________________________________
GHepParticle * GHepRecord::FindParticle(
                               int pdg, GHepStatus_t status, int start) const
{
// Returns the first GHepParticle with the input pdg-code and status
// starting from the specified position of the event record.

  int nentries = this->GetEntries();
  for(int i = start; i < nentries; i++) {
     GHepParticle * p = (GHepParticle *) (*this)[i];
     if(p->Status() == status && p->Pdg() == pdg) return p;
  }
  LOG("GHEP", pWARN)
        << "No GHepParticle found with: (pos >= " << start
        << ", pdg = " << pdg << ", ist = " << status << ") - Returning NULL";
  return 0;
}
//___________________________________________________________________________
int GHepRecord::ParticlePosition(
                               int pdg, GHepStatus_t status, int start) const
{
// Returns the position of the first GHepParticle with the input pdg-code
// and status starting from the specified position of the event record.

  int nentries = this->GetEntries();
  for(int i = start; i < nentries; i++) {
     GHepParticle * p = (GHepParticle *) (*this)[i];
     if(p->Status() == status && p->Pdg() == pdg) return i;
  }
  LOG("GHEP", pWARN) << "Returning invalid GHEP entry position";

  return -1;
}
//___________________________________________________________________________
int GHepRecord::ParticlePosition(GHepParticle * particle, int start) const
{
// Returns the position of the first match with the specified GHepParticle
// starting from the specified position of the event record.

  int nentries = this->GetEntries();
  for(int i = start; i < nentries; i++) {
     GHepParticle * p = (GHepParticle *) (*this)[i];
     if( p->Compare(particle) ) return i;
  }
  LOG("GHEP", pWARN) << "Returning invalid GHEP entry position";

  return -1;
}
//___________________________________________________________________________
vector<int> * GHepRecord::GetStableDescendants(int position) const
{
// Returns a list of all stable descendants of the GHEP entry in the input 
// slot. The user adopts the output vector.

  // return null if particle index is out of range
  int nentries = this->GetEntries();
  if(position<0 || position>=nentries) return 0;  

  vector<int> * descendants = new vector<int>;
  
  // return itself if it is a stable final state particle
  if(this->Particle(position)->Status() == kIStStableFinalState) {
    descendants->push_back(position);
    return descendants;
  }

  for(int i = 0; i < nentries; i++) {
    if(i==position) continue;
    GHepParticle * p = (GHepParticle *) (*this)[i];
    if(p->Status() != kIStStableFinalState) continue;
    bool is_descendant=false;
    int mom = p->FirstMother();
    while(mom>-1) {
      if(mom==position) is_descendant=true;
      if(is_descendant) {
	descendants->push_back(i);
        break;
      }
      mom = this->Particle(mom)->FirstMother();
    }
  }
  return descendants;
}
//___________________________________________________________________________
GHepParticle * GHepRecord::Probe(void) const
{
// Returns the GHepParticle representing the probe (neutrino, e,...).

  int ipos = this->ProbePosition();
  if(ipos>-1) return this->Particle(ipos);
  return 0;
}
//___________________________________________________________________________
GHepParticle * GHepRecord::TargetNucleus(void) const
{
// Returns the GHepParticle representing the target nucleus, or 0 if it does
// not exist.

  int ipos = this->TargetNucleusPosition();
  if(ipos>-1) return this->Particle(ipos);
  return 0;
}
//___________________________________________________________________________
GHepParticle * GHepRecord::RemnantNucleus(void) const
{
// Returns the GHepParticle representing the remnant nucleus, or 0 if it does
// not exist.

  int ipos = this->RemnantNucleusPosition();
  if(ipos>-1) return this->Particle(ipos);
  return 0;
}
//___________________________________________________________________________
GHepParticle * GHepRecord::HitNucleon(void) const
{
// Returns the GHepParticle representing the struck nucleon, or 0 if it does
// not exist.

  int ipos = this->HitNucleonPosition();
  if(ipos>-1) return this->Particle(ipos);
  return 0;
}
//___________________________________________________________________________
GHepParticle * GHepRecord::HitElectron(void) const
{
// Returns the GHepParticle representing the struck electron, or 0 if it does
// not exist.

  int ipos = this->HitElectronPosition();
  if(ipos>-1) return this->Particle(ipos);
  return 0;
}
//___________________________________________________________________________
GHepParticle * GHepRecord::FinalStatePrimaryLepton(void) const
{
// Returns the GHepParticle representing the final state primary lepton.

  int ipos = this->FinalStatePrimaryLeptonPosition();
  if(ipos>-1) return this->Particle(ipos);
  return 0;
}
//___________________________________________________________________________
GHepParticle * GHepRecord::FinalStateHadronicSystem(void) const
{
// Returns the GHepParticle representing the sum of the DIS pre-fragm f/s
// hadronic system, or 0 if it does not exist.

  int ipos = this->FinalStateHadronicSystemPosition();
  if(ipos>-1) return this->Particle(ipos);
  return 0;
}
//___________________________________________________________________________ 
int GHepRecord::ProbePosition(void) const
{
// Returns the GHEP position of the GHepParticle representing the probe 
// (neutrino, e,...).

  return 0; // The probe is *always* at slot 0.
}
//___________________________________________________________________________
int GHepRecord::TargetNucleusPosition(void) const
{
// Returns the GHEP position of the GHepParticle representing the target 
// nucleus - or -1 if the interaction takes place at a free nucleon.

  GHepParticle * p = this->Particle(1); // If exists, it will be at slot 1
  if(!p) return -1;

  if(p->IsNucleus() && p->Status()==kIStInitialState) return 1; 

  return -1;
}
//___________________________________________________________________________
int GHepRecord::RemnantNucleusPosition(void) const
{
// Returns the GHEP position of the GHepParticle representing the remnant
// nucleus - or -1 if the interaction takes place at a free nucleon.

  GHepParticle * p = this->TargetNucleus();
  if(!p) return -1;

  int dau1 = p->FirstDaughter();
  int dau2 = p->LastDaughter();

  if(dau1==-1 && dau2==-1) return -1;

  for(int i=dau1; i<=dau2; i++) {
    GHepParticle * dp = this->Particle(i);
    if(dp->IsNucleus() && dp->Status()==kIStStableFinalState) return i; 
  }
  return -1;
}
//___________________________________________________________________________
int GHepRecord::HitNucleonPosition(void) const
{
// Returns the GHEP position of the GHepParticle representing the hit nucleon.
// If a struck nucleon is set it will be at slot 2 (for scattering off nuclear
// targets) or at slot 1 (for free nucleon scattering).
// If the struck nucleon is not set (eg coherent scattering, ve- scattering) 
// it returns 0.

  GHepParticle * nucleus = this->TargetNucleus();

  int          ipos = (nucleus) ? 2 : 1;
  GHepStatus_t ist  = (nucleus) ? kIStNucleonTarget : kIStInitialState;

  GHepParticle * p = this->Particle(ipos);
  if(!p) return -1;

  bool isN = pdg::IsNeutronOrProton(p->Pdg());
  if(isN && p->Status()==ist) return ipos; 

  return -1;
}
//___________________________________________________________________________
int GHepRecord::HitElectronPosition(void) const
{
// Returns the GHEP position of the GHepParticle representing a hit electron.
// Same as above..

  GHepParticle * nucleus = this->TargetNucleus();

  int ipos = (nucleus) ? 2 : 1;

  GHepParticle * p = this->Particle(ipos);
  if(!p) return -1;

  bool ise = pdg::IsElectron(p->Pdg());
  if(ise && p->Status()==kIStInitialState) return ipos; 

  return -1;
}
//___________________________________________________________________________
int GHepRecord::FinalStatePrimaryLeptonPosition(void) const
{
// Returns the GHEP position GHepParticle representing the final state 
// primary lepton.

  GHepParticle * probe = this->Probe();
  if(!probe) return -1;

  int ifsl = probe->FirstDaughter();
  return ifsl;
}
//___________________________________________________________________________
int GHepRecord::FinalStateHadronicSystemPosition(void) const
{
  return this->ParticlePosition(
                        kPdgHadronicSyst,kIStDISPreFragmHadronicState,0);
}
//___________________________________________________________________________ 
unsigned int GHepRecord::NEntries(int pdg, GHepStatus_t ist, int start) const
{
  unsigned int nentries = 0;

  for(int i = start; i < this->GetEntries(); i++) {
     GHepParticle * p = (GHepParticle *) (*this)[i];
     if(p->Pdg()==pdg && p->Status()==ist) nentries++;
  }
  return nentries;
}
//___________________________________________________________________________
unsigned int GHepRecord::NEntries(int pdg, int start) const
{
  unsigned int nentries = 0;

  for(int i = start; i < this->GetEntries(); i++) {
     GHepParticle * p = (GHepParticle *) (*this)[i];
     if(p->Pdg()==pdg) nentries++;
  }
  return nentries;
}
//___________________________________________________________________________
void GHepRecord::AddParticle(const GHepParticle & p)
{
// Provides a simplified method for inserting entries in the TClonesArray

  unsigned int pos = this->GetEntries();

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("GHEP", pINFO)
    << "Adding particle with pdgc = " << p.Pdg() << " at slot = " << pos;
#endif
  new ((*this)[pos]) GHepParticle(p);

  // Update the mother's daughter list. If the newly inserted particle broke
  // compactification, then run CompactifyDaughterLists()
  this->UpdateDaughterLists();
}
//___________________________________________________________________________
void GHepRecord::AddParticle(
  int pdg, GHepStatus_t status, int mom1, int mom2, int dau1, int dau2,
                         const TLorentzVector & p, const TLorentzVector & v)
{
// Provides a 'simplified' method for inserting entries in the TClonesArray

  unsigned int pos = this->GetEntries();

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("GHEP", pINFO)
           << "Adding particle with pdgc = " << pdg << " at slot = " << pos;
#endif
  new ((*this)[pos]) GHepParticle(pdg,status, mom1,mom2,dau1,dau2, p, v);

  // Update the mother's daughter list. If the newly inserted particle broke
  // compactification, then run CompactifyDaughterLists()
  this->UpdateDaughterLists();
}
//___________________________________________________________________________
void GHepRecord::AddParticle(
  int pdg, GHepStatus_t status, int mom1, int mom2, int dau1, int dau2,
                               double px, double py, double pz, double E,
                                     double x, double y, double z, double t)
{
// Provides a 'simplified' method for inserting entries in the TClonesArray

  unsigned int pos = this->GetEntries();

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("GHEP", pINFO)
           << "Adding particle with pdgc = " << pdg << " at slot = " << pos;
#endif
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

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("GHEP", pINFO)
     << "Updating the daughter-list for the mother of particle at: " << pos;
#endif

  GHepParticle * p = this->Particle(pos);
  assert(p);

  int mom_pos = p->FirstMother();
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("GHEP", pINFO) << "Mother particle is at slot: " << mom_pos;
#endif
  if(mom_pos==-1) return; // may not have mom (eg init state)
  GHepParticle * mom = this->Particle(mom_pos);
  if(!mom) return; // may not have mom (eg init state)

  int dau1 = mom->FirstDaughter();
  int dau2 = mom->LastDaughter();

  // handles the case where the daughter list was initially empty
  if(dau1 == -1) {
     mom->SetFirstDaughter(pos);
     mom->SetLastDaughter(pos);
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
     LOG("GHEP", pINFO)
        << "Done! Daughter-list is compact: [" << pos << ", " << pos << "]";
#endif
     return;
  }
  // handles the case where the new daughter is added at the slot just before
  // an already compact daughter list
  if(pos == dau1-1) {
     mom->SetFirstDaughter(pos);
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
     LOG("GHEP", pINFO)
       << "Done! Daughter-list is compact: [" << pos << ", " << dau2 << "]";
#endif
     return;
  }
  // handles the case where the new daughter is added at the slot just after
  // an already compact daughter list
  if(pos == dau2+1) {
     mom->SetLastDaughter(pos);
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
     LOG("GHEP", pINFO)
       << "Done! Daughter-list is compact: [" << dau1 << ", " << pos << "]";
#endif
     return;
  }

  // If you are here, then the last particle insertion broke the daughter
  // list compactification - Run the compactifier
  LOG("GHEP", pNOTICE)
      << "Daughter-list is not compact - Running compactifier";
  this->CompactifyDaughterLists();
}
//___________________________________________________________________________
void GHepRecord::RemoveIntermediateParticles(void)
{
  LOG("GHEP", pNOTICE) << "Removing all intermediate particles from GHEP";
  this->Compress(); 

  int i=0;
  GHepParticle * p = 0;

  TIter iter(this);
  while( (p = (GHepParticle *)iter.Next()) ) {

    if(!p) continue;
    GHepStatus_t ist = p->Status();

    bool keep = (ist==kIStInitialState) || 
                (ist==kIStStableFinalState) || (ist==kIStNucleonTarget);
    if(keep) {
       p->SetFirstDaughter(-1);
       p->SetLastDaughter(-1);
       p->SetFirstMother(-1);
       p->SetLastMother(-1);
    } else {
       LOG("GHEP", pNOTICE) 
                   << "Removing: " << p->Name() << " from slot: " << i;
       this->RemoveAt(i);
    }
    i++;
  }
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("GHEP", pDEBUG) << "Compressing GHEP record to remove empty slots";
#endif
  this->Compress(); 
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
            GHepParticle * p = this->Particle(k);
            if(p->FirstMother() == i) {
                 ndau++;
                 this->SwapParticles(start,k);
                 start++;
            }
        }
        if(ndau>0) {
          this->Particle(i)->SetFirstDaughter(dau1);
          this->Particle(i)->SetLastDaughter(dau1+ndau);
        } else {
          this->Particle(i)->SetFirstDaughter(-1);
          this->Particle(i)->SetLastDaughter(-1);
        }
     } //!compact
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
     LOG("GHEP", pINFO)
        << "Done ompactifying daughter-list for particle at slot: " << i;
#endif
  }
  this->FinalizeDaughterLists();
}
//___________________________________________________________________________
bool GHepRecord::HasCompactDaughterList(int pos)
{
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("GHEP", pDEBUG) << "Examining daughter-list of particle at: " << pos;
#endif
  vector<int> daughters;
  GHepParticle * p = 0;
  TIter iter(this);
  int i=0;
  while( (p = (GHepParticle *)iter.Next()) ) {
    if(p->FirstMother() == pos) {
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
       LOG("GHEP", pDEBUG) << "Particle at: " << i << " is a daughter";
       daughters.push_back(i);
#endif
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
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("GHEP", pINFO)
      << "Daughter-list of particle at: " << pos << " is "
                            << (is_compact ? "" : "not") << " compact";
#endif
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
    if(ist != kIStInitialState && ist != kIStNucleonTarget) return pos;
    pos++;
  }
  return pos;
}
//___________________________________________________________________________
void GHepRecord::SwapParticles(int i, int j)
{
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("GHEP", pINFO) << "Swapping GHepParticles : " << i << " <--> " << j;
#endif
  int n = this->GetEntries();
  assert(i>=0 && j>=0 && i<n && j<n);

  if(i==j) return;

  GHepParticle * pi  = this->Particle(i);
  GHepParticle * pj  = this->Particle(j);
  GHepParticle * tmp = new GHepParticle(*pi);

  pi->Copy(*pj);
  pj->Copy(*tmp);

  delete tmp;

  // tell their daughters
  if(pi->HasDaughters()) {
    for(int k=pi->FirstDaughter(); k<=pi->LastDaughter(); k++)
      this->Particle(k)->SetFirstMother(j);
  }
  if(pj->HasDaughters()) {
    for(int k=pj->FirstDaughter(); k<=pj->LastDaughter(); k++)
      this->Particle(k)->SetFirstMother(i);
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
void GHepRecord::SetVertex(double x, double y, double z, double t)
{
  fVtx->SetXYZT(x,y,z,t);
}
//___________________________________________________________________________
void GHepRecord::SetVertex(const TLorentzVector & vtx)
{
  fVtx->SetXYZT(vtx.X(),vtx.Y(),vtx.Z(),vtx.T());
}
//___________________________________________________________________________
void GHepRecord::InitRecord(void)
{
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("GHEP", pDEBUG) << "Initializing GHepRecord";
#endif
  fInteraction = 0;
  fWeight      = 1.;
  fProb        = 1.;
  fXSec        = 0.;
  fDiffXSec    = 0.;
  fVtx         = new TLorentzVector(0,0,0,0);
  fEventFlags  = new TBits(GHepFlags::NFlags());

  this->SetOwner(true);
}
//___________________________________________________________________________
void GHepRecord::CleanRecord(void)
{
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("GHEP", pDEBUG) << "Cleaning up GHepRecord";
#endif
  this->Clear("C");
}
//___________________________________________________________________________
void GHepRecord::ResetRecord(void)
{
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("GHEP", pDEBUG) << "Reseting GHepRecord";
#endif
  this->CleanRecord();
  this->InitRecord();
}
//___________________________________________________________________________
void GHepRecord::Clear(Option_t * opt)
{
  if (fInteraction) delete fInteraction;
  fInteraction=0;

  if (fVtx) delete fVtx;
  fVtx=0;

  if(fEventFlags) delete fEventFlags;
  fEventFlags=0;

  TClonesArray::Clear(opt);

//  if (fInteraction) delete fInteraction;
//  delete fVtx;
//
//  delete fEventFlags;
//
//  TClonesArray::Clear(opt);
}
//___________________________________________________________________________
void GHepRecord::Copy(const GHepRecord & record)
{
  // clean up
  this->ResetRecord();

  // copy event record entries
  unsigned int ientry = 0;
  GHepParticle * p = 0;
  TIter ghepiter(&record);
  while ( (p = (GHepParticle *) ghepiter.Next()) )
                              new ( (*this)[ientry++] ) GHepParticle(*p);

  // copy summary
  fInteraction = new Interaction( *record.fInteraction );

  // copy flags
  *fEventFlags = *(record.EventFlags());

  // copy vtx position
  TLorentzVector * v = record.Vertex();
  fVtx->SetXYZT(v->X(),v->Y(),v->Z(),v->T());

  // copy weights & xsecs
  fWeight   = record.fWeight;
  fProb     = record.fProb;
  fXSec     = record.fXSec;
  fDiffXSec = record.fDiffXSec;
}
//___________________________________________________________________________
void GHepRecord::Print(ostream & stream) const
{
  // Check $GHEPPRINTLEVEL for the preferred GHEP printout detail level
  //  0 -> prints particle list
  //  1 -> prints particle list + event flags
  //  2 -> prints particle list + event flags + wght/xsec
  //  3 -> prints particle list + event flags + wght/xsec + summary
  // 10 -> as in level 0 but showing particle positions too
  // 11 -> as in level 1 but showing particle positions too
  // 12 -> as in level 2 but showing particle positions too
  // 13 -> as in level 3 but showing particle positions too

  int  printlevel = 1;
  bool showpos    = false; 
  if( gSystem->Getenv("GHEPPRINTLEVEL") ) {
     printlevel = atoi( gSystem->Getenv("GHEPPRINTLEVEL") );

     bool accept = (printlevel>= 0 && printlevel<= 3) ||
                   (printlevel>=10 && printlevel<=13);
     if(accept) {
       if(printlevel>=10) {
          printlevel-=10;
          showpos=true;
       }
     }
  }

  // start printing the record

  stream << "\n\n|";
  stream << setfill('-') << setw(115) << "|";

  stream << "\n|GENIE GHEP Event Record [shown using $GHEPPRINTLEVEL = "
         << printlevel << "]" << setfill(' ') << setw(58) << "|";

  stream << "\n|";
  stream << setfill('-') << setw(115) << "|";

  stream << "\n| ";
  stream << setfill(' ') << setw(6)  << "Idx | "
         << setfill(' ') << setw(16) << "Name | "
         << setfill(' ') << setw(6)  << "Ist | "
         << setfill(' ') << setw(13) << "PDG | "
         << setfill(' ') << setw(12) << "Mother  | "
         << setfill(' ') << setw(12) << "Daughter  | "
         << setfill(' ') << setw(10) << ((showpos) ? "Px(x) |" : "Px | ")
         << setfill(' ') << setw(10) << ((showpos) ? "Py(y) |" : "Py | ")
         << setfill(' ') << setw(10) << ((showpos) ? "Pz(z) |" : "Pz | ")
         << setfill(' ') << setw(10) << ((showpos) ?  "E(t) |" :  "E | ")
         << setfill(' ') << setw(10) << "m  | ";

  stream << "\n|";
  stream << setfill('-') << setw(115) << "|";

  GHepParticle * p = 0;
  TObjArrayIter piter(this);
  TVector3 polarization(0,0,0);

  unsigned int idx = 0;

  double sum_E  = 0;
  double sum_px = 0;
  double sum_py = 0;
  double sum_pz = 0;

  while( (p = (GHepParticle *) piter.Next()) ) {

     stream << "\n| ";
     stream << setfill(' ') << setw(3)  << idx++               << " | ";
     stream << setfill(' ') << setw(13) << p->Name()           << " | ";
     stream << setfill(' ') << setw(3)  << p->Status()         << " | ";
     stream << setfill(' ') << setw(10) << p->Pdg()            << " | ";
     stream << setfill(' ') << setw(3)  << p->FirstMother()    << " | ";
     stream << setfill(' ') << setw(3)  << p->LastMother()     << " | ";
     stream << setfill(' ') << setw(3)  << p->FirstDaughter()  << " | ";
     stream << setfill(' ') << setw(3)  << p->LastDaughter()   << " | ";
     stream << setiosflags(ios::fixed) << setprecision(3);
     stream << setfill(' ') << setw(7)  << p->Px()             << " | ";
     stream << setfill(' ') << setw(7)  << p->Py()             << " | ";
     stream << setfill(' ') << setw(7)  << p->Pz()             << " | ";
     stream << setfill(' ') << setw(7)  << p->E()              << " | ";

     if (p->IsOnMassShell())
        stream << setfill(' ') << setw(7)  << p->Mass()        << " | ";
     else
        stream << setfill('*') << setw(7)  << p->Mass()        << " | M = " 
               << p->GetP4()->M() << " ";

     if (p->PolzIsSet()) {
       p->GetPolarization(polarization);
       stream << "P = (" << polarization.x() << "," << polarization.y()
              << "," << polarization.z() << ")";
     }

     // plot particle position if requested
     if(showpos) {
       stream << "\n| ";
       stream << setfill(' ') << setw(6)  << " | ";
       stream << setfill(' ') << setw(16) << " | ";
       stream << setfill(' ') << setw(6)  << " | ";
       stream << setfill(' ') << setw(13) << " | ";
       stream << setfill(' ') << setw(6)  << " | ";
       stream << setfill(' ') << setw(6)  << " | ";
       stream << setfill(' ') << setw(6)  << " | ";
       stream << setfill(' ') << setw(6)  << " | ";
       stream << setiosflags(ios::fixed)  << setprecision(3);
       stream << setfill(' ') << setw(7)  << p->Vx()  << " | ";
       stream << setfill(' ') << setw(7)  << p->Vy()  << " | ";
       stream << setfill(' ') << setw(7)  << p->Vz()  << " | ";
       stream << setfill(' ') << setw(7)  << p->Vt()  << " | ";
       stream << setfill(' ') << setw(10) << " | ";
     }

     // compute P4Final - P4Initial
     //
     // Take into account real particles and fake (generator-specific)
     // particles (rootino, bindino, ...) used to record non-fake physics.
     // Ignore initial & final state ions (if any).

     if(p->Status() == kIStStableFinalState) {

          sum_E  += p->E();
          sum_px += p->Px();
          sum_py += p->Py();
          sum_pz += p->Pz();
     } else 
     if(p->Status() == kIStInitialState) {
       /*
     if(p->Status() == kIStInitialState || p->Status() == kIStNucleonTarget) {
       */
          sum_E  -= p->E();
          sum_px -= p->Px();
          sum_py -= p->Py();
          sum_pz -= p->Pz();
     }

  } // loop over particles

  stream << "\n|";
  stream << setfill('-') << setw(115) << "|";

  // Print SUMS
  stream << "\n| ";
  stream << setfill(' ') << setw(17) << "Fin-Init:| "
         << setfill(' ') << setw(6)  << "    | "
         << setfill(' ') << setw(18) << "    | "
         << setfill(' ') << setw(12) << "        | "
         << setfill(' ') << setw(12) << "          | ";
  stream << setiosflags(ios::fixed)  << setprecision(3);
  stream << setfill(' ') << setw(7)  << sum_px  << " | ";
  stream << setfill(' ') << setw(7)  << sum_py  << " | ";
  stream << setfill(' ') << setw(7)  << sum_pz  << " | ";
  stream << setfill(' ') << setw(7)  << sum_E   << " | ";
  stream << setfill(' ') << setw(10)  << "   | ";

  stream << "\n|";
  stream << setfill('-') << setw(115) << "|";

  // Print vertex

  GHepParticle * probe = this->Probe();
  if(probe){
    stream << "\n| ";
    stream << setfill(' ') << setw(17) << "Vertex:  | ";
    stream << setfill(' ') << setw(11)
                       << ((probe) ? probe->Name() : "unknown probe") << " @ (";

    stream << setiosflags(ios::fixed)  << setprecision(5);
    stream << "x = " << setfill(' ') << setw(11) << this->Vertex()->X() << " m, ";
    stream << "y = " << setfill(' ') << setw(11) << this->Vertex()->Y() << " m, ";
    stream << "z = " << setfill(' ') << setw(11) << this->Vertex()->Z() << " m, ";
    stream << setiosflags(ios::scientific) << setprecision(6);
    stream << "t = " << setfill(' ') << setw(15) << this->Vertex()->T() << " s) ";
    stream << setiosflags(ios::fixed)  << setprecision(3);
    stream << setfill(' ') << setw(2)  << "|";

    stream << "\n|";
    stream << setfill('-') << setw(115) << "|";
  }

  // Print FLAGS

  if(printlevel>=1) {
    stream << "\n| ";
    stream << setfill(' ') << setw(17) << "FLAGS:   | "
           << "UnPhys: " << setfill(' ') << setw(5) 
           << utils::print::BoolAsIOString(this->IsUnphysical()) << " |"
           << " ErrBits[" << fEventFlags->GetNbits() << "->0]:" 
           << *fEventFlags << " |" 
           << " 1stSet: " << setfill(' ') << setw(38) 
           << ( this->IsUnphysical() ? 
                 GHepFlags::Describe(GHepFlag_t(fEventFlags->FirstSetBit())) : 
                 "none") << "| ";
    stream << "\n|";
    stream << setfill('-') << setw(115) << "|";
  }

  if(printlevel>=2) {
    stream << "\n| ";
    stream << setiosflags(ios::scientific) << setprecision(5);

    stream << setfill(' ') << setw(17) << "XSC/WGT: | "
           << setfill(' ') << setw(17) << "XSec[Event] = "
           << fXSec/units::cm2 << " cm^2  |"
           << setfill(' ') << setw(17) << " XSec[Kinematics] = "
           << fDiffXSec/units::cm2 << " cm^2/{K}  |"
           << setfill(' ') << setw(17) << " Weight = "
           << setiosflags(ios::fixed) << setprecision(5)
           << fWeight   << " |";

    stream << "\n|";
    stream << setfill('-') << setw(115) << "|";
  }

  stream << "\n";
  stream << setfill(' ');

  if(printlevel==3) {
     if(fInteraction) stream << *fInteraction;
     else             stream << "NULL Interaction!" << endl;
  }
  stream << "\n";
}
//___________________________________________________________________________
