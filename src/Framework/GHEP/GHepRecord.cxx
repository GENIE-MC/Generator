//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab 

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Mar 08, 2008 - CA
   Modified the algorithm for compactifying the daughter lists. The new 
   algorithm is more robust and avoids failure modes of the old algorithm
   which appeared once simulation of nuclear de-excitations was enabled.
 @ Jun 20, 2008 - CA
   Fixed memory leak in Print()
 @ Jan 28, 2009 - CA
   When checking for energy / momentum conservation in Print(), use the new
   kIStFinalStateNuclearRemnant status code for nuclear remnants (previously
   marked as kIStStableFinalState).
 @ Sep 15, 2009 - CA
   IsNucleus() is no longer available in GHepParticle. Use pdg::IsIon().
   Print-out the particle 'rescattering code' (if it was set).
 @ May 05, 2010 - CR
   Adding special ctor for ROOT I/O purposes so as to avoid memory leak due to
   memory allocated in the default ctor when objects of this class are read by 
   the ROOT Streamer.
 @ Sep 26, 2011 - CA
   Demote a few messages from `warning' to `notice'.
 @ Nov 17, 2011 - CA
   Added `GEvGenMode_t EventGenerationMode(void) const'
 @ Nov 28, 2011 - CA
   HitNucleon() can return a hit nucleon cluster too, as needed for MEC.
 @ Jan 29, 2013 - CA
   Demote a few messages from `notice' to `info'.
 @ Jan 31, 2013 - CA
   Added static SetPrintLevel(int print_level) and corresponding static
   data member. $GHEPPRINTLEVEL env var is no longer used.
 @ Feb 01, 2013 - CA
   The GUNPHYSMASK env. var is no longer used. Added a private data member to
   store the bit-field mask in the GHEP record. Added GHepRecord::Accept() and
   GHepRecord::SetUnphysEventMask(). Tweaks in Print().
 @ Feb 06, 2013 - CA
   Added KinePhaseSpace_t fDiffXSecPhSp prov data members to specify which
   differential cross-section value is stored in fDiffXSec. Added method to 
   set it and tweaked Print() accordingly.
 @ May 02, 2013 - CA
   Added `KinePhaseSpace_t DiffXSecVars(void) const' to return fDiffXSecPhSp.

*/
//____________________________________________________________________________

#include <cassert>
#include <algorithm>
#include <iomanip>

#include <TLorentzVector.h>
#include <TVector3.h>
#include <TSystem.h>
#include <TRootIOCtor.h>

#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Units.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/GHEP/GHepStatus.h"
#include "Framework/GHEP/GHepFlags.h"
#include "Framework/GHEP/GHepUtils.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/Utils/PrintUtils.h"

using std::endl;
using std::setw;
using std::setprecision;
using std::setfill;
using std::ios;

using namespace genie;

ClassImp(GHepRecord)

int GHepRecord::fPrintLevel = 3;

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
GHepRecord::GHepRecord(TRootIOCtor*) :
TClonesArray("genie::GHepParticle"),
fInteraction(0),
fVtx(0), 
fEventFlags(0), 
fEventMask(0),
fWeight(0.),
fProb(0.),
fXSec(0.),
fDiffXSec(0.)
{

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
  LOG("GHEP", pINFO)
    << "No particle found with: (pos = " << position << ")";

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

  LOG("GHEP", pINFO)
    << "No particle found with: (pos >= " << start
    << ", pdg = " << pdg << ", ist = " << status << ")";

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

  LOG("GHEP", pINFO)
    << "No particle found with: (pos >= " << start
    << ", pdg = " << pdg << ", ist = " << status << ")";

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

  LOG("GHEP", pINFO)
    << "No particle found with pos >= " << start
    << " matching particle: " << *particle;

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
GEvGenMode_t GHepRecord::EventGenerationMode(void) const
{
  GHepParticle * p0 = this->Particle(0);
  if(!p0) return kGMdUnknown;
  GHepParticle * p1 = this->Particle(1);
  if(!p1) return kGMdUnknown;

  int p0pdg = p0->Pdg();
  GHepStatus_t p0st = p0->Status();
  int p1pdg = p1->Pdg();
  GHepStatus_t p1st = p1->Status();

  // In lepton+nucleon/nucleus mode, the 1st entry in the event record
  // is a charged or neutral lepton with status code = kIStInitialState
  if( pdg::IsLepton(p0pdg) && p0st == kIStInitialState )
  {
    return kGMdLeptonNucleus;
  }

  // In dark matter mode, the 1st entry in the event record is a dark
  // matter particle
  if( pdg::IsDarkMatter(p0pdg) && p0st == kIStInitialState )
  {
    return kGMdDarkMatterNucleus;
  }

  // In hadron+nucleon/nucleus mode, the 1st entry in the event record
  // is a hadron with status code = kIStInitialState and the 2nd entry
  // is a nucleon or nucleus with status code = kIStInitialState
  if( pdg::IsHadron(p0pdg) && p0st == kIStInitialState )
  {   
    if( (pdg::IsIon(p1pdg) || pdg::IsNucleon(p1pdg)) && p1st == kIStInitialState)
    {
       return kGMdHadronNucleus;
    }
  }

  // As above, with a photon as a probe
  if( p0pdg == kPdgGamma && p0st == kIStInitialState )
  {   
    if( (pdg::IsIon(p1pdg) || pdg::IsNucleon(p1pdg)) && p1st == kIStInitialState)
    {
       return kGMdPhotonNucleus;
    }
  }
      
  // In nucleon decay mode, 
  // - [if the decayed nucleon was a bound one] the 1st entry in the event
  //   record is a nucleus with status code = kIStInitialState and the
  //   2nd entry is a nucleon with code = kIStDecayedState
  // - [if the decayed nucleon was a free one] the first entry in the event
  //   record is a nucleon with status code = kIStInitialState and it has a
  //   single daughter which is a nucleon with status code = kIStDecayedState.
         
  if( pdg::IsIon(p0pdg)     && p0st == kIStInitialState &&
      pdg::IsNucleon(p1pdg) && p1st == kIStDecayedState)
  {
     return kGMdNucleonDecay;
  }
  if( pdg::IsNucleon(p0pdg) && p0st == kIStInitialState &&
      pdg::IsNucleon(p1pdg) && p1st == kIStDecayedState)
  {
     return kGMdNucleonDecay;
  }
         
  return kGMdUnknown;
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
// Returns the GHepParticle representing the target / initial state nucleus, 
// or 0 if it does not exist.

  int ipos = this->TargetNucleusPosition();
  if(ipos>-1) return this->Particle(ipos);
  return 0;
}
//___________________________________________________________________________
GHepParticle * GHepRecord::RemnantNucleus(void) const
{
// Returns the GHepParticle representing the remnant nucleus, 
// or 0 if it does not exist.

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

  // The probe is *always* at slot 0.
  GEvGenMode_t mode = this->EventGenerationMode();
  if(mode == kGMdLeptonNucleus || 
     mode == kGMdDarkMatterNucleus ||
     mode == kGMdHadronNucleus ||
     mode == kGMdPhotonNucleus) 
  {
    return 0;
  }
  return -1; 
}
//___________________________________________________________________________
int GHepRecord::TargetNucleusPosition(void) const
{
// Returns the GHEP position of the GHepParticle representing the target 
// nucleus - or -1 if the interaction takes place at a free nucleon.

  GEvGenMode_t mode = this->EventGenerationMode();

  if(mode == kGMdLeptonNucleus || 
     mode == kGMdDarkMatterNucleus ||
     mode == kGMdHadronNucleus ||
     mode == kGMdPhotonNucleus) 
  {
     GHepParticle * p = this->Particle(1); // If exists, it will be at slot 1
     if(!p) return -1;
     int pdgc = p->Pdg();
     if(pdg::IsIon(pdgc) && p->Status()==kIStInitialState) return 1; 
  }
  if(mode == kGMdNucleonDecay) {
     GHepParticle * p = this->Particle(0); // If exists, it will be at slot 0
     if(!p) return -1;
     int pdgc = p->Pdg();
     if(pdg::IsIon(pdgc) && p->Status()==kIStInitialState) return 0; 
  }

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
    int dpdgc = dp->Pdg();
    if(pdg::IsIon(dpdgc) && dp->Status()==kIStStableFinalState) return i; 
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

//  bool isN = pdg::IsNeutronOrProton(p->Pdg());
  bool isN = pdg::IsNucleon(p->Pdg()) || pdg::Is2NucleonCluster(p->Pdg()); 
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
  int n = this->GetEntries();
  if(n<1) return;

  int i = this->Particle(n-1)->FirstMother();
  if(i<0) return;

  //  for(int i=0; i<n; i++) {
     bool compact = this->HasCompactDaughterList(i);

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
     LOG("GHEP", pNOTICE) 
        << "Particle's " << i << " daughter list is " 
        << ((compact) ? "compact" : "__not__ compact");
#endif

     if(!compact) {
        GHepParticle * p = this->Particle(i);

        int dau1 = p->FirstDaughter();
        int dau2 = p->LastDaughter();
        int ndau = dau2-dau1+1;
        int ndp  = dau2+1;
        if(dau1==-1) {ndau=0;}

        int curr_pos = n-1;   
        while (curr_pos > ndp) {
           this->SwapParticles(curr_pos,curr_pos-1);
           curr_pos--;
        }
        if(ndau>0) {
          this->Particle(i)->SetFirstDaughter(dau1);
          this->Particle(i)->SetLastDaughter(dau2+1);
        } else {
          this->Particle(i)->SetFirstDaughter(-1);
          this->Particle(i)->SetLastDaughter(-1);
        }
        this->FinalizeDaughterLists();

     } //!compact

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
     LOG("GHEP", pINFO)
        << "Done ompactifying daughter-list for particle at slot: " << i;
#endif
     //  }
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
#endif
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
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
    LOG("GHEP", pINFO) 
      << pi->Name() << "(previously at pos: " << j 
      << ") is now at pos: " << i << " -> Notify daughters";
#endif
    for(int k=0; k<n; k++) {
      if(this->Particle(k)->FirstMother()==j) {
        this->Particle(k)->SetFirstMother(i);
      }
    }
  }

  if(pj->HasDaughters()) {
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
    LOG("GHEP", pINFO) 
      << pj->Name() << "(previously at pos: " << i 
      << ") is now at pos: " << j << " -> Notify daughters";
#endif
    for(int k=0; k<n; k++) {
      if(this->Particle(k)->FirstMother()==i) {
        this->Particle(k)->SetFirstMother(j);
      }
    }
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
  fInteraction  = 0;
  fWeight       = 1.;
  fProb         = 1.;
  fXSec         = 0.;
  fDiffXSec     = 0.;
  fDiffXSecPhSp = kPSNull;
  fVtx          = new TLorentzVector(0,0,0,0);

  fEventFlags  = new TBits(GHepFlags::NFlags());
  fEventFlags -> ResetAllBits(false);

  fEventMask   = new TBits(GHepFlags::NFlags());
//fEventMask  -> ResetAllBits(true);
  for(unsigned int i = 0; i < GHepFlags::NFlags(); i++) {
   fEventMask->SetBitNumber(i, true);
  }

  LOG("GHEP", pINFO)
    << "Initialised unphysical event mask (bits: " << GHepFlags::NFlags() - 1
    << " -> 0) : " << *fEventMask;

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

  if(fEventMask) delete fEventMask;
  fEventMask=0;

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

  // copy flags & mask
  *fEventFlags = *(record.EventFlags());
  *fEventMask  = *(record.EventMask());

  // copy vtx position
  TLorentzVector * v = record.Vertex();
  fVtx->SetXYZT(v->X(),v->Y(),v->Z(),v->T());

  // copy weights & xsecs
  fWeight       = record.fWeight;
  fProb         = record.fProb;
  fXSec         = record.fXSec;
  fDiffXSec     = record.fDiffXSec;
  fDiffXSecPhSp = record.fDiffXSecPhSp;
}
//___________________________________________________________________________
void GHepRecord::SetUnphysEventMask(const TBits & mask)
{
 *fEventMask = mask;

  LOG("GHEP", pINFO)
    << "Setting unphysical event mask (bits: " << GHepFlags::NFlags() - 1
    << " -> 0) : " << *fEventMask;
}
//___________________________________________________________________________
bool GHepRecord::Accept(void) const
{
  TBits flags = *fEventFlags;
  TBits mask  = *fEventMask;
  TBits bitwiseand = flags & mask;
  bool accept = (bitwiseand.CountBits() == 0);
  return accept;
}
//___________________________________________________________________________
void GHepRecord::SetPrintLevel(int print_level) 
{ 
  fPrintLevel = print_level; 
}
int  GHepRecord::GetPrintLevel()
{ 
  return fPrintLevel; 
}
//___________________________________________________________________________
void GHepRecord::Print(ostream & stream) const
{
  // Print levels:
  //  0 -> prints particle list
  //  1 -> prints particle list + event flags
  //  2 -> prints particle list + event flags + wght/xsec
  //  3 -> prints particle list + event flags + wght/xsec + summary
  // 10 -> as in level 0 but showing particle positions too
  // 11 -> as in level 1 but showing particle positions too
  // 12 -> as in level 2 but showing particle positions too
  // 13 -> as in level 3 but showing particle positions too

   bool accept_input_print_level = 
       (fPrintLevel >= 0 && fPrintLevel <= 3) ||
       (fPrintLevel >=10 && fPrintLevel <=13);
   
  int printlevel = (accept_input_print_level) ? fPrintLevel : 3;
  int printlevel_orig = printlevel;

  bool showpos = false; 
  if(printlevel >= 10) {
     printlevel-=10;
     showpos=true;
  }
  
  stream << "\n\n|";
  stream << setfill('-') << setw(115) << "|";

  stream << "\n|GENIE GHEP Event Record [print level: "
         << setfill(' ') << setw(3) << printlevel_orig << "]" 
         << setfill(' ') << setw(73) << "|";

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
     stream << std::fixed << setprecision(3);
     stream << setfill(' ') << setw(7)  << p->Px()             << " | ";
     stream << setfill(' ') << setw(7)  << p->Py()             << " | ";
     stream << setfill(' ') << setw(7)  << p->Pz()             << " | ";
     stream << setfill(' ') << setw(7)  << p->E()              << " | ";

     if (p->IsOnMassShell())
        stream << setfill(' ') << setw(7)  << p->Mass()        << " | ";
     else
        stream << setfill('*') << setw(7)  << p->Mass()        << " | M = " 
               << p->P4()->M() << " ";

     if (p->PolzIsSet()) {
       p->GetPolarization(polarization);
       stream << "P = (" << polarization.x() << "," << polarization.y()
              << "," << polarization.z() << ")";
     }

     if (p->RescatterCode() != -1) {
       stream << "FSI = " << p->RescatterCode();
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
       stream << std::fixed  << setprecision(3);
       stream << setfill(' ') << setw(7)  << p->Vx()  << " | ";
       stream << setfill(' ') << setw(7)  << p->Vy()  << " | ";
       stream << setfill(' ') << setw(7)  << p->Vz()  << " | ";
       stream << setfill(' ') << setw(7)  << p->Vt()  << " | ";
       stream << setfill(' ') << setw(10) << " | ";
     }

     // compute P4_{final} - P4_{nitial}

     if(p->Status() == kIStStableFinalState ||
        p->Status() == kIStFinalStateNuclearRemnant) {

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
  stream << setfill(' ') << setw(17) << "Fin-Init:  "
         << setfill(' ') << setw(6)  << "      "
         << setfill(' ') << setw(18) << "      "
         << setfill(' ') << setw(12) << "          "
         << setfill(' ') << setw(12) << "          | ";
  stream << std::fixed  << setprecision(3);
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
    stream << setfill(' ') << setw(17) << "Vertex:    ";
    stream << setfill(' ') << setw(11)
                       << ((probe) ? probe->Name() : "unknown probe") << " @ (";

    stream << std::fixed  << setprecision(5);
    stream << "x = " << setfill(' ') << setw(11) << this->Vertex()->X() << " m, ";
    stream << "y = " << setfill(' ') << setw(11) << this->Vertex()->Y() << " m, ";
    stream << "z = " << setfill(' ') << setw(11) << this->Vertex()->Z() << " m, ";
    stream << std::scientific << setprecision(6);
    stream << "t = " << setfill(' ') << setw(15) << this->Vertex()->T() << " s) ";
    stream << std::fixed  << setprecision(3);
    stream << setfill(' ') << setw(2)  << "|";

    stream << "\n|";
    stream << setfill('-') << setw(115) << "|";
  }

  // Print FLAGS

  if(printlevel>=1) {
    stream << "\n| ";
    stream << "Err flag [bits:" << fEventFlags->GetNbits()-1 << "->0] : " 
           << *fEventFlags << "    |  " 
           << "1st set: " << setfill(' ') << setw(56) 
           << ( this->IsUnphysical() ? 
                 GHepFlags::Describe(GHepFlag_t(fEventFlags->FirstSetBit())) : 
                 "none") << " | ";
    stream << "\n| ";
    stream << "Err mask [bits:" << fEventMask->GetNbits()-1 << "->0] : " 
           << *fEventMask << "    |  " 
           << "Is unphysical: " << setfill(' ') << setw(5) 
           << utils::print::BoolAsYNString(this->IsUnphysical()) << " |   "
           << "Accepted: " << setfill(' ') << setw(5) 
           << utils::print::BoolAsYNString(this->Accept()) 
           << "                          |";
    stream << "\n|";
    stream << setfill('-') << setw(115) << "|";
  }

  if(printlevel>=2) {
    stream << "\n| ";
    stream << std::scientific << setprecision(5);

    stream << "sig(Ev) = " 
           << setfill(' ') << setw(17) << fXSec/units::cm2 
           << " cm^2  |";

    switch(fDiffXSecPhSp) {
      case ( kPSyfE   ) :                                                                           
        stream << " dsig(y;E)/dy =          " << setfill(' ') << setw(13) << fDiffXSec/units::cm2 << " cm^2       |";
        break;
      case ( kPSxyfE  ) :
        stream << " d2sig(x,y;E)/dxdy =     " << setfill(' ') << setw(13) << fDiffXSec/units::cm2 << " cm^2       |";
        break;
      case ( kPSxytfE ) :
        stream << " d3sig(x,y,t;E)/dxdydt = " << setfill(' ') << setw(13) << fDiffXSec/units::cm2 << " cm^2/GeV^2 |";
        break;
      case ( kPSQ2fE  ) :
        stream << " dsig(Q2;E)/dQ2 =        " << setfill(' ') << setw(13) << fDiffXSec/units::cm2 << " cm^2/GeV^2 |";
        break;
      case ( kPSQ2vfE  ) :
        stream << " dsig(Q2,v;E)/dQ2dv =    " << setfill(' ') << setw(13) << fDiffXSec/units::cm2 << " cm^2/GeV^3 |";
        break;
      case ( kPSWQ2fE ) :
        stream << " d2sig(W,Q2;E)/dWdQ2 =   " << setfill(' ') << setw(13) << fDiffXSec/units::cm2 << " cm^2/GeV^3 |";
        break;
      default :
        stream << " dsig(Ev;{K_s})/dK   =   " << setfill(' ') << setw(13) << fDiffXSec/units::cm2 << " cm^2/{K}   |";
    }
    stream << " Weight = "
           << setfill(' ') << setw(16) 
           << std::fixed << setprecision(5)
           << fWeight   
           << " |";

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
