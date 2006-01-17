//____________________________________________________________________________
/*!

\class   genie::KineGeneratorWithCache

\brief   Abstract class. Provides a data caching mechanism for for concrete
         implementations of the EventRecordVisitorI interface, generating
         kinematics and wishing to cache maximum differential xsecs.

         This class provides some common implementation for handling
         (retrieving, creating, searching, adding to) the cache.
         The various super-classes should implement the ComputeMaxXSec(...)
         method for computing the maximum xsec in case it has not already
         being pushed into the cache at a previous iteration.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created December 15, 2004

*/
//____________________________________________________________________________

#include <sstream>

#include <TSQLResult.h>
#include <TSQLRow.h>
#include <TMath.h>

#include "EVGCore/EVGThreadException.h"
#include "EVGModules/KineGeneratorWithCache.h"
#include "GHEP/GHepRecord.h"
#include "Messenger/Messenger.h"
#include "Utils/Cache.h"
#include "Utils/MathUtils.h"

using std::ostringstream;

using namespace genie;

//___________________________________________________________________________
KineGeneratorWithCache::KineGeneratorWithCache() :
EventRecordVisitorI()
{

}
//___________________________________________________________________________
KineGeneratorWithCache::KineGeneratorWithCache(string name) :
EventRecordVisitorI(name)
{

}
//___________________________________________________________________________
KineGeneratorWithCache::KineGeneratorWithCache(string name, string config) :
EventRecordVisitorI(name, config)
{

}
//___________________________________________________________________________
KineGeneratorWithCache::~KineGeneratorWithCache()
{

}
//___________________________________________________________________________
double KineGeneratorWithCache::MaxXSec(GHepRecord * event_rec) const
{
  LOG("Kinematics", pINFO)
                << "Getting max. differential xsec for the rejection method";

  double xsec_max = -1;
  Interaction * interaction = event_rec->GetInteraction();

  LOG("Kinematics", pINFO)
                  << "Attempting to find a cached max{dxsec/dK} value";
  xsec_max = this->FindMaxXSec(interaction);
  if(xsec_max>0) return xsec_max;

  LOG("Kinematics", pINFO)
                  << "Attempting to compute the max{dxsec/dK} value";
  xsec_max = this->ComputeMaxXSec(interaction);
  if(xsec_max>0) {
     this->CacheMaxXSec(interaction, xsec_max);
     return xsec_max;
  }

  LOG("Kinematics", pNOTICE)
            << "Can not generate event kinematics {K} (max_xsec({K};E)<=0)";
  // xsec for selected kinematics = 0
  event_rec->SetDiffXSec(0);
  // switch on error flag = 0
  event_rec->SwitchGenericErrFlag(true);
  // reset 'trust' bits
  interaction->ResetBit(kISkipProcessChk);
  interaction->ResetBit(kISkipKinematicChk);
  // throw exception
  genie::exceptions::EVGThreadException exception;
  exception.SetReason("kinematics generation: max_xsec({K};E)<=0");
  exception.SwitchOnFastForward();
  throw exception;

  return 0;
}
//___________________________________________________________________________
double KineGeneratorWithCache::FindMaxXSec(
                                       const Interaction * interaction) const
{
// Find a cached max xsec for the specified xsec algorithm & interaction and
// close to the specified energy

  //-- get neutrino energy
  const InitialState & init_state = interaction -> GetInitialState();
  double E = init_state.GetProbeE(kRfStruckNucAtRest);
  LOG("Kinematics", pINFO) << "E = " << E;

  //-- access the the cache branch
  TNtuple * nt = this->AccessCacheBranch(interaction);
  assert(nt);

  // build the search rule
  double dE = TMath::Min(0.25, 0.05*E);
  ostringstream search;
  search << "(E-" << E << " < " << dE << ") && (E>=" << E << ")";

  // query for all the entries at a window around the current energy
  TSQLResult * result = nt->Query("E:xsec", search.str().c_str());
  int nrows = result->GetRowCount();
  LOG("Kinematics", pDEBUG)
            << "Found " << nrows << " rows with " << search.str();
  if(nrows <= 0) return -1;

  // and now select the entry with the closest energy
  double max_xsec = -1.0;
  double Ep       = 0;
  double dEmin    = 999;

  TSQLRow * row = 0;
  while( (row = result->Next()) ) {
     double cE    = atof( row->GetField(0) );
     double cxsec = atof( row->GetField(1) );
     double dE    = TMath::Abs(E-cE);
     if(dE < dEmin) {
        max_xsec = cxsec;
        Ep       = cE;
        dEmin    = TMath::Min(dE,dEmin);
     }
  }
  LOG("Kinematics", pINFO)
     << "\nRetrieved: max xsec = " << max_xsec << " cached at E = " << Ep;

  return max_xsec;
}
//___________________________________________________________________________
void KineGeneratorWithCache::CacheMaxXSec(
                     const Interaction * interaction, double max_xsec) const
{
  LOG("Kinematics", pINFO)
                       << "Adding the computed max{dxsec/dK} value to cache";

  //-- access the the cache branch
  TNtuple * nt = this->AccessCacheBranch(interaction);
  assert(nt);

  //-- get neutrino energy
  const InitialState & init_state = interaction -> GetInitialState();
  double E = init_state.GetProbeE(kRfStruckNucAtRest);

  if(max_xsec>0) nt->Fill(E, max_xsec);
}
//___________________________________________________________________________
TNtuple * KineGeneratorWithCache::AccessCacheBranch(
                                      const Interaction * interaction) const
{
// Returns the cache branch for this algorithm & interaction. If no branch is
// found then it is created.

  Cache * cache = Cache::Instance();

  string subbranch = SelectSubBranch(interaction);
  TNtuple * nt = cache->FindCacheBranchPtr(this, subbranch);

  if(!nt) {
    //-- create the cache branch at the first pass
    LOG("Kinematics", pINFO) << "Cache branch doesn't exist / creating";
    LOG("Kinematics", pINFO) << "Branch key = "
                               << cache->CacheBranchKey(this, subbranch);
    nt = cache->CreateCacheBranch(this, subbranch, "E:xsec");
  }
  return nt;
}
//___________________________________________________________________________
string KineGeneratorWithCache::SelectSubBranch(
                                       const Interaction * interaction) const
{
// Each algorithm has a main cache branch determined by its name and the name
// of its configuration set (the same algorithm in different configurations
// has different cache branches, as it should). The name of the sub-branch is
// set as the compactified interaction string which encodes initial state,
// process and final state information).

  return interaction->AsString();
}
//___________________________________________________________________________

