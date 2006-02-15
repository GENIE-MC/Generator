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
     LOG("Kinematics", pINFO) << "max{dxsec/dK} = " << xsec_max;
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

  // get neutrino energy
  double E = this->Energy(interaction);
  LOG("Kinematics", pINFO) << "E = " << E;

  // access the the cache branch
  TNtuple * cache_branch = this->AccessCacheBranch(interaction);

  // build the search rule
  double dE = TMath::Min(0.25, 0.05*E);
  ostringstream search;
  search << "(E-" << E << " < " << dE << ") && (E>=" << E << ")";

  // query for all the entries at a window around the current energy
  TSQLResult * result = cache_branch->Query("E:xsec", search.str().c_str());
  int nrows = result->GetRowCount();
  LOG("Kinematics", pDEBUG)
            << "Found " << nrows << " rows with " << search.str();
  if(nrows <= 0) {
     delete result;
     return -1;
  }

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
     delete row;
  }
  delete result;

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

  TNtuple * cache_branch = this->AccessCacheBranch(interaction);

  double E = this->Energy(interaction);
  if(max_xsec>0) cache_branch->Fill(E, max_xsec);
}
//___________________________________________________________________________
double KineGeneratorWithCache::Energy(const Interaction * interaction) const
{
// Returns the neutrino energy at the struck nucleon rest frame. Kinematic
// generators should override this method if they need to cache the max-xsec
// values for another energy value (eg kinematic generators for IMD or COH)

  const InitialState & init_state = interaction->GetInitialState();
  double E = init_state.GetProbeE(kRfStruckNucAtRest);
  return E;
}
//___________________________________________________________________________
TNtuple * KineGeneratorWithCache::AccessCacheBranch(
                                      const Interaction * interaction) const
{
// Returns the cache branch for this algorithm and this interaction. If no
// branch is found then one is created.

  Cache * cache = Cache::Instance();

  // build the cache branch key as: namespace::algorithm/config/interaction
  string algkey = this->Id().Key();
  string intkey = interaction->AsString();
  string key    = cache->CacheBranchKey(algkey, intkey);

  TNtuple * cache_branch = cache->FindCacheBranchPtr(key);
  if(!cache_branch) {
    //-- create the cache branch at the first pass
    LOG("Cache", pINFO) << "Max d^nXSec/d{K}^n cache branch doesn't exist";
    LOG("Cache", pINFO) << "Creating cache branch - key = " << key;

    string brdef = "E:xsec";
    cache_branch = cache->CreateCacheBranch(key, brdef);
  }

  assert(cache_branch);
  return cache_branch;
}
//___________________________________________________________________________
