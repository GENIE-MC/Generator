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

#include "EVGModules/KineGeneratorWithCache.h"
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
double KineGeneratorWithCache::MaxXSec(const Interaction * interaction) const
{
// Gets the max. cross section to be used with the rejection MC technique.
// The max. xsec cache branch for this algrothm is retrieved and searched.
// If an xsec is found then it is returned. If no xsec is found, then it is
// computed, added to the cache and returned.

  LOG("Kinematics", pINFO)
                << "Getting max. differential xsec for the rejection method";

  //-- get neutrino energy

  const InitialState & init_state = interaction -> GetInitialState();
  double E = init_state.GetProbeE(kRfStruckNucAtRest);

  LOG("Kinematics", pINFO) << "E = " << E;

  //-- create the cache branch at the first pass

  Cache * cache = Cache::Instance();

  string subbranch = SelectSubBranch(interaction);
  TNtuple * nt = cache->FindCacheBranchPtr(this, subbranch);

  if(!nt) {
    LOG("Kinematics", pINFO) << "Cache branch doesn't exist / creating";
    LOG("Kinematics", pINFO) << "Branch key = "
                               << cache->CacheBranchKey(this, subbranch);
    nt = cache->CreateCacheBranch(this, subbranch, "E:xsec");
  }

  //-- retrieve the cross section from the cache branch

  // query for all the entries at a window around the current energy
  ostringstream cut;
  //cut << "abs(E-" << E << ") < 0.5";
  cut << "(E-" << E << " < 0.5) && (E>" << E << ")";

  TSQLResult * result = nt->Query("E:xsec", cut.str().c_str());
  int nrows = result->GetRowCount();

  LOG("Kinematics", pDEBUG)
                    << "Found " << nrows << " rows with " << cut.str();

  // and now select the entry with the closest energy
  double Ep = 0; //
  double max_xsec = -1.0;
  if(nrows > 0) {
     TSQLRow * row = 0;
     double dEmin = 999;

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
     LOG("Kinematics", pDEBUG)
       << "\nRetrieved: max xsec = " << max_xsec << " cached at E = " << Ep;

     return max_xsec;
  }

  //-- if no good cached value was found, compute the max. cross section
  //   and add it to the cache ao that it can be used in a subsequent step.

  LOG("Kinematics", pINFO)
              << "No cached cross section value. Computing & caching one";

  max_xsec = ComputeMaxXSec(interaction);
  nt->Fill(E, max_xsec);

  LOG("Kinematics", pINFO) << "max xsec = " << max_xsec;

  return max_xsec;
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
