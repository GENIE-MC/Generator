//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab 

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Jun 25, 2008 - CA
   Using the new CacheBranchFx which avoids a significant memory leak.
 @ Feb 06, 2013 - CA
   When the value of the differential cross-section for the selected kinematics
   is set to the event, set the corresponding KinePhaseSpace_t value too.

*/
//____________________________________________________________________________

#include <sstream>
#include <cstdlib>
#include <map>

//#include <TSQLResult.h>
//#include <TSQLRow.h>
#include <TMath.h>

#include "Framework/EventGen/EVGThreadException.h"
#include "Physics/Common/KineGeneratorWithCache.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/GHEP/GHepFlags.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Utils/Cache.h"
#include "Framework/Utils/CacheBranchFx.h"
#include "Framework/Numerical/MathUtils.h"

using std::ostringstream;
using std::map;

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
  Interaction * interaction = event_rec->Summary();

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
  event_rec->SetDiffXSec(0,kPSNull);
  // switch on error flag 
  event_rec->EventFlags()->SetBitNumber(kKineGenErr, true);
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

  if(E < fEMin) {
     LOG("Kinematics", pINFO)
         << "Below minimum energy - Forcing explicit calculation";
     return -1.;
  }

  // access the the cache branch
  CacheBranchFx * cb = this->AccessCacheBranch(interaction);

  // if there are enough points stored in the cache buffer to build a
  // spline, then intepolate
  if( cb->Spl() ) {
     if( E >= cb->Spl()->XMin() && E <= cb->Spl()->XMax()) {
       double spl_max_xsec = cb->Spl()->Evaluate(E);
       LOG("Kinematics", pINFO)
          << "\nInterpolated: max xsec (E=" << E << ") = " << spl_max_xsec;
       return spl_max_xsec;      
     }
     LOG("Kinematics", pINFO)
          << "Outside spline boundaries - Forcing explicit calculation";
     return -1.;
  }

  // if there are not enough points at the cache buffer to have a spline,
  // look whether there is another point that is sufficiently close
  double dE = TMath::Min(0.25, 0.05*E);
  const map<double,double> & fmap = cb->Map();
  map<double,double>::const_iterator iter = fmap.lower_bound(E);
  if(iter != fmap.end()) {
     if(TMath::Abs(E - iter->first) < dE) return iter->second;
  }

  return -1;

/*
  // build the search rule
  double dE = TMath::Min(0.25, 0.05*E);
  ostringstream search;
  search << "(x-" << E << " < " << dE << ") && (x>=" << E << ")";

  // query for all the entries at a window around the current energy
  TSQLResult * result = cb->Ntuple()->Query("x:y", search.str().c_str());
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
  TSQLRow * row   = 0;
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
*/
}
//___________________________________________________________________________
void KineGeneratorWithCache::CacheMaxXSec(
                     const Interaction * interaction, double max_xsec) const
{
  LOG("Kinematics", pINFO)
                       << "Adding the computed max{dxsec/dK} value to cache";
  CacheBranchFx * cb = this->AccessCacheBranch(interaction);

  double E = this->Energy(interaction);
  if(max_xsec>0) cb->AddValues(E,max_xsec);

  if(! cb->Spl() ) {
    if( cb->Map().size() > 40 ) cb->CreateSpline();
  }

  if( cb->Spl() ) {
     if( E < cb->Spl()->XMin() || E > cb->Spl()->XMax() ) {
        cb->CreateSpline();
     }
  }
}
//___________________________________________________________________________
double KineGeneratorWithCache::Energy(const Interaction * interaction) const
{
// Returns the neutrino energy at the struck nucleon rest frame. Kinematic
// generators should override this method if they need to cache the max-xsec
// values for another energy value (eg kinematic generators for IMD or COH)

  const InitialState & init_state = interaction->InitState();
  double E = init_state.ProbeE(kRfHitNucRest);
  return E;
}
//___________________________________________________________________________
CacheBranchFx * KineGeneratorWithCache::AccessCacheBranch(
                                      const Interaction * interaction) const
{
// Returns the cache branch for this algorithm and this interaction. If no
// branch is found then one is created.

  Cache * cache = Cache::Instance();

  // build the cache branch key as: namespace::algorithm/config/interaction
  string algkey = this->Id().Key();
  string intkey = interaction->AsString();
  string key    = cache->CacheBranchKey(algkey, intkey);

  CacheBranchFx * cache_branch =
              dynamic_cast<CacheBranchFx *> (cache->FindCacheBranch(key));
  if(!cache_branch) {
    //-- create the cache branch at the first pass
    LOG("Kinematics", pINFO) << "No Max d^nXSec/d{K}^n cache branch found";
    LOG("Kinematics", pINFO) << "Creating cache branch - key = " << key;

    cache_branch = new CacheBranchFx("max[d^nXSec/d^n{K}] over phase space");
    cache->AddCacheBranch(key, cache_branch);
  }
  assert(cache_branch);

  return cache_branch;
}
//___________________________________________________________________________
void KineGeneratorWithCache::AssertXSecLimits(
         const Interaction * interaction, double xsec, double xsec_max) const
{
  // check the computed cross section for the current kinematics against the
  // maximum cross section used in the rejection MC method for the current
  // interaction at the current energy.
  if(xsec>xsec_max) {
    double f = 200*(xsec-xsec_max)/(xsec_max+xsec);
    if(f>fMaxXSecDiffTolerance) {
       LOG("Kinematics", pFATAL) 
    	  << "xsec: (curr) = " << xsec 
      	         << " > (max) = " << xsec_max << "\n for " << *interaction;
       LOG("Kinematics", pFATAL) 
 	  << "*** Exceeding estimated maximum differential cross section";
       exit(1);
    } else {
       LOG("Kinematics", pWARN) 
    	  << "xsec: (curr) = " << xsec 
      	         << " > (max) = " << xsec_max << "\n for " << *interaction;
       LOG("Kinematics", pWARN) 
  	    << "*** The fractional deviation of " << f << " % was allowed";
    }
  }

  // this should never happen - print an error mesg just in case...
  if(xsec<0) {
    LOG("Kinematics", pERROR) 
     << "Negative cross section for current kinematics!! \n" << *interaction;
  }
}
//___________________________________________________________________________
