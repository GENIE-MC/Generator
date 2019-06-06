//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab 

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Mar 05, 2010 - CA
   User options regarding whether to pass-through unphysical events take
   precedence over processing instructions read from caught exception.
 @ Feb 01, 2013 - CA
   The GUNPHYSMASK env. var is no longer used. The bit-field mask is stored
   in the GHEP record and GHepRecord::Accept() is now checked.
*/
//____________________________________________________________________________

#include <cassert>
#include <sstream>
#include <cstdlib>
#include <algorithm>

#include <TMath.h>
#include <TBits.h>
#include <TSystem.h>
#include <TStopwatch.h>

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/EventGen/EventGenerator.h"
#include "Framework/EventGen/InteractionListGeneratorI.h"
#include "Framework/EventGen/EVGThreadException.h"
#include "Framework/EventGen/GVldContext.h"
#include "Framework/GHEP/GHepVirtualListFolder.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/GHEP/GHepFlags.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Utils/PrintUtils.h"

using std::ostringstream;

using namespace genie;
using namespace genie::utils;
using namespace genie::controls;
using namespace genie::exceptions;

//___________________________________________________________________________
EventGenerator::EventGenerator() :
EventGeneratorI("genie::EventGenerator")
{
  this->Init();
}
//___________________________________________________________________________
EventGenerator::EventGenerator(string config) :
EventGeneratorI("genie::EventGenerator", config)
{
  this->Init();
}
//___________________________________________________________________________
EventGenerator::~EventGenerator()
{
  delete fWatch;
  delete fFiltUnphysMask;

  if(fEVGModuleVec) delete fEVGModuleVec;
  if(fEVGTime)      delete fEVGTime;
  if(fVldContext)   delete fVldContext;
}
//___________________________________________________________________________
void EventGenerator::ProcessEventRecord(GHepRecord * event_rec) const
{
  LOG("EventGenerator", pNOTICE) << "Generating Event...";

  //-- Clear previous virtual list folder
  LOG("EventGenerator", pNOTICE) << "Clearing the GHepVirtualListFolder";
  GHepVirtualListFolder * vlfolder = GHepVirtualListFolder::Instance();
  vlfolder->Clear();

  //-- Clean previous history + add the bootstrap record in the buffer
  fRecHistory.PurgeHistory();
  fRecHistory.AddSnapshot(-1, event_rec);

  //-- Initialize evg thread control flags
  bool ffwd = false;
  unsigned int nexceptions = 0;

  //-- Reset stop-watch
  fWatch->Reset();

  string mesgh = "Event generation thread: " + this->Id().Key() + 
                 " -> Running module: ";

  //-- Loop over the event record processing modules
  int istep=0;
  vector<const EventRecordVisitorI *>::const_iterator miter;

  for(miter = fEVGModuleVec->begin();
      miter != fEVGModuleVec->end(); ++miter)
  {
    const EventRecordVisitorI * visitor = *miter; // generation module

    string mesg = mesgh + visitor->Id().Key();
    LOG("EventGenerator", pNOTICE)
                 << utils::print::PrintFramedMesg(mesg,0,'~');
    if(ffwd) {
      LOG("EventGenerator", pNOTICE)
           << "Fast Forward flag was set - Skipping processing step!";
      continue;
    }
    try
    {
      fWatch->Start();
      visitor->ProcessEventRecord(event_rec);
      fWatch->Stop();
      fRecHistory.AddSnapshot(istep, event_rec);
      (*fEVGTime)[istep] = fWatch->CpuTime(); // sec
    }
    catch (EVGThreadException exception)
    {
      LOG("EventGenerator", pNOTICE)
           << "An exception was thrown and caught by EventGenerator!";
      LOG("EventGenerator", pNOTICE) << exception;

      nexceptions++;
      if ( nexceptions > kMaxEVGThreadExceptions ) {
         LOG("EventGenerator", pFATAL)
           << "Caught max allowed number (" << kMaxEVGThreadExceptions
           << ") of EVGThreadExceptions/thread. Aborting";
         LOG("EventGenerator", pFATAL) << "Event : \n" << *event_rec ;
         exit(1);
      }

      //
      // What should I do with this exception?
      // Check whether the user wants to get this unphysical event anyway
      // If not, follow instructions coming with the expection thrown:
      // Either step back to a particular processing step and try re-generating
      // the event or pass it through
      //

      // check whether user wants this type of unphysical events passed through
      bool accept = event_rec->Accept();
      if(accept) {
        LOG("EventGenerator", pWARN)
          << "An unphysical event was generated and was accepted by the user";
        break;
      } else {
        LOG("EventGenerator", pWARN)
          << "An unphysical event was generated and was rejected";
      }

      // now, follow instructions contained in the exception itself

      // make sure we are not asked to go at both directions...
      assert( !(exception.FastForward() && exception.StepBack()) );

      ffwd = exception.FastForward();
      if(exception.StepBack()) {

         // get return step (if return_step > current_step just ignore it)
         if(exception.ReturnStep() >= 0 && exception.ReturnStep() <= istep) {

           int rstep = exception.ReturnStep();
           LOG("EventGenerator", pNOTICE)
               << "Return at processing step " << rstep;
           advance(miter, rstep-istep-1);
           istep = rstep;

           // restore the event record as it was just before the processing
           // step we are about to return to
           LOG("EventGenerator", pNOTICE)
                  << "Restoring GHEP as it was just before the return step";
           event_rec->ResetRecord();
           istep--;
           GHepRecord * snapshot = fRecHistory[istep];
           fRecHistory.PurgeRecentHistory(istep+1);
           event_rec->Copy(*snapshot);
         } // valid-return-step
      } // step-back
    } // catch exception

    istep++;
  }

  LOG("EventGenerator", pNOTICE)
              << utils::print::PrintFramedMesg("Thread Summary",0,'*');
  LOG("EventGenerator", pNOTICE)
           << "The EventRecord was visited by all EventRecordVisitors";

  LOG("EventGenerator", pINFO) << "** Event generation timing info **";
  istep=0;
  for(miter = fEVGModuleVec->begin();
                               miter != fEVGModuleVec->end(); ++miter){
    const EventRecordVisitorI * visitor = *miter;

    BLOG("EventGenerator", pINFO)
       << "module " << visitor->Id().Key() << " -> ~"
                        << TMath::Max(0.,(*fEVGTime)[istep++]) << " s";
  }
  LOG("EventGenerator", pNOTICE) << "Done generating event!";
}
//___________________________________________________________________________
const InteractionListGeneratorI * EventGenerator::IntListGenerator(void) const
{
  return fIntListGen;
}
//___________________________________________________________________________
const XSecAlgorithmI * EventGenerator::CrossSectionAlg(void) const
{
  return fXSecModel;
}
//___________________________________________________________________________
void EventGenerator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//___________________________________________________________________________
void EventGenerator::Configure(string param_set)
{
  Algorithm::Configure(param_set);

  AddLowRegistry( AlgConfigPool::Instance() -> GlobalParameterList(), false ) ;

  this->LoadConfig();
}
//___________________________________________________________________________
const GVldContext & EventGenerator::ValidityContext(void) const
{
  return *fVldContext;
}
//___________________________________________________________________________
void EventGenerator::Init(void)
{
  fWatch        = new TStopwatch;
  fVldContext   = 0;
  fEVGModuleVec = 0;
  fEVGTime      = 0;
  fXSecModel    = 0;
  fIntListGen   = 0;

  fFiltUnphysMask = new TBits(GHepFlags::NFlags());
  fFiltUnphysMask->ResetAllBits(false);

  if (gSystem->Getenv("GUNPHYSMASK")) {
     unsigned int i=0;
     const char * bitfield = gSystem->Getenv("GUNPHYSMASK");

     while(bitfield[i]) {
        bool flag = (bitfield[i]=='1');

        if(i<GHepFlags::NFlags()) fFiltUnphysMask->SetBitNumber(i,flag);
        i++;
     }
  }
  LOG("EventGenerator", pDEBUG)
    << "Initializing unphysical event mask (" << GHepFlags::NFlags()
    << "->0) = " << *fFiltUnphysMask;
}
//___________________________________________________________________________
void EventGenerator::LoadConfig(void)
{
  if(fEVGModuleVec) delete fEVGModuleVec;
  if(fEVGTime)      delete fEVGTime;
  if(fVldContext)   delete fVldContext;

  LOG("EventGenerator", pDEBUG) << "Loading the generator validity context";

  string encoded_vld_context ;
  GetParam("VldContext", encoded_vld_context ) ;
  fVldContext = new GVldContext;
  fVldContext->Decode( encoded_vld_context );

  LOG("EventGenerator", pDEBUG) << "Loading the event generation modules";

  int nsteps ;
  GetParam("NModules", nsteps) ;
  if(nsteps == 0) {
    LOG("EventGenerator", pFATAL)
         << "EventGenerator configuration declares null visitor list!";
  }
  assert(nsteps>0);

  fEVGModuleVec = new vector<const EventRecordVisitorI *> (nsteps);
  fEVGTime      = new vector<double>(nsteps);

  for(int istep = 0; istep < nsteps; istep++) {

    ostringstream keystream;
    keystream << "Module-" << istep;
    RgKey key = keystream.str();

    RgAlg temp_alg ;
    GetParam( key, temp_alg ) ;

    SLOG("EventGenerator", pINFO)
        << " -- Loading module " << istep << " : " << temp_alg ;

    const EventRecordVisitorI * visitor =
               dynamic_cast<const EventRecordVisitorI *>(this->SubAlg(key));

    (*fEVGModuleVec)[istep] = visitor;
    (*fEVGTime)[istep]      = 0;
  }

  //-- load the interaction list generator
  RgKey ikey = "ILstGen";
  RgAlg ialg ;
  GetParam( ikey, ialg ) ;
  LOG("EventGenerator", pINFO) 
      << " -- Loading the interaction list generator: " << ialg;
  fIntListGen = 
      dynamic_cast<const InteractionListGeneratorI *> (this->SubAlg(ikey));
  assert(fIntListGen);

  //-- load the cross section model
  RgKey xkey    = "XSecModel@" + this->Id().Key();
  RgAlg xalg ;
  GetParam( xkey, xalg ) ;
  LOG("EventGenerator", pINFO) 
     << " -- Loading the cross section model: " << xalg;
  fXSecModel =
    dynamic_cast<const XSecAlgorithmI *> (
      this -> SubAlg( xkey ) ) ;
  assert(fXSecModel);
}
//___________________________________________________________________________


