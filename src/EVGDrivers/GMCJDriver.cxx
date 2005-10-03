//____________________________________________________________________________
/*!

\class   genie::GMCJDriver

\brief   GENIE MC Job Driver (event generation for the input flux & geometry)

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created May 25, 2005

*/
//____________________________________________________________________________

#include <cassert>

#include <TLorentzVector.h>
#include <TVector3.h>

#include "Algorithm/AlgConfigPool.h"
#include "Conventions/Constants.h"
#include "Conventions/Units.h"
#include "EVGDrivers/GMCJDriver.h"
#include "EVGDrivers/GEVGDriver.h"
#include "EVGDrivers/GEVGPool.h"
#include "EVGCore/EventRecord.h"
#include "FluxDrivers/GFluxI.h"
#include "Geo/GeomAnalyzerI.h"
#include "Geo/PathLengthList.h"
#include "Interaction/InitialState.h"
#include "Messenger/Messenger.h"
#include "Numerical/RandomGen.h"
#include "PDG/PDGCodeList.h"
#include "Utils/PrintUtils.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
GMCJDriver::GMCJDriver()
{
  this->Initialize();
}
//___________________________________________________________________________
GMCJDriver::~GMCJDriver()
{
  if (fGPool) delete fGPool;
}
//___________________________________________________________________________
void GMCJDriver::UseFluxDriver(GFluxI * flux_driver)
{
  fFluxDriver = flux_driver;
}
//___________________________________________________________________________
void GMCJDriver::UseGeomAnalyzer(GeomAnalyzerI * geom_analyzer)
{
  fGeomAnalyzer = geom_analyzer;
}
//___________________________________________________________________________
void GMCJDriver::CreateSplines(bool useLogE)
{
  LOG("GMCJDriver", pNOTICE)
            << "Creating missing (not already loaded) cross section splines";

  fUseSplines = true;
  fUseLogE    = useLogE;

  GEVGPool::const_iterator iter = fGPool->begin();
  for(; iter != fGPool->end(); ++iter) {
    GEVGDriver * evgdriver = iter->second;
    if(evgdriver) {
       evgdriver->CreateSplines(useLogE);
    } else {
       LOG("GMCJDriver", pERROR) << "NULL GEVGDriver";
    }
  }
}
//___________________________________________________________________________
void GMCJDriver::Initialize(void)
{
  fFluxDriver   = 0;
  fGeomAnalyzer = 0;
  fGPool        = 0;
  fPmax         = 0;
  fUseSplines   = false;

  assert( Messenger::Instance()     ); // send prolific printout on top
  assert( AlgConfigPool::Instance() ); //           --//--
}
//___________________________________________________________________________
void GMCJDriver::Configure(void)
{
  LOG("GMCJDriver", pINFO)
                  << print_utils::PrintFramedMesg("Configuring GMCJDriver");

  //-- Ask the input GFluxI for a list of all neutrino-types it will be using
  LOG("GMCJDriver", pINFO)
              << "Asking input GFluxI for the list of flux particles";
  const PDGCodeList & nulist = fFluxDriver->FluxParticles();
  LOG("GMCJDriver", pINFO) << "Flux particles: " << nulist;

  //-- Ask the ROOTGeomHandler for a list of all target materials
  LOG("GMCJDriver", pINFO)
     << "Asking input GeomAnalyzerI for the list of target materials";
  const PDGCodeList & tgtlist = fGeomAnalyzer->ListOfTargetNuclei();
  LOG("GMCJDriver", pINFO) << "Target materials: " << tgtlist;

  //-- Ask the ROOTGeomHandler for a the max. path length for each material
  //   in the list of all target materials (to be used in calculating the max.
  //   interaction probability, Pmax,  to scale all interaction probabilities)
  LOG("GMCJDriver", pINFO)
             << "Asking input GeomAnalyzerI for the max path-lengths";
  const PathLengthList & plmax = fGeomAnalyzer->ComputeMaxPathLengths();
  LOG("GMCJDriver", pINFO) << "Maximum path lengths: " << plmax;

  //-- Ask the input GFluxI for the max. neutrino energy (to compute Pmax)
  LOG("GMCJDriver", pINFO)
            << "Asking input GFluxI for maximum flux neutrino energy";
  double Emax = fFluxDriver->MaxEnergy();
  LOG("GMCJDriver", pINFO) << "Maximum flux neutrino energy = " << Emax;

  //-- Create all possible initial states and initialize/store a
  //   GEVGDriver driver object for each one of them.

  LOG("GMCJDriver", pINFO)
        << "Creating GEVGPool & adding a GEVGDriver object per init-state";

  if (fGPool) delete fGPool;
  fGPool = new GEVGPool;

  PDGCodeList::const_iterator nuiter;
  PDGCodeList::const_iterator tgtiter;

  for(nuiter = nulist.begin(); nuiter != nulist.end(); ++nuiter) {
   for(tgtiter = tgtlist.begin(); tgtiter != tgtlist.end(); ++tgtiter) {

     int target_pdgc   = *tgtiter;
     int neutrino_pdgc = *nuiter;

     InitialState init_state(target_pdgc, neutrino_pdgc);

     LOG("GMCJDriver", pINFO)
       << "\n\n ---- Creating a GEVGDriver object configured for init-state: "
       << init_state.AsString() << " ----\n\n";

     GEVGDriver * evgdriver = new GEVGDriver;
     evgdriver->SetInitialState(init_state);

     LOG("GMCJDriver", pDEBUG) << "Adding new GEVGDriver object to GEVGPool";
     fGPool->insert( GEVGPool::value_type(init_state.AsString(), evgdriver) );
   } // targets
  } // neutrinos

  LOG("GMCJDriver", pINFO)
             << "All necessary GEVGDriver object were pushed into GEVGPool\n";


  // If requested splines then coordinate spline creation from all GEVGDriver
  // objects pushed into GEVGPool. This will create all xsec splines needed
  // for all simulated processes involving the particles in the input flux
  // and geometry. Spline creation will be skipped for every spline that has
  // been pre-loaded into the the XSecSplineList

  if(fUseSplines) {
     LOG("GMCJDriver", pINFO) << "Creating cross section splines";

     for(nuiter = nulist.begin(); nuiter != nulist.end(); ++nuiter) {
       for(tgtiter = tgtlist.begin(); tgtiter != tgtlist.end(); ++tgtiter) {

         int target_pdgc   = *tgtiter;
         int neutrino_pdgc = *nuiter;

         InitialState init_state(target_pdgc, neutrino_pdgc);

         LOG("GMCJDriver", pINFO)
             << "Computing all splines needed for init-state: "
                                               << init_state.AsString();
         GEVGDriver * evgdriver = fGPool->FindDriver(init_state);
         evgdriver->CreateSplines(fUseLogE);
       } // targets
     } // neutrinos
     LOG("GMCJDriver", pINFO) << "Finished creating cross section splines\n";
  }

  // Compute the max. interaction probability to scale all interaction
  // probabilities to be computed by this driver

  LOG("GMCJDriver", pINFO)
    << "Computing the max. interaction probability (probability scale)";

  fPmax = 0; // maximum interaction probability
  TLorentzVector nup4(0,0,Emax,Emax);

  for(nuiter = nulist.begin(); nuiter != nulist.end(); ++nuiter) {
   for(tgtiter = tgtlist.begin(); tgtiter != tgtlist.end(); ++tgtiter) {

     int target_pdgc   = *tgtiter;
     int neutrino_pdgc = *nuiter;

     InitialState init_state(target_pdgc, neutrino_pdgc);

     LOG("GMCJDriver", pINFO)
           << "Computing Pmax for init-state: " << init_state.AsString();

     GEVGDriver * evgdriver = fGPool->FindDriver(init_state); // get the appropriate driver

     double xsec_sum    = evgdriver->SumCrossSection(nup4);  // interaction xsec
     double path_length = plmax.PathLength(target_pdgc); // max L for material

     double Pmax = this->PInt(xsec_sum, path_length);
     fPmax += Pmax;

     LOG("GMCJDriver", pINFO)
              << "Pmax[" << init_state.AsString()
                               << ", for given flux & geometry] = " << Pmax;
   } // targets
  } // neutrinos

  LOG("GMCJDriver", pINFO) << "Total Pmax[interaction] = " << fPmax << "\n";

  LOG("GMCJDriver", pINFO) << "Finished configuring GMCJDriver\n\n";
}
//___________________________________________________________________________
EventRecord * GMCJDriver::GenerateEvent(void)
{
  LOG("GMCJDriver", pINFO) << "Generating next event...";

  //-- generate a neutrino using the input GFluxI & get current pdgc/p4/x4

  LOG("GMCJDriver", pINFO) << "Generating a flux neutrino";

  bool ok = fFluxDriver->GenerateNext();
  if(!ok) {
     LOG("GMCJDriver", pWARN)
       << "Couldn't generate a flux neutrino - Returning NULL EventRecord";
     return 0;
  }
  int                    nupdg = fFluxDriver -> PdgCode  ();
  const TLorentzVector & nup4  = fFluxDriver -> Momentum ();
  const TLorentzVector & nux4  = fFluxDriver -> Position ();

  LOG("GMCJDriver", pINFO)
     << "\n [-] Generated flux neutrino: "
     << "\n  |----o PDG-code   : " << nupdg
     << "\n  |----o 4-momentum : " << print_utils::P4AsString(&nup4)
     << "\n  |----o 4-position : " << print_utils::X4AsString(&nux4);

  //-- get the list of detector material and the v path-length in each one,
  //   staring from nux4 and travelling along the direction of nup4

  const PathLengthList & pl = fGeomAnalyzer->ComputePathLengths(nux4, nup4);

  if(pl.size() == 0) {
     LOG("GMCJDriver", pWARN)
           << "Got an empty PathLengthList - Returning NULL EventRecord";
     return 0;
  }
  LOG("GMCJDriver", pINFO) << "Path lengths for flux neutrino direction: " << pl;

  //-- compute the total (scaled) interaction probabilities for each material
  //   (for the selected nu type

  LOG("GMCJDriver", pINFO)
             << "Computing interaction probabilities for each material";
  GEVGDriver * evgdriver = 0;
  double probsum=0;
  map<int,double> probm;

  PathLengthList::const_iterator pliter;
  for(pliter = pl.begin(); pliter != pl.end(); ++pliter) {

     int    mpdg = pliter->first;  // material PDG code
     double pL   = pliter->second; // density x path-length

     // find the GEVGDriver object that is handling the current init state
     InitialState init_state(mpdg, nupdg);
     LOG("GMCJDriver", pINFO) << "Current init state: " << init_state.AsString();
     evgdriver = fGPool->FindDriver(init_state);
     assert(evgdriver);

     // compute the interaction xsec and the interaction probability
     double xsec  = evgdriver->SumCrossSection(nup4);
     double prob  = this->PInt(xsec,pL); // interaction probability
     double probn = prob/fPmax; // normalized interaction probability

     LOG("GMCJDriver", pINFO)
         << "TotXSec = " << xsec/cm2 << " cm^2, "
         << "Prob = " << prob << ", Norm.Prob = " << 100*probn << "%";

     probsum += probn;
     probm.insert(map<int,double>::value_type(mpdg,probsum));

  } // materials crossed by neutrino direction

  double Pesc = 1-probsum;
  LOG("GMCJDriver", pINFO)
          << "The 'no interaction' probability is: " << 100*Pesc << " %";

  //-- select a detector material
  LOG("GMCJDriver", pINFO)
             << "Computing interaction probabilities for each material";

  RandomGen * rnd = RandomGen::Instance();
  double R = rnd->Random2().Rndm();
  LOG("GMCJDriver", pDEBUG) << "Rndm [0,1] = " << R;

  //if the neutrino does not interact go into a recursive mode until
  //an interaction does take place.
  if(R>=1-Pesc) {
     LOG("GMCJDriver", pINFO) << "Flux neutrino didn't interact - Retrying!";
     return this->GenerateEvent();
  }

  //an interaction does happen: select material
  LOG("GMCJDriver", pINFO)
                   << "The neutrino does interact - Selecting material";
  int tgtpdg = 0;
  map<int,double>::const_iterator probiter;
  for(probiter = probm.begin(); probiter != probm.end(); ++probiter) {
     double prob = probiter->second;
     if(R<prob) {
         tgtpdg = probiter->first;
         break;
     }
  }
  if(!tgtpdg) {
     LOG("GMCJDriver", pWARN)
       << "Couldn't select target material - Returning NULL EventRecord";
     return 0;
  }
  LOG("GMCJDriver", pINFO) << "Selected target material = " << tgtpdg;

  //-- find the GEVGDriver object that generates interactions for the given
  //   initial state (neutrino + target)
  InitialState init_state(tgtpdg, nupdg);
  LOG("GMCJDriver", pWARN)
      << "Searching GEVGDriver configured with: " << init_state.AsString();

  evgdriver = fGPool->FindDriver(init_state);
  if(!evgdriver) {
     LOG("GMCJDriver", pWARN)
       << "No GEVGDriver object found - Returning NULL EventRecord";
     return 0;
  }

  //-- ask the GEVGDriver object to select and generate an interaction for the
  //   selected initial state & neutrino 4-momentum
  LOG("GMCJDriver", pINFO)
              << "Asking the selected GEVGDriver object to generate an event";
  EventRecord * event = evgdriver->GenerateEvent(nup4);

  //-- generate an 'interaction position' in the selected material, along
  //   the direction of nup4
  LOG("GMCJDriver", pINFO)
                  << "Asking the geometry analyzer to generate a vertex";
  const TVector3 & vtx = fGeomAnalyzer->GenerateVertex(nux4, nup4, tgtpdg);

  //-- the GEVGDriver object generates events at (x=0,y=0,z=0,t=0) / shift the
  //   event record entries accrording to the selected interaction vertex

  TVector3 origin(nux4.X(), nux4.Y(), nux4.Z());
  origin-=vtx; // computes vector dr = origin - vtx

  double dL = origin.Mag();
  double c  = kLightSpeed /(units::meter/units::second);
  double dt = dL/c;

  LOG("GMCJDriver", pINFO)
        << "|vtx - origin|: dL = " << dL << " m, dt = " << dt << " sec";

  TLorentzVector vtx4(vtx.x(), vtx.y(), vtx.z(), nux4.T() + dt);
  event->ShiftVertex(vtx4);

  return event;
}
//___________________________________________________________________________
double GMCJDriver::PInt(double xsec, double path_length)
{
// tmp - currently pl (path-length) absorbs other factors

  return (xsec*path_length);
}
//___________________________________________________________________________
