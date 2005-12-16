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
#include <TSystem.h>
#include <TVector3.h>

#include "Algorithm/AlgConfigPool.h"
#include "Conventions/Constants.h"
#include "Conventions/Units.h"
#include "EVGCore/EventRecord.h"
#include "EVGDrivers/GMCJDriver.h"
#include "EVGDrivers/GEVGDriver.h"
#include "EVGDrivers/GEVGPool.h"
#include "EVGDrivers/GFluxI.h"
#include "EVGDrivers/GeomAnalyzerI.h"
#include "EVGDrivers/PathLengthList.h"
#include "Interaction/InitialState.h"
#include "Messenger/Messenger.h"
#include "Numerical/RandomGen.h"
#include "PDG/PDGCodeList.h"
#include "Utils/PrintUtils.h"
#include "Utils/XSecSplineList.h"

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
void GMCJDriver::UseSplines(bool useLogE)
{
  fUseSplines = true;
  fUseLogE    = useLogE;
}
//___________________________________________________________________________
void GMCJDriver::UseMaxPathLengths(string xml_filename)
{
// If you supply the maximum path lengths for all materials in your geometry
// (eg for ROOT/GEANT geometries by running GENIE#s gmxpl application - see
// $GENIE/src/stdapp/gMaxPathLengths.cxx - ) you can speed up the driver init
// time (quite a bit for complex geometries).

  fMaxPlXmlFilename = xml_filename;

  bool is_accessible = !(gSystem->AccessPathName(fMaxPlXmlFilename.c_str()));

  if (!is_accessible) fUseExtMaxPl = false;
  else                fUseExtMaxPl = true;
}
//___________________________________________________________________________
void GMCJDriver::AllowRecursiveMode(bool allow)
{
  LOG("GMCJDriver", pNOTICE)
        << "[Allow recursive mode] flag is set to: "
                             << utils::print::BoolAsYNString(allow);
  fAllowRecursMode = allow;
}
//___________________________________________________________________________
void GMCJDriver::FilterUnphysical(bool filter)
{
  LOG("GMCJDriver", pNOTICE)
        << "[Filter unphysical] flag is set to: "
                             << utils::print::BoolAsYNString(filter);
  fFilterUnphysical = filter;

  //if Configure() was run first configure all GEVGDrivers now
  if(fGPool) {
    GEVGPool::const_iterator giter;
    for(giter = fGPool->begin(); giter != fGPool->end(); ++giter) {
      GEVGDriver * driver = giter->second;
      driver->FilterUnphysical(filter);
    }
  }
}
//___________________________________________________________________________
void GMCJDriver::Initialize(void)
{
  fFluxDriver       = 0;  // <-- flux driver
  fGeomAnalyzer     = 0;  // <-- geometry driver
  fGPool            = 0;  // <-- pool of GEVGDriver event generation drivers
  fPmax             = 0;  // <-- (scaled) maximum interaction probability
  fMaxPlXmlFilename = ""; // <-- XML file with external path lengths
  fUseExtMaxPl      = false;
  fUseSplines       = false;

  // Go into recursive mode when it does not generate an event (neutrino
  // does not cross the detector, does not interact etc...) so that it never
  // returns NULL (except when in error)
  this->AllowRecursiveMode(true);

  // Allow the selected GEVGDriver to go into recursive mode and regenerate
  // an interaction that turns out to be unphysical.
  this->FilterUnphysical(true);

  // Force early initialization of singleton objects that are typically
  // would be initialized at their first use later on.
  // This is purely cosmetic and I do it to send the banner and some prolific
  // initialization printout at the top.
  assert( Messenger::Instance()     );
  assert( AlgConfigPool::Instance() );

  // Autoload splines (from the XML file pointed at the $GSPLOAD env. var.,
  // if the env. var. has been set);
  XSecSplineList * xspl = XSecSplineList::Instance();
  xspl->AutoLoad();
}
//___________________________________________________________________________
void GMCJDriver::Configure(void)
{
  LOG("GMCJDriver", pNOTICE)
                  << utils::print::PrintFramedMesg("Configuring GMCJDriver");

  //-- Ask the input GFluxI for a list of all neutrino-types it will be using
  LOG("GMCJDriver", pNOTICE)
                  << "Asking the flux driver for the list of flux particles";
  const PDGCodeList & nulist = fFluxDriver->FluxParticles();

  //-- Ask the ROOTGeomHandler for a list of all target materials
  LOG("GMCJDriver", pNOTICE)
            << "Asking the geometry driver for the list of target materials";
  const PDGCodeList & tgtlist = fGeomAnalyzer->ListOfTargetNuclei();

  //-- Print particle lists
  LOG("GMCJDriver", pNOTICE) << "Flux particles: "   << nulist;
  LOG("GMCJDriver", pNOTICE) << "Target materials: " << tgtlist;

  //-- Ask the ROOTGeomHandler for a the max. path length for each material
  //   in the list of all target materials (to be used in calculating the max
  //   interaction probability which will scale all interaction probabilities)

  PathLengthList plmax;

  if(fUseExtMaxPl) {
     LOG("GMCJDriver", pNOTICE)
               << "Loading external max path-length list for input geometry";
     plmax.LoadFromXml(fMaxPlXmlFilename);

  } else {
     LOG("GMCJDriver", pNOTICE)
         << "Asking the geometry driver to compute the max path-length list";
     plmax = fGeomAnalyzer->ComputeMaxPathLengths();
  }
  //-- Print maximum path lengths & neutrino energy
  LOG("GMCJDriver", pNOTICE) << "Maximum path length list: " << plmax;

  //-- Ask the input GFluxI for the max. neutrino energy (to compute Pmax)
  LOG("GMCJDriver", pNOTICE)
       << "Asking the flux driver for the maximum energy of flux neutrinos";
  double Emax = fFluxDriver->MaxEnergy();
  LOG("GMCJDriver", pNOTICE) << "Maximum flux neutrino energy = " << Emax;

  //-- Create all possible initial states and initialize/store a
  //   GEVGDriver driver object for each one of them.

  LOG("GMCJDriver", pDEBUG)
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

     LOG("GMCJDriver", pNOTICE)
       << "\n\n ---- Creating a GEVGDriver object configured for init-state: "
       << init_state.AsString() << " ----\n\n";

     GEVGDriver * evgdriver = new GEVGDriver;
     evgdriver->SetInitialState(init_state);
     evgdriver->FilterUnphysical(fFilterUnphysical);

     LOG("GMCJDriver", pDEBUG) << "Adding new GEVGDriver object to GEVGPool";
     fGPool->insert( GEVGPool::value_type(init_state.AsString(), evgdriver) );
   } // targets
  } // neutrinos

  LOG("GMCJDriver", pNOTICE)
             << "All necessary GEVGDriver object were pushed into GEVGPool\n";

  // If requested splines then coordinate spline creation from all GEVGDriver
  // objects pushed into GEVGPool. This will create all xsec splines needed
  // for all simulated processes involving the particles in the input flux
  // and geometry. Spline creation will be skipped for every spline that has
  // been pre-loaded into the the XSecSplineList

  if(fUseSplines) {
     LOG("GMCJDriver", pNOTICE) << "Creating cross section splines";

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

     // get the appropriate driver
     GEVGDriver * evgdriver = fGPool->FindDriver(init_state);

     double xsec_sum    = evgdriver->SumCrossSection(nup4); // sum{xsec}
     double path_length = plmax.PathLength(target_pdgc);    // max{L*density}

     double Pmax = this->PInt(xsec_sum, path_length);
     fPmax += Pmax;

     LOG("GMCJDriver", pNOTICE)
              << "Pmax[" << init_state.AsString()
                               << ", for given flux & geometry] = " << Pmax;
   } // targets
  } // neutrinos

  LOG("GMCJDriver", pNOTICE) << "Total Pmax[interaction] = "<< fPmax << "\n";
  LOG("GMCJDriver", pNOTICE) << "Finished configuring GMCJDriver\n\n";
}
//___________________________________________________________________________
EventRecord * GMCJDriver::GenerateEvent(void)
{
  LOG("GMCJDriver", pNOTICE) << "Generating next event...";

  //-- Generate a neutrino using the input GFluxI & get current pdgc/p4/x4

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

  LOG("GMCJDriver", pNOTICE)
     << "\n [-] Generated flux neutrino: "
     << "\n  |----o PDG-code   : " << nupdg
     << "\n  |----o 4-momentum : " << utils::print::P4AsString(&nup4)
     << "\n  |----o 4-position : " << utils::print::X4AsString(&nux4);

  //-- Get the list of detector material and the v path-length in each one,
  //   starting from nux4 and travelling along the direction of nup4

  const PathLengthList & pl = fGeomAnalyzer->ComputePathLengths(nux4, nup4);

  if(pl.size() == 0) {
     LOG("GMCJDriver", pERROR)
        << "Got an empty PathLengthList - No material found in geometry?";
     LOG("GMCJDriver", pERROR) << "Returning NULL EventRecord!";
     return 0;
  }

  //-- Check that everything is ok
  //-- If all pathlengths are 0 then the neutrino didn't cross the detector
  //   so do not bother continuing...

  LOG("GMCJDriver", pNOTICE) << "Path lengths for flux v direction: " << pl;

  if(pl.AreAllZero()) {
    LOG("GMCJDriver", pNOTICE)
              << "The flux v doesn't even enter the detector";

    if(fAllowRecursMode) {
       LOG("GMCJDriver", pNOTICE)
          << "In recursive mode - Attermting to regenerate the event...";
       return this->GenerateEvent(); // enter in reccursive mode...
    } else {
       LOG("GMCJDriver", pNOTICE)
          << "Recursive mode not allowed  - Returning NULL EventRecord!";
       return 0;
    }
  }

  //-- Compute the total (scaled) interaction probabilities for each
  //   material (for the selected neutrino type

  LOG("GMCJDriver", pNOTICE)
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

     // compute the interaction xsec and probability (if path-length>0)
     double xsec  = 0.; // sum of interaction xsecs for given init state
     double prob  = 0.; // interaction probability
     double probn = 0.; // normalized interaction probability
     if(pL>0.) {
        xsec  = evgdriver->SumCrossSection(nup4);
        prob  = this->PInt(xsec,pL);
        probn = prob/fPmax;
     }

     LOG("GMCJDriver", pINFO)
         << "TotXSec = " << xsec/cm2 << " cm^2, "
         << "Prob = " << prob << ", Norm.Prob = " << 100*probn << "%";

     probsum += probn;
     probm.insert(map<int,double>::value_type(mpdg,probsum));

  } // materials crossed by neutrino direction

  double Pesc = 1-probsum;
  LOG("GMCJDriver", pINFO)
          << "The 'no interaction' probability is: " << 100*Pesc << " %";

  //-- Select a detector material
  LOG("GMCJDriver", pINFO)
        << "Deciding whether the neutrino interacts and on which target";

  RandomGen * rnd = RandomGen::Instance();
  double R = rnd->Random2().Rndm();
  LOG("GMCJDriver", pDEBUG) << "Rndm [0,1] = " << R;

  //-- If the neutrino does not interact go into a recursive mode until
  //-- an interaction does take place.
  if(R>=1-Pesc) {
     LOG("GMCJDriver", pINFO) << "Flux neutrino didn't interact - Retrying!";

     if(fAllowRecursMode) {
         LOG("GMCJDriver", pINFO)
             << "In recursive mode - Attermting to regenerate the event...";
         return this->GenerateEvent(); // enter in reccursive mode...
      } else {
         LOG("GMCJDriver", pINFO)
            << "Recursive mode not allowed  - Returning NULL EventRecord!";
         return 0;
      }
  }

  //-- If an interaction does happen then select target material
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

  //-- Find the GEVGDriver object that generates interactions for the
  //   given initial state (neutrino + target)
  InitialState init_state(tgtpdg, nupdg);
  LOG("GMCJDriver", pWARN)
      << "Searching GEVGDriver configured with: " << init_state.AsString();

  evgdriver = fGPool->FindDriver(init_state);
  if(!evgdriver) {
     LOG("GMCJDriver", pWARN)
       << "No GEVGDriver object found - Returning NULL EventRecord";
     return 0;
  }

  //-- Ask the GEVGDriver object to select and generate an interaction for
  //   the selected initial state & neutrino 4-momentum
  LOG("GMCJDriver", pINFO)
              << "Asking the selected GEVGDriver object to generate an event";
  EventRecord * event = evgdriver->GenerateEvent(nup4);

  //-- Generate an 'interaction position' in the selected material, along
  //   the direction of nup4
  LOG("GMCJDriver", pINFO)
                  << "Asking the geometry analyzer to generate a vertex";
  const TVector3 & vtx = fGeomAnalyzer->GenerateVertex(nux4, nup4, tgtpdg);

  //-- The GEVGDriver object generates events at (x=0,y=0,z=0,t=0) / shift the
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
