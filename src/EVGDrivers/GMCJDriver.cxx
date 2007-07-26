//____________________________________________________________________________
/*
 Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         STFC, Rutherford Appleton Laboratory - May 25, 2005

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

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
#include "GHEP/GHepFlags.h"
#include "Interaction/InitialState.h"
#include "Messenger/Messenger.h"
#include "Numerical/RandomGen.h"
#include "Numerical/Spline.h"
#include "PDG/PDGUtils.h"
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
// (eg for ROOT/GEANT geometries they can be computed running GENIE's gmxpl 
// application, see $GENIE/src/stdapp/gMaxPathLengths.cxx ) you can speed up 
// the driver init phase by quite a bit (especially for complex geometries).

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
void GMCJDriver::FilterUnphysical(const TBits & unphysmask)
{
  fUnphysMask = unphysmask;

  //if Configure() was run first configure all GEVGDrivers now
  if(fGPool) {
    GEVGPool::const_iterator giter;
    for(giter = fGPool->begin(); giter != fGPool->end(); ++giter) {
      GEVGDriver * driver = giter->second;
      driver->FilterUnphysical(fUnphysMask);
    }
  }
}
//___________________________________________________________________________
void GMCJDriver::Initialize(void)
{
  fFluxDriver       = 0;  // <-- flux driver
  fGeomAnalyzer     = 0;  // <-- geometry driver
  fGPool            = 0;  // <-- pool of GEVGDriver event generation drivers
  fEmax             = 0;  // <-- maximum neutrino energy
  fPmax             = 0;  // <-- maximum interaction probability (prob scale)
  fMaxPlXmlFilename = ""; // <-- XML file with external path lengths
  fUseExtMaxPl      = false;
  fUseSplines       = false;
  fNFluxNeutrinos   = 0;  // <-- number of flux neutrinos thrown so far

  fSelTgtPdg        = 0;
  fCurEvt           = 0;
  fCurVtx.SetXYZT(0.,0.,0.,0.);

  // Go into recursive mode when it does not generate an event (neutrino
  // does not cross the detector, does not interact etc...) so that it never
  // returns NULL (except when in error)
  this->AllowRecursiveMode(true);

  // Allow the selected GEVGDriver to go into recursive mode and regenerate
  // an interaction that turns out to be unphysical.
  TBits unphysmask(GHepFlags::NFlags());
  unphysmask.ResetAllBits(false); 
  this->FilterUnphysical(unphysmask);

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

  // Clear the target and neutrino lists
  fNuList.clear();
  fTgtList.clear();

  // Clear the maximum path length list
  fMaxPathLengths.clear();
  fCurPathLengths.clear();
}
//___________________________________________________________________________
void GMCJDriver::Configure(void)
{
  LOG("GMCJDriver", pNOTICE)
                  << utils::print::PrintFramedMesg("Configuring GMCJDriver");

  //-- Get the list of neutrino types from the input flux driver and the list
  //   of target materials from the input geometry driver
  this->GetParticleLists();

  //-- Ask the input geometry driver to compute the max. path length for each
  //   material in the list of target materials (or load a precomputed list)
  this->GetMaxPathLengthList();

  //-- Ask the input GFluxI for the max. neutrino energy (to compute Pmax)
  this->GetMaxFluxEnergy();

  //-- Create all possible initial states and initialize/store a GEVGDriver
  //   driver object for each one of them.
  this->CreateGEVGDriverPool();

  // If requested splines then coordinate spline creation from all GEVGDriver
  // objects pushed into GEVGPool. This will create all xsec splines needed
  // for all simulated processes involving the particles in the input flux
  // and geometry. Spline creation will be skipped for every spline that has
  // been pre-loaded into the the XSecSplineList
  this->CreateXSecSplines();

  // Create cross section splines describing the total interaction xsec
  // for a given initial state (Create them by summing all xsec splines
  // for each possible initial state)
  this->CreateXSecSumSplines();

  // Compute the max. interaction probability to scale all interaction
  // probabilities to be computed by this driver
  this->ComputeMaxIntProb();

  LOG("GMCJDriver", pNOTICE) << "Finished configuring GMCJDriver\n\n";
}
//___________________________________________________________________________
void GMCJDriver::GetParticleLists(void)
{
  //-- get the list of flux neutrinos from the flux driver
  LOG("GMCJDriver", pNOTICE)
                    << "Asking the flux driver for its list of neutrinos";
  fNuList = fFluxDriver->FluxParticles();

  LOG("GMCJDriver", pNOTICE) << "Flux particles: " << fNuList;

  // gets the list of target materials from the geometry driver
  LOG("GMCJDriver", pNOTICE)
                  << "Asking the geometry driver for its list of targets";
  fTgtList = fGeomAnalyzer->ListOfTargetNuclei();

  LOG("GMCJDriver", pNOTICE) << "Target materials: " << fTgtList;
}
//___________________________________________________________________________
void GMCJDriver::GetMaxPathLengthList(void)
{
  if(fUseExtMaxPl) {
     LOG("GMCJDriver", pNOTICE)
               << "Loading external max path-length list for input geometry";
     fMaxPathLengths.LoadFromXml(fMaxPlXmlFilename);

  } else {
     LOG("GMCJDriver", pNOTICE)
         << "Asking the geometry driver to compute the max path-length list";
     fMaxPathLengths = fGeomAnalyzer->ComputeMaxPathLengths();
  }
  //-- Print maximum path lengths & neutrino energy
  LOG("GMCJDriver", pNOTICE)
                          << "Maximum path length list: " << fMaxPathLengths;
}
//___________________________________________________________________________
void GMCJDriver::GetMaxFluxEnergy(void)
{
  LOG("GMCJDriver", pNOTICE)
       << "Asking the flux driver for the maximum energy of flux neutrinos";
  fEmax = fFluxDriver->MaxEnergy();

  LOG("GMCJDriver", pNOTICE) << "Maximum flux neutrino energy = " << fEmax;
}
//___________________________________________________________________________
void GMCJDriver::CreateGEVGDriverPool(void)
{
  LOG("GMCJDriver", pDEBUG)
          << "Creating GEVGPool & adding a GEVGDriver object per init-state";

  if (fGPool) delete fGPool;
  fGPool = new GEVGPool;

  PDGCodeList::const_iterator nuiter;
  PDGCodeList::const_iterator tgtiter;

  for(nuiter = fNuList.begin(); nuiter != fNuList.end(); ++nuiter) {
   for(tgtiter = fTgtList.begin(); tgtiter != fTgtList.end(); ++tgtiter) {

     int target_pdgc   = *tgtiter;
     int neutrino_pdgc = *nuiter;

     InitialState init_state(target_pdgc, neutrino_pdgc);

     LOG("GMCJDriver", pNOTICE)
       << "\n\n ---- Creating a GEVGDriver object configured for init-state: "
       << init_state.AsString() << " ----\n\n";

     GEVGDriver * evgdriver = new GEVGDriver;
     evgdriver->Configure(init_state);
     evgdriver->FilterUnphysical(fUnphysMask);
     evgdriver->UseSplines(); // check if all splines needed are loaded

     LOG("GMCJDriver", pDEBUG) << "Adding new GEVGDriver object to GEVGPool";
     fGPool->insert( GEVGPool::value_type(init_state.AsString(), evgdriver) );
   } // targets
  } // neutrinos

  LOG("GMCJDriver", pNOTICE)
             << "All necessary GEVGDriver object were pushed into GEVGPool\n";
}
//___________________________________________________________________________
void GMCJDriver::CreateXSecSplines(void)
{
  if(!fUseSplines) return;

  LOG("GMCJDriver", pNOTICE) << "Creating cross section splines";

  PDGCodeList::const_iterator nuiter;
  PDGCodeList::const_iterator tgtiter;

  for(nuiter = fNuList.begin(); nuiter != fNuList.end(); ++nuiter){
     for(tgtiter = fTgtList.begin(); tgtiter != fTgtList.end(); ++tgtiter) {

       int target_pdgc   = *tgtiter;
       int neutrino_pdgc = *nuiter;

       InitialState init_state(target_pdgc, neutrino_pdgc);

       LOG("GMCJDriver", pINFO)
             << "Computing all splines needed for init-state: "
                                                  << init_state.AsString();
       GEVGDriver * evgdriver = fGPool->FindDriver(init_state);
       evgdriver->CreateSplines(-1,-1,fUseLogE);
     } // targets
  } // neutrinos
  LOG("GMCJDriver", pINFO) << "Finished creating cross section splines\n";
}
//___________________________________________________________________________
void GMCJDriver::CreateXSecSumSplines(void)
{
  LOG("GMCJDriver", pNOTICE)
         << "Creating {sum interaction cross section | init state} splines";

  GEVGPool::iterator diter;
  for(diter = fGPool->begin(); diter != fGPool->end(); ++diter) {
    string       init_state = diter->first;
    GEVGDriver * evgdriver  = diter->second;

    assert(evgdriver);

    LOG("GMCJDriver", pNOTICE)
             << "**** Summing xsec splines for init-state = " << init_state;

    Range1D_t rE = evgdriver->ValidEnergyRange();
    assert(fEmax<rE.max && fEmax>rE.min);

    // decide the energy range for the sum spline - extend the spline a little
    // bit above the maximum beam energy (but below the maximum valid energy)
    // to avoid the evaluation of the cubic spline around the viscinity of
    // knots with zero y values (although the GENIE Spline object handles it)
    double dE  = fEmax/10.;
    double min = rE.min;
    double max = (fEmax+dE < rE.max) ? fEmax+dE : rE.max;
    evgdriver->CreateXSecSumSpline(100,min,max,true);
  }
  LOG("GMCJDriver", pNOTICE)
        << "Finished summing all interaction xsec splines per initial state";
}
//___________________________________________________________________________
void GMCJDriver::ComputeMaxIntProb(void)
{
  LOG("GMCJDriver", pINFO)
    << "Computing the max. interaction probability (probability scale)";

  fPmax = 0; // maximum interaction probability
  TLorentzVector nup4(0,0,fEmax,fEmax);

  PDGCodeList::const_iterator nuiter;
  PDGCodeList::const_iterator tgtiter;

  for(nuiter = fNuList.begin(); nuiter != fNuList.end(); ++nuiter) {
   for(tgtiter = fTgtList.begin(); tgtiter != fTgtList.end(); ++tgtiter) {

     int target_pdgc   = *tgtiter;
     int neutrino_pdgc = *nuiter;

     InitialState init_state(target_pdgc, neutrino_pdgc);

     LOG("GMCJDriver", pINFO)
           << "Computing Pmax for init-state: " << init_state.AsString();

     // get the appropriate driver
     GEVGDriver * evgdriver = fGPool->FindDriver(init_state);

     double sxsec = evgdriver->XSecSumSpline()->Evaluate(fEmax); // sum{xsec}
     double plmax = fMaxPathLengths.PathLength(target_pdgc); // max{L*density}
     int    A     = pdg::IonPdgCodeToA(target_pdgc);

     double Pmax = this->PInt(sxsec, plmax, A);
     fPmax += Pmax;

     LOG("GMCJDriver", pNOTICE)
              << "Pmax[" << init_state.AsString()
                               << ", for given flux & geometry] = " << Pmax;
   } // targets
  } // neutrinos

  LOG("GMCJDriver", pNOTICE) << "Total Pmax[interaction] = "<< fPmax << "\n";
}
//___________________________________________________________________________
void GMCJDriver::InitEventGeneration(void)
{
  fCurPathLengths.clear();
  fSelTgtPdg = 0;
  fCurVtx.SetXYZT(0.,0.,0.,0.);
}
//___________________________________________________________________________
EventRecord * GMCJDriver::GenerateEvent(void)
{
  LOG("GMCJDriver", pNOTICE) << "Generating next event...";

  this->InitEventGeneration();

  //-- Generate flux neutrinos until you have got one crossing some geometry
  //   material
  while(fCurPathLengths.AreAllZero()) {
    // Generate a neutrino using the input GFluxI & get current pdgc/p4/x4
    assert( this->GenerateFluxNeutrino() );
    // Compute (pathLength x density x weight fraction) for all materials
    // in the input geometry, for the neutrino generated by the flux driver
    assert( this->ComputePathLengths() );
  }

  //-- The neutrino enters the detector - Select target material
  fSelTgtPdg = this->SelectTargetMaterial();

  //-- If the neutrino did not interact enter in recursive mode to regenerate
  //   it (if allowed)
  if(fSelTgtPdg==0) {
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

  //-- Ask the GEVGDriver object to select and generate an interaction and
  //   its kinematics for the selected initial state & neutrino 4-momentum
  this->GenerateEventKinematics();

  //-- Generate an 'interaction position' in the selected material (in the
  //   detector coord system), along the direction of nup4 & set it 
  this->GenerateVertex();

  //-- Set the event probability (probability for this event to happen given
  //   the detector setup & the selected flux neutrino)
  //   Note for users: 
  //   The above probability is stored at GHepRecord::Probability()
  //   For normalization purposes make sure that you take into account the
  //   GHepRecord::Weight() -if event generation is weighted-, and
  //   GFluxI::Weight() -if beam simulation is weighted-.
  this->ComputeEventProbability();

  return fCurEvt;
}
//___________________________________________________________________________
bool GMCJDriver::GenerateFluxNeutrino(void)
{
// Ask the neutrino flux driver to generate a flux neutrino and make sure
// that things look ok...
//
  LOG("GMCJDriver", pINFO) << "Generating a flux neutrino";

  bool ok = fFluxDriver->GenerateNext();
  if(!ok) {
     LOG("GMCJDriver", pFATAL)
              << "*** The flux driver couldn't generate a flux neutrino!!";
     return false;
  }

  fNFluxNeutrinos++;
  int                    nupdg = fFluxDriver -> PdgCode  ();
  const TLorentzVector & nup4  = fFluxDriver -> Momentum ();
  const TLorentzVector & nux4  = fFluxDriver -> Position ();

  LOG("GMCJDriver", pNOTICE)
     << "\n [-] Generated flux neutrino: "
     << "\n  |----o PDG-code   : " << nupdg
     << "\n  |----o 4-momentum : " << utils::print::P4AsString(&nup4)
     << "\n  |----o 4-position : " << utils::print::X4AsString(&nux4);

  if(nup4.Energy() > fEmax) {
     LOG("GMCJDriver", pFATAL)
       << "\n *** Flux driver error ***"
       << "\n Generated flux v with E = " << nup4.Energy() << " GeV"
       << "\n Max v energy (declared by flux driver) = " << fEmax << " GeV"
       << "\n My interaction probability scaling is invalidated!!";
     return false;
  }
  if(!fNuList.ExistsInPDGCodeList(nupdg)) {
     LOG("GMCJDriver", pFATAL)
       << "\n *** Flux driver error ***"
       << "\n Generated flux v with pdg = " << nupdg
       << "\n It does not belong to the declared list of flux neutrinos"
       << "\n I was not configured to handle this!!";
     return false;
  }
  return true;
}
//___________________________________________________________________________
bool GMCJDriver::ComputePathLengths(void)
{
// Ask the geometry driver to compute (pathLength x density x weight frac.)
// for all detector materials for the neutrino generated by the flux driver
// and make sure that things look ok...

  fCurPathLengths.clear();

  const TLorentzVector & nup4  = fFluxDriver -> Momentum ();
  const TLorentzVector & nux4  = fFluxDriver -> Position ();

  fCurPathLengths = fGeomAnalyzer->ComputePathLengths(nux4, nup4);

  LOG("GMCJDriver", pNOTICE) << fCurPathLengths;

  if(fCurPathLengths.size() == 0) {
     LOG("GMCJDriver", pFATAL)
       << "\n *** Geometry driver error ***"
       << "\n Got an empty PathLengthList - No material found in geometry?";
     return false;
  }

  if(fCurPathLengths.AreAllZero()) {
         LOG("GMCJDriver", pNOTICE)
                 << "current flux v doesn't cross any geometry material...";
  }
  return true;
}
//___________________________________________________________________________
int GMCJDriver::SelectTargetMaterial(void)
{
//-- Compute the total (scaled) interaction probabilities for each
//   material (for the selected neutrino type) and select one. If the
//   flux neutrino does not interact return 0

  LOG("GMCJDriver", pNOTICE)
       << "Computing relative interaction probabilities for each material";

  int                    nupdg = fFluxDriver->PdgCode();
  const TLorentzVector & nup4  = fFluxDriver->Momentum();

  map<int,double> probm;
  double probsum=0;
  PathLengthList::const_iterator pliter;
  for(pliter = fCurPathLengths.begin();
                            pliter != fCurPathLengths.end(); ++pliter) {
     int    mpdg  = pliter->first;            // material PDG code
     double pl    = pliter->second;           // density x path-length
     int    A     = pdg::IonPdgCodeToA(mpdg);
     double xsec  = 0.;                       // sum of xsecs for given init state
     double prob  = 0.;                       // interaction probability
     double probn = 0.;                       // normalized interaction probability

     // find the GEVGDriver object that is handling the current init state
     InitialState init_state(mpdg, nupdg);
     GEVGDriver * evgdriver = fGPool->FindDriver(init_state);
     if(!evgdriver) {
       LOG("GMCJDriver", pFATAL)
        << "No GEVGDriver object for init state: " << init_state.AsString();
       exit(1);
     }
     // compute the interaction xsec and probability (if path-length>0)
     if(pl>0.) {
        xsec  = evgdriver->XSecSum(nup4);
        prob  = this->PInt(xsec,pl,A);

        // scale the interaction probability to the maximum one so as not
        // to have to throw few billions of flux neutrinos before getting
        // an interaction...
        probn = prob/fPmax;
     }

     LOG("GMCJDriver", pNOTICE)
         << "tgt: " << mpdg << " -> TotXSec = "
         << xsec/units::cm2 << " cm^2, Norm.Prob = " << 100*probn << "%";

     probsum += probn;
     probm.insert(map<int,double>::value_type(mpdg,probsum));
  }

  double Pno = 1-probsum;
  LOG("GMCJDriver", pNOTICE)
        << "The 'no interaction' probability is: " << 100*Pno << " %";

  //-- Select a detector material
  LOG("GMCJDriver", pINFO)
        << "Deciding whether the neutrino interacts and on which target";

  RandomGen * rnd = RandomGen::Instance();
  double R = rnd->RndEvg().Rndm();
  LOG("GMCJDriver", pDEBUG) << "Rndm [0,1] = " << R;

  //-- Check whether the neutrino interacts or not
  if(R>=1-Pno) return 0;

  //-- If an interaction does happen then select target material
  LOG("GMCJDriver", pINFO)
                   << "The neutrino does interact - Selecting material";
  int tgtpdg = 0;
  map<int,double>::const_iterator probiter;
  for(probiter = probm.begin(); probiter != probm.end(); ++probiter) {
     double prob = probiter->second;
     if(R<prob) {
        tgtpdg = probiter->first;
        LOG("GMCJDriver", pINFO) << "Selected target material = " << tgtpdg;
        return tgtpdg;
     }
  }

  LOG("GMCJDriver", pFATAL)
           << "Could not select target material for an interacting neutrino";
  exit(1);
  return 0;
}
//___________________________________________________________________________
void GMCJDriver::GenerateEventKinematics(void)
{
  int                    nupdg = fFluxDriver->PdgCode();
  const TLorentzVector & nup4  = fFluxDriver->Momentum();

  //-- Find the GEVGDriver object that generates interactions for the
  //   given initial state (neutrino + target)
  InitialState init_state(fSelTgtPdg, nupdg);
  GEVGDriver * evgdriver = fGPool->FindDriver(init_state);
  if(!evgdriver) {
     LOG("GMCJDriver", pFATAL)
       << "No GEVGDriver object for init state: " << init_state.AsString();
     exit(1);
  }

  //-- Ask the GEVGDriver object to select and generate an interaction for
  //   the selected initial state & neutrino 4-momentum
  LOG("GMCJDriver", pNOTICE)
          << "Asking the selected GEVGDriver object to generate an event";
  fCurEvt = evgdriver->GenerateEvent(nup4);
}
//___________________________________________________________________________
void GMCJDriver::GenerateVertex(void)
{
  //-- Generate an 'interaction position' in the selected material, along
  //   the direction of nup4
  LOG("GMCJDriver", pINFO)
                  << "Asking the geometry analyzer to generate a vertex";

  const TLorentzVector & p4 = fFluxDriver->Momentum ();
  const TLorentzVector & x4 = fFluxDriver->Position ();

  const TVector3 & vtx = fGeomAnalyzer->GenerateVertex(x4, p4, fSelTgtPdg);

  TVector3 origin(x4.X(), x4.Y(), x4.Z());
  origin-=vtx; // computes vector dr = origin - vtx

  double dL = origin.Mag();
  double c  = kLightSpeed /(units::meter/units::second);
  double dt = dL/c;

  LOG("GMCJDriver", pINFO)
        << "|vtx - origin|: dL = " << dL << " m, dt = " << dt << " sec";

  fCurVtx.SetXYZT(vtx.x(), vtx.y(), vtx.z(), x4.T() + dt);

  fCurEvt->SetVertex(fCurVtx);
}
//___________________________________________________________________________
void GMCJDriver::ComputeEventProbability(void)
{
// Compute event probability for the given flux neutrino & detector geometry

  // interaction cross section
  // (convert from physical units -> cm^2)
  double xsec = fCurEvt->XSec() / units::cm2; 

  // path length in detector along v direction for specified target material
  // (convert from kgr/m2 to gr/cm2)
  PathLengthList::const_iterator pliter = fCurPathLengths.find(fSelTgtPdg);
  double path_length = pliter->second;
  path_length *= ((units::kilogram/units::m2)/(units::gram/units::cm2));

  // target material mass number
  // (in gr)
  int A = pdg::IonPdgCodeToA(fSelTgtPdg);

  // interaction probability
  double P = this->PInt(xsec, path_length, A);

  fCurEvt->SetProbability(P);
}
//___________________________________________________________________________
double GMCJDriver::PInt(double xsec, double path_length, int A)
{
  return (xsec*path_length)/A;
}
//___________________________________________________________________________

