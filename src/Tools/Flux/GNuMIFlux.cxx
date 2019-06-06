//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Robert Hatcher <rhatcher@fnal.gov>
         Fermi National Accelerator Laboratory

         Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Jun 27, 2008 - CA
   The first implementation of this concrete flux driver was first added in
   the development version 2.5.1
 @ Aug 26, 2008 - RH
   A start at getting GNuMIFlux working.  Major problems included the wrong 
   type for fLf_Ntype; a mismatch on capitalization for some of the branches
   (Fortran has things like "NWtFar", but in the root file it ended up at 
   "Nwtfar"); and a lack of the use of "Nimpwt" (importance weight) as a overall 
   weight beyond the x-y position weights. Also lowered the max energy from 
   200 GeV to 125 -- beam is 120 GeV so no neutrinos should be generated above 
   that.  Begin generalizing so that one can have off axis x-y weights; not so
   important for FarDet, but not ignorable for NearDet, especially in the case 
   of rock events or higher beam energies.
 @ Aug 29, 2008 - RH
   Expand  GNuMIFluxPassThroughInfo to have all elements of the ntuple.	Use
   copy of GNuMIFluxPassThroughInfo held by the GNuMIFlux as the repository of 
   the leave elements. This avoids a bunch of duplication and generally makes 
   the code more readable.
 @ Dec 15, 2008 - RH
   Progress on managing numi flux.  Make use of .h and .C files generated from 
   the ntuples themselves via MakeClass to automate setting branches and 
   creating  the proper leaves. Unfortunately the geant3 and geant4 GNuMI 
   ntuples have different formats.  All the most important variables are there 
   in each, but upper/lower case varies. A critical difference is the switch 
   from Float_t to Double_t which means one can't simply use the same structure 
   an change branch assignments. The common format GNuMIFluxPassThrough needs a 
   Copy() function defined for each type of ntuple.
   The g3 case has been implemented, but the g4 version is a skeleton. The retained 
   variables could do with a review to ensure that we're not carrying more than 
   necessary. The helper classes g3numi and g4numi srouce (.h + .C) live in a 
   subdirectory GNuMINtuple that needn't otherwise be exposed to any GENIE user
   -- the actual code gets into the library by having both #include "gXnumi.h" 
   and #include "gXnumi.C".
   The x-y reweight function has been integrated into GNuMIFluxPassThrough. This 
   has been tested against the fortran version (this set of code still retains 
   some of that testing code which should be purged in the next iteration).  
   Generally the values are comparable at the expected precision (g3 32-bit
   float ~2.5e-5 or better in general for the x-y weight, even better for the 
   neutrino energy). There are some large outliers when calculating x-y weights 
   for muon decays because the algorithm, as written, depends on taking the 
   difference in two largish numbers (e.g. ~6 - ~6 = ~0.002) which results in a 
   large loss of fractional error precision when starting with floats.
   Some interfaces have been put in place for handling coordinate transformations 
   between the detector and beam, but these are not fully implemented and/or tested.
 @ Mar 13, 2009 - RH
   Lots of changes. XML parsing of configuration file. Coord transformations.
   All exchanges current in *user* coords.
   XML config file might be in $GENIE/src/FluxDrivers/GNuMINtuple
 @ Mar 27, 2009 - RH 
   gNuMIExptEvGen expect flux driver method LoadBeamSimData to take two args 
   (flux file, det config) not just one it did previously.  Added second arg
   that then calls LoadConfig().  
   Also add bogus POT_1cycle() method gNuMIExptEvGen expects ... I'm not sure 
   what the function is exactly supposed to return so it's bogus but at least 
   now the  EvtGen builds.
 @ Apr 01, 2009 - RH 
   Call ScanForMaxWeight() in GenerateNext() automatically if the user hasn't 
   already done so. Comment out annoying debug messages deep in inner loop.
   When calculating a starting point of the neutrino ray using the flux window 
   vectors store it into fgX4 (beam coord position) NOT fgP4 (beam coord p4).
   Don't try to store -1 in a size_t variable (though only some versions of 
   gcc warn about this).  This was only relevant if the flux file was given 
   without any path (ie. no "/").
 @ Apr 02, 2009 - RH 
   Improved scheme for estimating maximum weight - no longer depend solely on 
   existing near/far weights, but calculate weights for some (configurable) 
   number of entries and apply a (configurable) fudge factor. This allows 
   off-axis detectors (NOvA-IPND) to get something more reasonable. 
   Lower reported maximum energy from 125 to 120; no reason really to fudge this 
   up as it just adds to the rejection fraction.
   Accept GXMLPATH or GXMLPATHS as specifying locations.
 @ Apr 03, 2009 - RH 
   Internalize End() condition Remove SetFilePOT function (intent not mappable 
   to GNuMI?). SetNumOfCycles() optional 2nd arg to allow immediate reuse of 
   entries.  Do ScanForMaxWeight() at the end of the config so that MaxEv will 
   be set *before* any generation of neutrinos (for GMCJDriver). Allow MaxEnergy 
   (fMaxEv) to be set during ScanForMaxWeight() if the scan finds a value (*fudge 
   factor) higher than previously set. New <enumax> in XML config allows user 
   to set estimated enu maximum	and the fudge factor to use during the scan 
   (1.0=use exactly scanned max). Remove GNuMIFluxXMLHelper::TokenizeString() 
   in favor of existing genie::utils::str::Split() which I didn't know about.
 @ Apr 10, 2009 - RH 
   Fix coord transform code so that unit conversion doesn't screw it up.
   Add PrintConfig() method for dumping current config/state.
 @ Apr 13, 2009 - RH 
   Generally make "meters" the default 'user' units -- genie expects this.  
   Allow re-use of ntuple entry (don't reset it until moving on).  
   Provide means of determining distance between ray origin and dk vertex.  
   First entry depends on random # (ie. not always first ntuple entry).
   Resetting unit scale w/out resetting window/transform should work. Best 
   current guess for g4numi unpacking; currently still some unset variables 
   in the passthrough class, but they don't look critical. Protection in x-y 
   weight calculation for case where parent particle came to a stop 
   (parentp==0) from Trish Vahle.  Remove POT_1cycle() method.
 @ Apr 14, 2009 - RH 
   Add public MoveToZ0(double) method for pushing ray origin to specified
   user coordinate z Automatically call MoveToZ0() if SetUpstreamZ() has been 
   called with sensible (abs(z) < 1.0e30) value. Split SetEntryReuse(int) 
   function off from SetNumOfCycles(). In our case SetNumOfCycles probably 
   is going to be deprecated. Rename GNuMIFluxPassThrougInfo::Copy() to 
   ::MakeCopy()   so that we  don't confuse the issue w/ TObject::Copy() which 
   has	completelydifferent symantics and to avoid a annoying compiler warning. 
   XML parsing for <upstreamz> and <reuse> tags.
 @ Apr 22, 2009 - RH 
   Spin off AddFile() method where one can try to determine how many POTs each 
   file represents. Change some vars to Long64_t.
 @ Apr 23, 2009 - RH 
   First attempt at proton-on-target accounting (POTs); should be right	for
   unweighted neutrinos but probably isn't for weighted ones. Moved some 
   generated entry info from GNuMIFlux class into the GNuMIFluxPassThroughInfo 
   class so that it can be passed out for users to record or use.  
   Some general cleanup and reordering.
 @ May 13, 2009 - RH 
   Calculate flux window area correctly zero out fSumWeight,fNNeutrinos,
   fAccumPOTs after weight scan. Initialize fNuTot,fFilePots rather than 
   accept random garbage. Initialize w/ SetUpstreamZ such that the default 
   is the flux window. 
   Make Print() signature look like ROOT's typical (const Option_t* opt=""). 
   Tweak printout formats.
 @ Jul 22, 2009 - RH 
   New FLUGG flux ntuples have neutrinos from Omega parents; x-y reweight 
   function needs to know about those as well.
 @ Aug 25, 2009 - CA
   Adapt code to use the new utils::xml namespace.
 @ Aug 28, 2009 - RH
   New XML tag <using_param_set> allows one configuration to include another.
 @ Sep 17, 2009 - RH
   Fix EffPOTsPerNu calculation so it doesn't get lost in integer arithmetic.
   If fMaxWeight is bumped print new value.
   Make debug class xypartials output-able to ostream; print as part of 
   GNuMIFluxPathThrough ostreaming if GNUMI_TEST_XY_WGT defined.
   Still looking for why flux files with low E cuts tend to walk the MaxWeight 
   up to very high values (and thus severly lower the efficiency).
 @ Sep 17, 2009 - RH
   Support "flugg" ntuples as well as "g3numi" and "g4numi".  
   Each ntuple file/tree has a slightly different layout (element sizes, 
   capitalization in the names, mix of elements).  All share sufficient basic 
   info in order to calculate flavor, p4, x4 and weight.
 @ Feb 04, 2010 - RH
   New methods: 
   GetEntryNumber() - current entry # in the TChain.
   GetFileList() - get list of expanded file names used in TChain.
   Also, initialize fFlugg pointer; tweak debug messages and ostream<< op.
 @ Feb 22, 2011 - JD
   Implemented dummy versions of the new GFluxI::Clear, GFluxI::Index and 
   GFluxI::GenerateWeighted methods needed for pre-generation of flux
   interaction probabilities in GMCJDriver. 
 @ Mar 14, 2014 - TD
   Prevent an infinite loop in GenerateNext() when the flux driver has not been
   properly configured by exiting within GenerateNext_weighted().

*/
//____________________________________________________________________________

#include <cstdlib>
#include <fstream>
#include <vector>
#include <sstream>
#include <cassert>
#include <climits>

#include "libxml/xmlmemory.h"
#include "libxml/parser.h"

#include "Framework/Utils/XmlParserUtils.h"
#include "Framework/Utils/StringUtils.h"

#include <TFile.h>
#include <TChain.h>
#include <TChainElement.h>
#include <TSystem.h>
#include <TStopwatch.h>

#include "Framework/Conventions/Units.h"
#include "Framework/Conventions/GBuild.h"

#include "Tools/Flux/GNuMIFlux.h"
#include "Tools/Flux/GNuMINtuple/g3numi.h"
#include "Tools/Flux/GNuMINtuple/g3numi.C"
#include "Tools/Flux/GNuMINtuple/g4numi.h"
#include "Tools/Flux/GNuMINtuple/g4numi.C"
#include "Tools/Flux/GNuMINtuple/flugg.h"
#include "Tools/Flux/GNuMINtuple/flugg.C"

#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGCodeList.h"
#include "Framework/Numerical/MathUtils.h"
#include "Framework/Utils/PrintUtils.h"
#include "Framework/Utils/UnitUtils.h"

using std::endl;

#include <vector>
#include <algorithm>
#include <iomanip>
#include "TRegexp.h"
#include "TString.h"

#include "Tools/Flux/GFluxDriverFactory.h"
FLUXDRIVERREG4(genie,flux,GNuMIFlux,genie::flux::GNuMIFlux)

#ifdef  GNUMI_TEST_XY_WGT
static genie::flux::xypartials gpartials;  // global one used by CalcEnuWgt()
#endif

using namespace genie;
using namespace genie::flux;

// declaration of helper class
namespace genie {
  namespace flux  {
    class GNuMIFluxXMLHelper {
    public:
       GNuMIFluxXMLHelper(GNuMIFlux* gnumi) : fVerbose(0), fGNuMI(gnumi) { ; }
      ~GNuMIFluxXMLHelper() { ; }
       bool LoadConfig(std::string cfg);

      // these should go in a more general package
       std::vector<double>   GetDoubleVector(std::string str);
       std::vector<long int> GetIntVector(std::string str);

    private:
      bool     LoadParamSet(xmlDocPtr&, std::string cfg);
      void     ParseParamSet(xmlDocPtr&, xmlNodePtr&);
      void     ParseBeamDir(xmlDocPtr&, xmlNodePtr&);
      void     ParseBeamPos(std::string);
      void     ParseRotSeries(xmlDocPtr&, xmlNodePtr&);
      void     ParseWindowSeries(xmlDocPtr&, xmlNodePtr&);
      void     ParseEnuMax(std::string);
      TVector3 AnglesToAxis(double theta, double phi, std::string units = "deg");
      TVector3 ParseTV3(const std::string& );

      int         fVerbose;  ///< how noisy to be when parsing XML
      // who to apply these changes to
      GNuMIFlux*  fGNuMI;

      // derived offsets/rotations
      TVector3    fBeamPosXML;
      TRotation   fBeamRotXML;
      TVector3    fFluxWindowPtXML[3];
    };
  }
}

ClassImp(GNuMIFluxPassThroughInfo)

// some nominal positions used in the original g3 ntuple
const TLorentzVector kPosCenterNearBeam(0.,0.,  1039.35,0.);
const TLorentzVector kPosCenterFarBeam (0.,0.,735340.00,0.);

//____________________________________________________________________________
GNuMIFlux::GNuMIFlux() :
  GFluxExposureI(genie::flux::kPOTs)
{
  this->Initialize();
}
//___________________________________________________________________________
GNuMIFlux::~GNuMIFlux()
{
  this->CleanUp();
}

//___________________________________________________________________________
double GNuMIFlux::GetTotalExposure() const
{
  // complete the GFluxExposureI interface
  return UsedPOTs();
}

//___________________________________________________________________________
long int GNuMIFlux::NFluxNeutrinos(void) const
{
  /// number of flux neutrinos ray generated so far 
  return fNNeutrinos; 
}

//___________________________________________________________________________
bool GNuMIFlux::GenerateNext(void)
{
// Get next (unweighted) flux ntuple entry on the specified detector location
//
  RandomGen* rnd = RandomGen::Instance();
  while ( true ) {
     // Check for end of flux ntuple
     bool end = this->End();
     if ( end ) {
       LOG("Flux", pNOTICE) << "GenerateNext signaled End() ";
       return false;
     }

     // Get next weighted flux ntuple entry
     bool nextok = this->GenerateNext_weighted();
     if ( fGenWeighted ) return nextok;
     if ( ! nextok ) continue;

     /* RWH - debug purposes
     if ( fNCycles == 0 ) {
       LOG("Flux", pNOTICE)
          << "Got flux entry: " << fIEntry
          << " - Cycle: "<< fICycle << "/ infinite";
     } else {
       LOG("Flux", pNOTICE)
          << "Got flux entry: "<< fIEntry
          << " - Cycle: "<< fICycle << "/"<< fNCycles;
     }
     */

     // Get fractional weight & decide whether to accept curr flux neutrino
     double f = this->Weight() / fMaxWeight;
     //LOG("Flux", pNOTICE)
     //   << "Curr flux neutrino fractional weight = " << f;
     if (f > 1.) {
       fMaxWeight = this->Weight() * fMaxWgtFudge; // bump the weight
       LOG("Flux", pERROR)
         << "** Fractional weight = " << f 
         << " > 1 !! Bump fMaxWeight estimate to " << fMaxWeight
         << PassThroughInfo();
     }
     double r = (f < 1.) ? rnd->RndFlux().Rndm() : 0;
     bool accept = ( r < f );
     if ( accept ) {

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
       LOG("Flux", pNOTICE)
         << "Generated beam neutrino: "
         << "\n pdg-code: " << fCurEntry->fgPdgC
         << "\n p4: " << utils::print::P4AsShortString(&(fCurEntry->fgP4))
         << "\n x4: " << utils::print::X4AsString(&(fCurEntry->fgX4))
         << "\n p4: " << utils::print::P4AsShortString(&(fCurEntry->fgP4User))
         << "\n x4: " << utils::print::X4AsString(&(fCurEntry->fgX4User));
#endif

       fWeight = 1.;
       return true;
     }

     //LOG("Flux", pNOTICE)
     //  << "** Rejecting current flux neutrino based on the flux weight only";
  }
  return false;
}
//___________________________________________________________________________
bool GNuMIFlux::GenerateNext_weighted(void)
{
// Get next (weighted) flux ntuple entry on the specified detector location
//

  // Check whether a flux ntuple has been loaded
  if ( ! fG3NuMI && ! fG4NuMI && ! fFlugg ) {
     LOG("Flux", pFATAL)
          << "The flux driver has not been properly configured";
     //return false; //  don't do this - creates an infinite loop!
     exit(1);	
  }

  // Reuse an entry?
  //std::cout << " ***** iuse " << fIUse << " nuse " << fNUse
  //          << " ientry " << fIEntry << " nentry " << fNEntries
  //          << " icycle " << fICycle << " ncycle " << fNCycles << std::endl;
  if ( fIUse < fNUse && fIEntry >= 0 ) {
    // Reuse this entry
    fIUse++;
  } else {
    // Reset previously generated neutrino code / 4-p / 4-x
    this->ResetCurrent();
    // Move on, read next flux ntuple entry
    fIEntry++;
    if ( fIEntry >= fNEntries ) {
      // Ran out of entries @ the current cycle of this flux file
      // Check whether more (or infinite) number of cycles is requested
      if ( fICycle < fNCycles || fNCycles == 0 ) {
        fICycle++;
        fIEntry=0;
      } else {
        LOG("Flux", pWARN)
          << "No more entries in input flux neutrino ntuple, cycle "
          << fICycle << " of " << fNCycles;
        fEnd = true;
        // assert(0);
        return false;	
      }
    }
    
    if ( fG3NuMI ) {
      fG3NuMI->GetEntry(fIEntry); 
      fCurEntry->MakeCopy(fG3NuMI); 
    } else if ( fG4NuMI ) { 
      fG4NuMI->GetEntry(fIEntry); 
      fCurEntry->MakeCopy(fG4NuMI); 
    } else if ( fFlugg ) { 
      fFlugg->GetEntry(fIEntry); 
      fCurEntry->MakeCopy(fFlugg); 
    } else {
      LOG("Flux", pERROR) << "No ntuple configured";
      fEnd = true;
      //assert(0);
      return false;	
    }
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("Flux",pDEBUG) 
    << "got " << fNNeutrinos << " new fIEntry " << fIEntry 
    << " evtno " << fCurEntry->evtno;
#endif

    fIUse = 1; 
    fCurEntry->pcodes = 0;  // fetched entry has geant codes
    fCurEntry->units  = 0;  // fetched entry has original units

    // Convert the current gnumi neutrino flavor mode into a neutrino pdg code
    // Also convert other particle codes in GNuMIFluxPassThroughInfo to PDG
    fCurEntry->ConvertPartCodes();
    // here we might want to do flavor oscillations or simple mappings
    fCurEntry->fgPdgC = fCurEntry->ntype;
  }

  // Check neutrino pdg against declared list of neutrino species declared
  // by the current instance of the NuMI neutrino flux driver.
  // No undeclared neutrino species will be accepted at this point as GENIE
  // has already been configured to handle the specified list.
  // Make sure that the appropriate list of flux neutrino species was set at
  // initialization via GNuMIFlux::SetFluxParticles(const PDGCodeList &)

  // update the # POTs, number of neutrinos 
  // do this HERE (before rejecting flavors that users might be weeding out)
  // in order to keep the POT accounting correct.  This allows one to get
  // the right normalization for generating only events from the intrinsic
  // nu_e entries.
  fAccumPOTs += fEffPOTsPerNu / fMaxWeight;
  fNNeutrinos++;

  if ( ! fPdgCList->ExistsInPDGCodeList(fCurEntry->fgPdgC) ) {
     /// user might modify list via SetFluxParticles() in order to reject certain
     /// flavors, even if they're found in the file.  So don't make a big fuss.
     /// Spit out a single message and then stop reporting that flavor as problematic.
     int badpdg = fCurEntry->fgPdgC;
     if ( ! fPdgCListRej->ExistsInPDGCodeList(badpdg) ) {
       fPdgCListRej->push_back(badpdg);
       LOG("Flux", pWARN)
         << "Encountered neutrino specie (" << badpdg 
         << " pcodes=" << fCurEntry->pcodes << ")"
         << " that wasn't in SetFluxParticles() list, "
         << "\nDeclared list of neutrino species: " << *fPdgCList;
     }
     return false;	
  }

  // Update the curr neutrino weight and energy

  // Check current neutrino energy against the maximum flux neutrino energy 
  // declared by the current instance of the NuMI neutrino flux driver.
  // No flux neutrino exceeding that maximum energy will be accepted at this 
  // point as that maximum energy has already been used for normalizing the
  // interaction probabilities.
  // Make sure that the appropriate maximum flux neutrino energy was set at
  // initialization via GNuMIFlux::SetMaxEnergy(double Ev)

  fCurEntry->fgX4 = fFluxWindowBase;

  double Ev = 0;
  double& wgt_xy = fCurEntry->fgXYWgt;
  switch ( fUseFluxAtDetCenter ) {
  case -1:  // near detector
    wgt_xy   = fCurEntry->nwtnear;
    Ev       = fCurEntry->nenergyn;
    break;
  case +1:  // far detector
    wgt_xy   = fCurEntry->nwtfar;
    Ev       = fCurEntry->nenergyf;
    break;
  default:  // recalculate on x-y window
    RandomGen * rnd = RandomGen::Instance();
    fCurEntry->fgX4 += ( rnd->RndFlux().Rndm()*fFluxWindowDir1 +
                         rnd->RndFlux().Rndm()*fFluxWindowDir2   );
    fCurEntry->CalcEnuWgt(fCurEntry->fgX4,Ev,wgt_xy);
    break;
  }

  if (Ev > fMaxEv) {
     LOG("Flux", pWARN)
          << "Flux neutrino energy exceeds declared maximum neutrino energy"
          << "\nEv = " << Ev << "(> Ev{max} = " << fMaxEv << ")";
  }

  // Set the current flux neutrino 4-momentum
  // this is in *beam* coordinates
  fgX4dkvtx = TLorentzVector( fCurEntry->vx,
                              fCurEntry->vy,
                              fCurEntry->vz, 0.);
  // don't use TLorentzVector here for Mag() due to - on metric
  // normalize direction to 1.0
  TVector3 dirNu = (fCurEntry->fgX4.Vect() - fgX4dkvtx.Vect()).Unit();
  fCurEntry->fgP4.SetPxPyPzE( Ev*dirNu.X(), 
                              Ev*dirNu.Y(),
                              Ev*dirNu.Z(), Ev);

  // calculate the weight, potentially includes effect from tilted window
  // must be done *after* neutrino direction is determined
  fWeight = fCurEntry->nimpwt * fCurEntry->fgXYWgt;  // full weight
  if ( fApplyTiltWeight ) {
    double tiltwgt = dirNu.Dot( FluxWindowNormal() );
    fWeight *= TMath::Abs( tiltwgt );
  }

  // update sume of weights
  fSumWeight += this->Weight();

  // Set the current flux neutrino 4-position, direction in user coord
  Beam2UserP4(fCurEntry->fgP4,fCurEntry->fgP4User);
  Beam2UserPos(fCurEntry->fgX4,fCurEntry->fgX4User);

  // if desired, move to user specified user coord z
  if ( TMath::Abs(fZ0) < 1.0e30 ) this->MoveToZ0(fZ0);

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("Flux", pINFO)
    << "Generated neutrino: " << fIEntry << " " << fCurEntry->evtno
    << " nenergyn " << fCurEntry->nenergyn
    << "\n pdg-code: " << fCurEntry->fgPdgC
    << "\n p4 beam: " << utils::print::P4AsShortString(&fCurEntry->fgP4)
    << "\n x4 beam: " << utils::print::X4AsString(&fCurEntry->fgX4)
    << "\n p4 user: " << utils::print::P4AsShortString(&(fCurEntry->fgP4User))
    << "\n x4 user: " << utils::print::X4AsString(&(fCurEntry->fgX4User));
#endif
  if ( Ev > fMaxEv ) {
    LOG("Flux", pFATAL)
      << "Generated neutrino had E_nu = " << Ev << " > " << fMaxEv 
      << " maximum ";
    assert(0);
  }

  return true;
}
//___________________________________________________________________________
double GNuMIFlux::GetDecayDist() const
{
  // return distance (user units) between dk point and start position
  // these are in beam units
  TVector3 x3diff = fCurEntry->fgX4.Vect() - fgX4dkvtx.Vect();
  return x3diff.Mag() * fLengthScaleB2U;
}
//___________________________________________________________________________
void GNuMIFlux::MoveToZ0(double z0usr)
{
  // move ray origin to specified user z0
  // move beam coord entry correspondingly

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("Flux", pNOTICE)
    << "MoveToZ0 (z0usr=" << z0usr << ") before:"
    << "\n p4 user: " << utils::print::P4AsShortString(&(fCurEntry->fgP4User))
    << "\n x4 user: " << utils::print::X4AsString(&(fCurEntry->fgX4User));
#endif

  double pzusr    = fCurEntry->fgP4User.Pz();
  if ( TMath::Abs(pzusr) < 1.0e-30 ) {
    // neutrino is moving almost entirely in x-y plane
    LOG("Flux", pWARN)
      << "MoveToZ0(" << z0usr << ") not possible due to pz_usr (" << pzusr << ")";
    return;
  }

  double scale = (z0usr - fCurEntry->fgX4User.Z()) / pzusr; 
  fCurEntry->fgX4User += (scale*fCurEntry->fgP4User);
  fCurEntry->fgX4     += ((fLengthScaleU2B*scale)*fCurEntry->fgP4);
  // this scaling works for distances, but not the time component
  fCurEntry->fgX4.SetT(0);
  fCurEntry->fgX4User.SetT(0);

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("Flux", pNOTICE)
    << "MoveToZ0 (" << z0usr << ") after:"
    << "\n x4 user: " << utils::print::X4AsString(&(fCurEntry->fgX4User));
#endif

}

//___________________________________________________________________________
void GNuMIFlux::CalcEffPOTsPerNu()
{
  // do this if flux window changes or # of files changes

  if (!fNuFluxTree) return;  // not yet fully configured

  // effpots = mc_pots * (wgtfunction-area) / window-area / wgt-max-est
  //   wgtfunction-area = pi * radius-det-element^2 = pi * (100.cm)^2

  // this should match what is used in the CalcEnuWgt()
  const double kRDET = 100.0;   // set to flux per 100 cm radius
  const double kRDET2 = kRDET * kRDET;
  double flux_area = fFluxWindowDir1.Vect().Cross(fFluxWindowDir2.Vect()).Mag();
  LOG("Flux",pNOTICE) << "in CalcEffPOTsPerNu, area = " << flux_area;

  if ( flux_area < 1.0e-30 ) {
    LOG("Flux", pWARN)
          << "CalcEffPOTsPerNu called with flux window area effectively zero";
    flux_area = 1;
  }
  double area_ratio = TMath::Pi() * kRDET2 / flux_area;
  fEffPOTsPerNu = area_ratio * ( (double)fFilePOTs / (double)fNEntries );
}

//___________________________________________________________________________
double GNuMIFlux::UsedPOTs(void) const
{
// Compute current number of flux POTs

  if (!fNuFluxTree) {
     LOG("Flux", pWARN)
          << "The flux driver has not been properly configured";
     return 0;	
  }
  return fAccumPOTs;
}

//___________________________________________________________________________
double GNuMIFlux::POT_curr(void) { 
  // RWH: Not sure what POT_curr is supposed to represent I'll guess for
  // now that that it means what I mean by UsedPOTs().
  return UsedPOTs(); 
}

//___________________________________________________________________________
void GNuMIFlux::LoadBeamSimData(const std::vector<std::string>& patterns,
                                const std::string&              config )
{
// Loads in a gnumi beam simulation root file (converted from hbook format)
// into the GNuMIFlux driver.

  bool found_cfg = this->LoadConfig(config);
  if ( ! found_cfg ) {
    LOG("Flux", pFATAL) 
      << "LoadBeamSimData could not find XML config \"" << config << "\"\n";
    exit(1);
  }

  fNuFluxFilePatterns = patterns;
  std::vector<int> nfiles_from_pattern;

  // create a (sorted) set of file names
  // this also helps ensure that the same file isn't listed multiple times
  std::set<std::string> fnames;

  LOG("Flux", pINFO) << "LoadBeamSimData was passed " << patterns.size()
                     << " patterns";

  for (size_t ipatt = 0; ipatt < patterns.size(); ++ipatt ) {
    string pattern = patterns[ipatt];
    nfiles_from_pattern.push_back(0);

    LOG("Flux", pNOTICE)
      << "Loading gnumi flux tree from ROOT file pattern ["
      << std::setw(3) << ipatt << "] \"" << pattern << "\"";

    // !WILDCARD only works for file name ... NOT directory
    string dirname = gSystem->UnixPathName(gSystem->WorkingDirectory());
    size_t slashpos = pattern.find_last_of("/");
    size_t fbegin;
    if ( slashpos != std::string::npos ) {
      dirname = pattern.substr(0,slashpos);
      LOG("Flux", pINFO) << "Look for flux using directory " << dirname;
      fbegin = slashpos + 1;
    } else { fbegin = 0; }

    void* dirp = gSystem->OpenDirectory(gSystem->ExpandPathName(dirname.c_str()));
    if ( dirp ) {
      std::string basename = 
      pattern.substr(fbegin,pattern.size()-fbegin);
      TRegexp re(basename.c_str(),kTRUE);
      const char* onefile;
      while ( ( onefile = gSystem->GetDirEntry(dirp) ) ) {
        TString afile = onefile;
        if ( afile=="." || afile==".." ) continue;
        if ( basename!=afile && afile.Index(re) == kNPOS ) continue;
        std::string fullname = dirname + "/" + afile.Data();
        fnames.insert(fullname);
        nfiles_from_pattern[ipatt]++;
      }
      gSystem->FreeDirectory(dirp);
    } // legal directory
  } // loop over patterns

  size_t indx = 0;
  std::set<string>::const_iterator sitr = fnames.begin();
  for ( ; sitr != fnames.end(); ++sitr, ++indx ) {
    string filename = *sitr;
    //std::cout << "  [" << std::setw(3) << indx << "]  \"" 
    //          << filename << "\"" << std::endl;
    bool isok = true; 
    // this next test only works for local files, so we can't do that
    // if we want to stream via xrootd
    // ! (gSystem->AccessPathName(filename.c_str()));
    if ( isok ) {
      // open the file to see what it contains
      // h10    => g3numi _or_ flugg
      // nudata => g4numi
      // for now distinguish between g3numi/flugg using file name
      TFile* tf = TFile::Open(filename.c_str(),"READ");
      int isflugg = ( filename.find("flugg") != string::npos ) ? 1 : 0;
      const std::string tnames[] = { "h10", "nudata" };
      const std::string gnames[] = { "g3numi","g4numi","flugg","g4flugg"};
      for (int j = 0; j < 2 ; ++j ) { 
        TTree* atree = (TTree*)tf->Get(tnames[j].c_str());
        if ( atree ) {
          const std::string tname_this = tnames[j];
          const std::string gname_this = gnames[j+2*isflugg];
          // create chain if none exists
          if ( ! fNuFluxTree ) {
            this->SetTreeName(tname_this);
            fNuFluxTree = new TChain(fNuFluxTreeName.c_str());
            fNuFluxGen = gname_this;
            // here we should scan for estimated POTs/file
            // also maybe the check on max weight
          }
          // sanity check for mixing g3/g4/flugg files
          if ( fNuFluxTreeName !=  tname_this ||
               fNuFluxGen      !=  gname_this    ) {
            LOG("Flux", pFATAL) 
              << "Inconsistent flux file types\n"
              << "The input gnumi flux file \"" << filename 
              << "\"\ncontains a '" << tname_this << "' " << gname_this
              << "numi ntuple " 
              << "but a '" << fNuFluxTreeName << "' " << fNuFluxGen
              << " numi ntuple has alread been seen in the chain";
            exit(1);
          } // sanity mix/match g3/g4
            // add the file to the chain
          this->AddFile(atree,filename);
        } // found a tree
      } // loop over either g3 or g4
      tf->Close();
      delete tf;
    } // loop over tree type
  } // loop over sorted file names

  if ( fNuFluxTreeName == "" ) {
    LOG("Flux", pFATAL)
     << "The input gnumi flux file doesn't exist! Initialization failed!";
    exit(1);
  }
  if ( fNuFluxGen == "g3numi" ) fG3NuMI = new g3numi(fNuFluxTree);
  if ( fNuFluxGen == "g4numi" ) fG4NuMI = new g4numi(fNuFluxTree);
  if ( fNuFluxGen == "flugg"  ) fFlugg  = new flugg(fNuFluxTree);

  // this will open all files and read header!!
  fNEntries = fNuFluxTree->GetEntries();

  if ( fNEntries == 0 ) {
    LOG("Flux", pERROR)
      << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!";
    LOG("Flux", pERROR)
      << "Loaded flux tree contains " <<  fNEntries << " entries";
    LOG("Flux", pERROR)
      << "Was passed the file patterns: ";
    for (size_t ipatt = 0; ipatt < patterns.size(); ++ipatt ) {
      string pattern = patterns[ipatt];
      LOG("Flux", pERROR)
        << "  [" << std::setw(3) << ipatt <<"] " << pattern;
    }
    LOG("Flux", pERROR)
      << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!";
  } else {
    LOG("Flux", pNOTICE)
      << "Loaded flux tree contains " <<  fNEntries << " entries"
      << " from " << fnames.size() << " unique files";
    for (size_t ipatt = 0; ipatt < patterns.size(); ++ipatt ) {
      string pattern = patterns[ipatt];
      LOG("Flux", pINFO)
        << " pattern: " << pattern << " yielded "
        << nfiles_from_pattern[ipatt] << " files";
    }
  }

  // we have a file we can work with
  if (!fDetLocIsSet) {
     LOG("Flux", pERROR)
       << "LoadBeamSimData left detector location unset";
  }
  if (fMaxWeight<=0) {
     LOG("Flux", pINFO)
       << "Run ScanForMaxWeight() as part of LoadBeamSimData";
     this->ScanForMaxWeight();	
  }

  // current ntuple cycle # (flux ntuples may be recycled)
  fICycle =  0;
  // pick a starting entry index [0:fNEntries-1]
  // pretend we just used up the the previous one
  RandomGen* rnd = RandomGen::Instance();
  fIUse   =  9999999;
  fIEntry = rnd->RndFlux().Integer(fNEntries) - 1;
  
  // don't count things we used to estimate max weight
  fSumWeight  = 0;
  fNNeutrinos = 0;
  fAccumPOTs  = 0;

  LOG("Flux",pNOTICE) << "about to CalcEffPOTsPerNu";
  this->CalcEffPOTsPerNu();
  
}
//___________________________________________________________________________
void GNuMIFlux::GetBranchInfo(std::vector<std::string>& branchNames,
                              std::vector<std::string>& branchClassNames,
                              std::vector<void**>&      branchObjPointers)
{
  // allow flux driver to report back current status and/or ntuple entry 
  // info for possible recording in the output file by supplying
  // the class name, and a pointer to the object that will be filled
  // as well as a suggested name for the branch.
  
  branchNames.push_back("flux");
  branchClassNames.push_back("genie::flux::GNuMIFluxPassThroughInfo");
  branchObjPointers.push_back((void**)&fCurEntry);

}
TTree* GNuMIFlux::GetMetaDataTree() { return 0; } // there is none

//___________________________________________________________________________
void GNuMIFlux::ScanForMaxWeight(void)
{
  if (!fDetLocIsSet) {
     LOG("Flux", pERROR)
       << "Specify a detector location before scanning for max weight";
     return;	
  }

  // scan for the maximum weight
  int ipos_estimator = fUseFluxAtDetCenter;
  if ( ipos_estimator == 0 ) {
    // within 100m of a known point?
    double zbase = fFluxWindowBase.Z();
    if ( TMath::Abs(zbase-103648.837) < 10000. ) ipos_estimator = -1; // use NearDet
    if ( TMath::Abs(zbase-73534000. ) < 10000. ) ipos_estimator = +1; // use FarDet
  }
  if ( ipos_estimator != 0 ) {

    //// one can't really be sure which Nwtfar/Nwtnear this refers to
    //// some gnumi files have "NOvA" weights
    const char* ntwgtstrv[4] = { "Nimpwt*Nwtnear", 
                                 "Nimpwt*Nwtfar",
                                 "Nimpwt*NWtNear[0]",
                                 "Nimpwt*NWtFar[0]"  };
    int strindx = 0;
    if ( ipos_estimator > 0 ) strindx = 1;
    if ( fG4NuMI ) strindx += 2;
    // set upper limit on how many entries to scan
    Long64_t nscan = TMath::Min(fNEntries,200000LL);
    
    fNuFluxTree->Draw(ntwgtstrv[strindx],"","goff",nscan);
    //std::cout << " Draw \"" << ntwgtstrv[strindx] << "\"" << std::endl;
    //std::cout << " SelectedRows " << fNuFluxTree->GetSelectedRows()
    //          << " V1 " << fNuFluxTree->GetV1() << std::endl;

    Long64_t idx = TMath::LocMax(fNuFluxTree->GetSelectedRows(),
                                 fNuFluxTree->GetV1());
    //std::cout << "idx " << idx << " of " << fNuFluxTree->GetSelectedRows() << std::endl;
    fMaxWeight = fNuFluxTree->GetV1()[idx];
    LOG("Flux", pNOTICE) << "Maximum flux weight from Nwt in ntuple = " 
                         << fMaxWeight;
    if ( fMaxWeight <= 0 ) {
      LOG("Flux", pFATAL) << "Non-positive maximum flux weight!";
      exit(1);
    }
  }
  // the above works only for things close to the MINOS stored weight
  // values.  otherwise we need to work out our own estimate.
  double wgtgenmx = 0, enumx = 0;
  TStopwatch t;
  t.Start();
  for (int itry=0; itry < fMaxWgtEntries; ++itry) {
    this->GenerateNext_weighted();
    double wgt = this->Weight();
    if ( wgt > wgtgenmx ) wgtgenmx = wgt;
    double enu = fCurEntry->fgP4.Energy();
    if ( enu > enumx ) enumx = enu;
  }
  t.Stop();
  t.Print("u");
  LOG("Flux", pNOTICE) << "Maximum flux weight for spin = " 
                       << wgtgenmx << ", energy = " << enumx
                       << " (" << fMaxWgtEntries << ")";

  if (wgtgenmx > fMaxWeight ) fMaxWeight = wgtgenmx;
  // apply a fudge factor to estimated weight
  fMaxWeight *= fMaxWgtFudge;
  // adjust max energy?
  if ( enumx*fMaxEFudge > fMaxEv ) {
    LOG("Flux", pNOTICE) << "Adjust max: was=" << fMaxEv
                         << " now " << enumx << "*" << fMaxEFudge
                         << " = " << enumx*fMaxEFudge;
    fMaxEv = enumx * fMaxEFudge;
  }

  LOG("Flux", pNOTICE) << "Maximum flux weight = " << fMaxWeight 
                       << ", energy = " << fMaxEv;

}
//___________________________________________________________________________
void GNuMIFlux::SetMaxEnergy(double Ev)
{
  fMaxEv = TMath::Max(0.,Ev);

  LOG("Flux", pINFO)
    << "Declared maximum flux neutrino energy: " << fMaxEv;
}
//___________________________________________________________________________
void GNuMIFlux::SetEntryReuse(long int nuse)
{
// With nuse > 1 then the same entry in the file is used "nuse" times
// before moving on to the next entry in the ntuple

  fNUse    = TMath::Max(1L, nuse);
}
//___________________________________________________________________________
void GNuMIFlux::SetTreeName(string name)
{
  fNuFluxTreeName = name;
}
//___________________________________________________________________________
void GNuMIFlux::UseFluxAtNearDetCenter(void)
{
  SetLengthUnits(genie::utils::units::UnitFromString("m"));
  fFluxWindowBase      = kPosCenterNearBeam;
  fFluxWindowDir1      = TLorentzVector();   // no extent
  fFluxWindowDir2      = TLorentzVector();
  fUseFluxAtDetCenter  = -1;
  fDetLocIsSet         = true;
}
//___________________________________________________________________________
void GNuMIFlux::UseFluxAtFarDetCenter(void)
{
  SetLengthUnits(genie::utils::units::UnitFromString("m"));
  fFluxWindowBase      = kPosCenterFarBeam;
  fFluxWindowDir1      = TLorentzVector();  // no extent
  fFluxWindowDir2      = TLorentzVector();
  fUseFluxAtDetCenter  = +1;
  fDetLocIsSet         = true;
}
//___________________________________________________________________________
bool GNuMIFlux::SetFluxWindow(GNuMIFlux::StdFluxWindow_t stdwindow, double padding)
{
  // set some standard flux windows
  // rwh:  should also set detector coord transform
  // rwh:  padding allows add constant padding to pre-existing set
  double padbeam = padding / fLengthScaleB2U;  // user might set different units
  LOG("Flux",pNOTICE)
    << "SetBeamFluxWindow " << (int)stdwindow << " padding " << padbeam << " cm";


  switch ( stdwindow ) {
#ifdef THIS_IS_NOT_YET_IMPLEMENTED
  case kMinosNearDet:
    SetBeamFluxWindow(103648.837);
    break;
  case kMinosFarDear:
    SetBeamFluxWindow(73534000.);
    break;
  case kMinosNearRock:
    SetBeamFluxWindow();
    break;
  case kMinosFarRock:
    SetBeamFluxWindow();
    break;
#endif
  case kMinosNearCenter:
    {
      SetLengthUnits(genie::utils::units::UnitFromString("m"));
      fFluxWindowBase = kPosCenterNearBeam;
      fFluxWindowDir1 = TLorentzVector();  // no extent
      fFluxWindowDir2 = TLorentzVector();
      TLorentzVector usrpos;
      Beam2UserPos(fFluxWindowBase, usrpos);
      fFluxWindowPtUser[0] = usrpos.Vect();
      fFluxWindowPtUser[1] = fFluxWindowPtUser[0];
      fFluxWindowPtUser[2] = fFluxWindowPtUser[0];
      fFluxWindowLen1 = 0;
      fFluxWindowLen2 = 0;
      break;
    }
  case kMinosFarCenter:
    {
      SetLengthUnits(genie::utils::units::UnitFromString("m"));
      fFluxWindowBase = kPosCenterFarBeam;
      fFluxWindowDir1 = TLorentzVector();  // no extent
      fFluxWindowDir2 = TLorentzVector();
      TLorentzVector usrpos;
      Beam2UserPos(fFluxWindowBase, usrpos);
      fFluxWindowPtUser[0] = usrpos.Vect();
      fFluxWindowPtUser[1] = fFluxWindowPtUser[0];
      fFluxWindowPtUser[2] = fFluxWindowPtUser[0];
      fFluxWindowLen1 = 0;
      fFluxWindowLen2 = 0;
      break;
    }
  default:
    LOG("Flux", pERROR)
      << "SetBeamFluxWindow - StdFluxWindow " << stdwindow 
      << " not yet implemented";
    return false;
  }
  LOG("Flux",pNOTICE) << "about to CalcEffPOTsPerNu";
  this->CalcEffPOTsPerNu();
  return true;
}
//___________________________________________________________________________
void GNuMIFlux::SetFluxWindow(TVector3 p0, TVector3 p1, TVector3 p2)
                             // bool inDetCoord)  future extension
{
  // set flux window
  // NOTE: internally these are in "cm", but user might have set a preference
  fUseFluxAtDetCenter  = 0;
  fDetLocIsSet         = true;

  fFluxWindowPtUser[0] = p0;
  fFluxWindowPtUser[1] = p1;
  fFluxWindowPtUser[2] = p2;

  // convert from user to beam coord and from 3 points to base + 2 directions
  // apply units conversion
  TLorentzVector ptbm0, ptbm1, ptbm2;
  User2BeamPos(TLorentzVector(fFluxWindowPtUser[0],0),ptbm0);
  User2BeamPos(TLorentzVector(fFluxWindowPtUser[1],0),ptbm1);
  User2BeamPos(TLorentzVector(fFluxWindowPtUser[2],0),ptbm2);

  fFluxWindowBase = ptbm0;
  fFluxWindowDir1 = ptbm1 - ptbm0;
  fFluxWindowDir2 = ptbm2 - ptbm0;

  fFluxWindowLen1 = fFluxWindowDir1.Mag();
  fFluxWindowLen2 = fFluxWindowDir2.Mag();
  fWindowNormal = fFluxWindowDir1.Vect().Cross(fFluxWindowDir2.Vect()).Unit();

  double dot = fFluxWindowDir1.Dot(fFluxWindowDir2);
  if ( TMath::Abs(dot) > 1.0e-8 ) 
    LOG("Flux",pWARN) << "Dot product between window direction vectors was "
                      << dot << "; please check for orthoganality";
  
  LOG("Flux",pNOTICE) << "about to CalcEffPOTsPerNu";
  this->CalcEffPOTsPerNu();
}

//___________________________________________________________________________
void GNuMIFlux::GetFluxWindow(TVector3& p0, TVector3& p1, TVector3& p2) const
{
  // return flux window points
  p0 = fFluxWindowPtUser[0];
  p1 = fFluxWindowPtUser[1];
  p2 = fFluxWindowPtUser[2];
  
}
//___________________________________________________________________________
void GNuMIFlux::SetBeamRotation(TRotation beamrot)
{
  // rotation is really only 3-d vector, but we'll be operating on LorentzV's
  fBeamRot    = TLorentzRotation(beamrot);
  fBeamRotInv = fBeamRot.Inverse();
}

void GNuMIFlux::SetBeamCenter(TVector3 beam0)
{
  // set coord transform between detector and beam
  // NOTE: internally these are in "cm", but user might have set a preference
  fBeamZero = TLorentzVector(beam0,0);  // no time shift
}

//___________________________________________________________________________
TRotation GNuMIFlux::GetBeamRotation() const
{
  // rotation is really only 3-d vector, but we'll be operating on LorentzV's
  // give people back the original TRotation ... not pretty
  // ... it think this is right
  TRotation rot3;
  const TLorentzRotation& rot4 = fBeamRot;
  TVector3 newX(rot4.XX(),rot4.XY(),rot4.XZ());
  TVector3 newY(rot4.YX(),rot4.YY(),rot4.YZ());
  TVector3 newZ(rot4.ZX(),rot4.ZY(),rot4.ZZ());
  rot3.RotateAxes(newX,newY,newZ);
  return rot3.Inverse();
}
TVector3 GNuMIFlux::GetBeamCenter() const
{
  TVector3 beam0 = fBeamZero.Vect();
  return beam0;
}

//___________________________________________________________________________
//void GNuMIFlux::SetCoordTransform(TVector3 beam0, TRotation beamrot)
//{
//  // set coord transform between detector and beam
//  // NOTE: internally these are in "cm", but user might have set a preference
//
//  beam0 *= (1./fLengthScaleB2U);
//  fDetectorZero = TLorentzVector(beam0,0);  // no time shift
//  fDetectorRot  = TLorentzRotation(beamrot);
//
//}
//___________________________________________________________________________
//void GNuMIFlux::GetDetectorCoord(TLorentzVector& det0, TLorentzRotation& detrot) const
//{
//  // get coord transform between detector and beam
//  // NOTE: internally these are in "cm", but user might have set a preference
//
//  det0 = fDetectorZero;
//  det0 *= fLengthScaleB2U;
//  detrot = fDetectorRot;
//
//}
//___________________________________________________________________________

void GNuMIFlux::Beam2UserPos(const TLorentzVector& beamxyz, 
                                   TLorentzVector& usrxyz) const
{
  usrxyz = fLengthScaleB2U*(fBeamRot*beamxyz) + fBeamZero;
}
void GNuMIFlux::Beam2UserDir(const TLorentzVector& beamdir, 
                                   TLorentzVector& usrdir) const
{
  usrdir = fLengthScaleB2U*(fBeamRot*beamdir);
}
void GNuMIFlux::Beam2UserP4 (const TLorentzVector& beamp4, 
                                   TLorentzVector& usrp4 ) const
{
  usrp4 = fBeamRot*beamp4;
}

void GNuMIFlux::User2BeamPos(const TLorentzVector& usrxyz,
                                   TLorentzVector& beamxyz) const
{
  beamxyz = fLengthScaleU2B*(fBeamRotInv*(usrxyz-fBeamZero));
}
void GNuMIFlux::User2BeamDir(const TLorentzVector& usrdir,
                                   TLorentzVector& beamdir) const
{
  beamdir = fLengthScaleU2B*(fBeamRotInv*usrdir);
}
void GNuMIFlux::User2BeamP4 (const TLorentzVector& usrp4,
                                   TLorentzVector& beamp4) const
{
  beamp4 = fBeamRotInv*usrp4;
}

//___________________________________________________________________________
void GNuMIFlux::PrintCurrent(void)
{
  LOG("Flux", pNOTICE) << "CurrentEntry:" << *fCurEntry;
}
//___________________________________________________________________________
void GNuMIFlux::Clear(Option_t * opt)
{
  // Clear the driver state
  //
  LOG("Flux", pWARN) << "GSimpleNtpFlux::Clear(" << opt << ") called";
  // do it in all cases, but EVGDriver/GMCJDriver will pass "CycleHistory"

  fICycle     = 0;

  fSumWeight  = 0;
  fNNeutrinos = 0;
  fAccumPOTs  = 0;
}
//___________________________________________________________________________
void GNuMIFlux::GenerateWeighted(bool gen_weighted)
{
  // Set whether to generate weighted rays
  //
  fGenWeighted = gen_weighted;
}
//___________________________________________________________________________
void GNuMIFlux::Initialize(void)
{
  LOG("Flux", pNOTICE) << "Initializing GNuMIFlux driver";

  fMaxEv           =  0;
  fEnd             =  false;

  fCurEntry        = new GNuMIFluxPassThroughInfo;

  fNuFluxTree      =  0;
  fG3NuMI          =  0;
  fG4NuMI          =  0;
  fFlugg           =  0;
  fNuFluxTreeName  = "";
  fNuFluxGen       = "";
  fNFiles          =  0;

  fNEntries        =  0;
  fIEntry          = -1;
  fICycle          =  0;
  fNUse            =  1;
  fIUse            =  999999;

  fNuTot           = 0;
  fFilePOTs        = 0;

  fMaxWeight       = -1;
  fMaxWgtFudge     =  1.05;
  fMaxWgtEntries   = 2500000;
  fMaxEFudge       =  0;

  fSumWeight       =  0;
  fNNeutrinos      =  0;
  fEffPOTsPerNu    =  0;
  fAccumPOTs       =  0;

  fGenWeighted     = false;
  fApplyTiltWeight = true;
  fUseFluxAtDetCenter = 0;
  fDetLocIsSet        = false;
  // by default assume user length is cm
  SetLengthUnits(genie::utils::units::UnitFromString("m"));

  this->SetDefaults();
  this->ResetCurrent();
}
//___________________________________________________________________________
void GNuMIFlux::SetDefaults(void)
{
// - Set default neutrino species list (nue, nuebar, numu, numubar) and
//   maximum energy (120 GeV).
//   These defaults can be overwritten by user calls (at the driver init) to
//   GNuMIlux::SetMaxEnergy(double Ev) and
//   GNuMIFlux::SetFluxParticles(const PDGCodeList & particles)
// - Set the default file normalization to 1E+21 POT
// - Set the default flux neutrino start z position at -5m (z=0 is the
//   detector centre).
// - Set number of cycles to 1

  LOG("Flux", pNOTICE) << "Setting default GNuMIFlux driver options";

  PDGCodeList particles;
  particles.push_back(kPdgNuMu);
  particles.push_back(kPdgAntiNuMu);
  particles.push_back(kPdgNuE);
  particles.push_back(kPdgAntiNuE);

  this->SetFluxParticles (particles);
  this->SetMaxEnergy     (120./*GeV*/);  // was 200, but that would be wasteful
  this->SetUpstreamZ     (-3.4e38); // way upstream ==> use flux window
  this->SetNumOfCycles   (0);
  this->SetEntryReuse    (1);

  this->SetXMLFileBase("GNuMIFlux.xml");
}
//___________________________________________________________________________
void GNuMIFlux::ResetCurrent(void)
{
// reset running values of neutrino pdg-code, 4-position & 4-momentum
// and the input ntuple leaves

  fCurEntry->ResetCurrent();
  fCurEntry->ResetCopy();
}
//___________________________________________________________________________
void GNuMIFlux::CleanUp(void)
{
  LOG("Flux", pNOTICE) << "Cleaning up...";

  if (fPdgCList)    delete fPdgCList;
  if (fPdgCListRej) delete fPdgCListRej;
  if (fCurEntry)    delete fCurEntry;

  if ( fG3NuMI )    delete fG3NuMI;
  if ( fG4NuMI )    delete fG4NuMI;
  if ( fFlugg  )    delete fFlugg;

  LOG("Flux", pNOTICE)
    << " flux file cycles: " << fICycle << " of " << fNCycles 
    << ", entry " << fIEntry << " use: " << fIUse << " of " << fNUse;
}

//___________________________________________________________________________
void GNuMIFlux::AddFile(TTree* thetree, string fname)
{
  // Add a file to the chain

  ULong64_t nentries = thetree->GetEntries();

  // first/last "evtno" are the proton # of the first/last proton
  // that generated a neutrino ... not necessarily true first/last #
  // estimate we're probably not off by more than 100 ...
  // !Important change
  // Some files (due to some intermediate flugg issues) list evtno==0
  // when it isn't really true, we need to scan nearby values in case the
  // last entry is one of these (otherwise file contributes 0 POTs).  
  // Also assume quantization of 500 (instead of 100).

  thetree->SetMakeClass(1); // need full ntuple decomposition for
  // the SetBranchAddress to work on g4numi ntuples.  Works fine
  // without it on gnumi (geant3) and flugg ntuples [each with their
  // own slight differences] but shouldn't harm them either.

  Int_t evtno = 0;
  TBranch* br_evtno = 0;
  thetree->SetBranchAddress("evtno",&evtno, &br_evtno);
  Int_t evt_1 = 0x7FFFFFFF;
  Int_t evt_N = 1;
#define OLDGUESS
#ifdef OLDGUESS
  for (int j=0; j<50; ++j) {
    thetree->GetEntry(j);
    if (evtno != 0) evt_1 = TMath::Min(evtno,evt_1);
    thetree->GetEntry(nentries-1 -j );
    if (evtno != 0) evt_N = TMath::Max(evtno,evt_N);
  }
  // looks like low counts are due to "evtno" not including
  // protons that miss the actual target (hit baffle, etc)
  // this will vary with beam conditions parameters
  // so we should round way up, those generating flugg files
  // aren't going to quantize less than 1000
  // though 500 would probably be okay, 100 would not.
  const Int_t    nquant = 1000; // 500;  // 100
  const Double_t rquant = nquant;
#else
  for (int j=0; j<50; ++j) {
    thetree->GetEntry(j);
    if (evtno != 0) evt_1 = TMath::Min(evtno,evt_1);
    std::cout << "[" << j << "] evtno=" << evtno << " evt_1=" << evt_1 << std::endl;
  }
  for (int j=0; j<50; ++j) {
    thetree->GetEntry(nentries-1 -j );
    if (evtno != 0) evt_N = TMath::Max(evtno,evt_N);
    std::cout << "[" << (nentries-1-j) << "] evtno=" << evtno << " evt_N=" << evt_N << std::endl;
  }

  Int_t    nquant = 1000; // 500;
  Double_t rquant = nquant;
#endif

  Int_t est_1 = (TMath::FloorNint(evt_1/rquant))*nquant + 1;
  Int_t est_N = (TMath::FloorNint((evt_N-1)/rquant)+1)*nquant;
  ULong64_t npots = est_N - est_1 + 1;
  LOG("Flux",pNOTICE) //INFO)
    << fNuFluxTreeName << "->AddFile() of " << nentries << " entries ["
    << evt_1 << ":" << evt_N << "%" << nquant 
    << "(" <<  est_1 << ":" << est_N << ")=" 
    << npots <<" POTs] in {" << fNuFluxGen << "} file: " << fname;
  fNuTot    += nentries;
  fFilePOTs += npots;
  fNFiles++;

  fNuFluxTree->AddFile(fname.c_str());
}

//___________________________________________________________________________
void GNuMIFlux::SetLengthUnits(double user_units)
{
  // Set the scale factor for lengths going from beam (cm) to user coordinates

  // GNuMIFlux uses "cm" as the length unit consistently internally (this is 
  // the length units used by both the g3 and g4 ntuples).  User interactions 
  // setting the beam-to-detector coordinate transform, flux window, and the 
  // returned position might need to be in other units.  Use:
  //     double scale = genie::utils::units::UnitFromString("cm");
  // ( #include "Utils/UnitUtils.h for declaration )
  // to get the correct scale factor to pass in.

  double rescale = fLengthUnits / user_units;
  fLengthUnits = user_units;
  double cm = genie::utils::units::UnitFromString("cm");
  fLengthScaleB2U = cm / user_units;
  fLengthScaleU2B = user_units / cm;

  // in case we're changing units without resetting transform/window
  // not recommended, but should work
  fCurEntry->fgX4User  *= rescale;
  fBeamZero            *= rescale;
  fFluxWindowPtUser[0] *= rescale;
  fFluxWindowPtUser[1] *= rescale;
  fFluxWindowPtUser[2] *= rescale;

  // case GNuMIFlux::kmeter:  fLengthScaleB2U =   0.01  ; break;
  // case GNuMIFlux::kcm:     fLengthScaleB2U =   1.    ; break;
  // case GNuMIFlux::kmm:     fLengthScaleB2U =  10.    ; break;
  // case GNuMIFlux::kfm:     fLengthScaleB2U =   1.e13 ; break;  // 10e-2m -> 10e-15m

}

//___________________________________________________________________________
double GNuMIFlux::LengthUnits(void) const
{
  // Return the scale factor for lengths the user is getting
  double cm = genie::utils::units::UnitFromString("cm");
  return fLengthScaleB2U * cm ;
}

//___________________________________________________________________________
GNuMIFluxPassThroughInfo::GNuMIFluxPassThroughInfo()
  : TObject()
{ 
  ResetCopy(); 
  ResetCurrent();
}

//___________________________________________________________________________
void GNuMIFluxPassThroughInfo::ResetCopy()
{
  pcodes   = -1;
  units    = -1;
  
  run      = -1;
  evtno    = -1;
  ndxdz    = 0.;
  ndydz    = 0.;
  npz      = 0.;
  nenergy  = 0.;
  ndxdznea = 0.;
  ndydznea = 0.;
  nenergyn = 0.;
  nwtnear  = 0.;
  ndxdzfar = 0.;
  ndydzfar = 0.;
  nenergyf = 0.;
  nwtfar   = 0.;
  norig    = 0;
  ndecay   = 0;
  ntype    = 0;
  vx       = 0.;
  vy       = 0.;
  vz       = 0.;
  pdpx     = 0.;
  pdpy     = 0.;
  pdpz     = 0.;
  ppdxdz   = 0.;
  ppdydz   = 0.;
  pppz     = 0.;
  ppenergy = 0.;
  ppmedium = 0 ;
  ptype    = 0 ;
  ppvx     = 0.;
  ppvy     = 0.;
  ppvz     = 0.;
  muparpx  = 0.;
  muparpy  = 0.;
  muparpz  = 0.;
  mupare   = 0.;
  necm     = 0.;
  nimpwt   = 0.;
  xpoint   = 0.;
  ypoint   = 0.;
  zpoint   = 0.;
  tvx      = 0.;
  tvy      = 0.;
  tvz      = 0.;
  tpx      = 0.;
  tpy      = 0.;
  tpz      = 0.;
  tptype   = 0 ;
  tgen     = 0 ;
  tgptype  = 0 ;
  tgppx    = 0.;
  tgppy    = 0.;
  tgppz    = 0.;
  tprivx   = 0.;
  tprivy   = 0.;
  tprivz   = 0.;
  beamx    = 0.;
  beamy    = 0.;
  beamz    = 0.;
  beampx   = 0.;
  beampy   = 0.;
  beampz   = 0.;

#ifndef SKIP_MINERVA_MODS
  //=========================================
  // The following was inserted by MINERvA
  //=========================================
  ntrajectory = 0;
  overflow    = false;

  for ( unsigned int itclear = 0; itclear < MAX_N_TRAJ; itclear++ ) {
     pdgcode[itclear]  = 0;
     trackId[itclear]  = 0;
     parentId[itclear] = 0;
     proc[itclear]     = 0;
     ivol[itclear]     = 0;
     fvol[itclear]     = 0;
     startx[itclear]   = 0;
     starty[itclear]   = 0;
     startz[itclear]   = 0;
     startpx[itclear]  = 0;
     startpy[itclear]  = 0;
     startpz[itclear]  = 0;
     stopx[itclear]    = 0;
     stopy[itclear]    = 0;
     stopz[itclear]    = 0;
     stoppx[itclear]   = 0;
     stoppy[itclear]   = 0;
     stoppz[itclear]   = 0;
     pprodpx[itclear]  = 0;
     pprodpy[itclear]  = 0;
     pprodpz[itclear]  = 0;
  }
  //END of minerva additions
#endif
}

//___________________________________________________________________________
void GNuMIFluxPassThroughInfo::ResetCurrent()
{
  // reset the state of the "generated" entry
  fgPdgC  = 0;
  fgXYWgt = 0;
  fgP4.SetPxPyPzE(0.,0.,0.,0.);
  fgX4.SetXYZT(0.,0.,0.,0.);
}

//___________________________________________________________________________
/*******

ALLOW DEFAULT COPY CONSTRUCTOR
GNuMIFluxPassThroughInfo::GNuMIFluxPassThroughInfo(
                        const GNuMIFluxPassThroughInfo & info) :
TObject(),
pcodes   ( info.pcodes   ),
units    ( info.units    ),
run      ( info.run      ),
evtno    ( info.evtno    ),
ndxdz    ( info.ndxdz    ),
ndydz    ( info.ndydz    ),
npz      ( info.npz      ),
nenergy  ( info.nenergy  ),
ndxdznea ( info.ndxdznea ),
ndydznea ( info.ndydznea ),
nenergyn ( info.nenergyn ),
nwtnear  ( info.nwtnear  ),
ndxdzfar ( info.ndxdzfar ),
ndydzfar ( info.ndydzfar ),
nenergyf ( info.nenergyf ),
nwtfar   ( info.nwtfar   ),
norig    ( info.norig    ),
ndecay   ( info.ndecay   ),
ntype    ( info.ntype    ),
vx       ( info.vx       ),
vy       ( info.vy       ),
vz       ( info.vz       ),
pdPx     ( info.pdPx     ),
pdPy     ( info.pdPy     ),
pdPz     ( info.pdPz     ),
ppdxdz   ( info.ppdxdz   ),
ppdydz   ( info.ppdydz   ),
pppz     ( info.pppz     ),
ppenergy ( info.ppenergy ),
ppmedium ( info.ppmedium ),
ptype    ( info.ptype    ),
ppvx     ( info.ppvx     ),
ppvy     ( info.ppvy     ),
ppvz     ( info.ppvz     ),
muparpx  ( info.muparpx  ),
muparpy  ( info.muparpy  ),
muparpz  ( info.muparpz  ),
mupare   ( info.mupare   ),
necm     ( info.necm     ),
nimpwt   ( info.nimpwt   ),
xpoint   ( info.xpoint   ),
ypoint   ( info.ypoint   ),
zpoint   ( info.zpoint   ),
tvx      ( info.tvx      ),
tvy      ( info.tvy      ),
tvz      ( info.tvz      ),
tpx      ( info.tpx      ),
tpy      ( info.tpy      ),
tpz      ( info.tpz      ),
tptype   ( info.tptype   ),
tgen     ( info.tgen     ),
tgptype  ( info.tgptype  ),
tgppx    ( info.tgppx    ),
tgppy    ( info.tgppy    ),
tgppz    ( info.tgppz    ),
tprivx   ( info.tprivx   ),
tprivy   ( info.tprivy   ),
tprivz   ( info.tprivz   ),
beamx    ( info.beamx    ),
beamy    ( info.beamy    ),
beamz    ( info.beamz    ),
beampx   ( info.beampx   ),
beampy   ( info.beampy   ),
beampz   ( info.beampz   )
{

}
*************/

//___________________________________________________________________________
void GNuMIFluxPassThroughInfo::ConvertPartCodes()
{
  if ( pcodes == 0 ) {
    pcodes = 1;  // flag that conversion has been made
    switch ( ntype ) {
    case 56: ntype = kPdgNuMu;     break;
    case 55: ntype = kPdgAntiNuMu; break;
    case 53: ntype = kPdgNuE;      break;
    case 52: ntype = kPdgAntiNuE;  break;
    default:
      LOG("Flux", pNOTICE)
        << "ConvertPartCodes saw ntype " << ntype << " -- unknown ";
    }
    if ( ptype   != 0 ) ptype   = pdg::GeantToPdg(ptype);
    if ( tptype  != 0 ) tptype  = pdg::GeantToPdg(tptype);
    if ( tgptype != 0 ) tgptype = pdg::GeantToPdg(tgptype);
  } else if ( pcodes != 1 ) {
    // flag as unknown state ...
    LOG("Flux", pNOTICE)
      << "ConvertPartCodes saw pcodes flag as " << pcodes;
  }

}

//___________________________________________________________________________
void GNuMIFluxPassThroughInfo::Print(const Option_t* /* opt */ ) const
{
  std::cout << *this << std::endl;
}

//___________________________________________________________________________
void GNuMIFluxPassThroughInfo::MakeCopy(const g3numi* g3 )
{
  run      = g3->run;
  evtno    = g3->evtno;
  ndxdz    = g3->Ndxdz;
  ndydz    = g3->Ndydz;
  npz      = g3->Npz;
  nenergy  = g3->Nenergy;
  ndxdznea = g3->Ndxdznea;
  ndydznea = g3->Ndydznea;
  nenergyn = g3->Nenergyn;
  nwtnear  = g3->Nwtnear;
  ndxdzfar = g3->Ndxdzfar;
  ndydzfar = g3->Ndydzfar;
  nenergyf = g3->Nenergyf;
  nwtfar   = g3->Nwtfar;
  norig    = g3->Norig;
  ndecay   = g3->Ndecay;
  ntype    = g3->Ntype;
  vx       = g3->Vx;
  vy       = g3->Vy;
  vz       = g3->Vz;
  pdpx     = g3->pdpx;
  pdpy     = g3->pdpy;
  pdpz     = g3->pdpz;
  ppdxdz   = g3->ppdxdz;
  ppdydz   = g3->ppdydz;
  pppz     = g3->pppz;
  ppenergy = g3->ppenergy;
  ppmedium = g3->ppmedium;
  ptype    = g3->ptype;
  ppvx     = g3->ppvx;
  ppvy     = g3->ppvy;
  ppvz     = g3->ppvz;
  muparpx  = g3->muparpx;
  muparpy  = g3->muparpy;
  muparpz  = g3->muparpz;
  mupare   = g3->mupare;

  necm     = g3->Necm;
  nimpwt   = g3->Nimpwt;
  xpoint   = g3->xpoint;
  ypoint   = g3->ypoint;
  zpoint   = g3->zpoint;

  tvx      = g3->tvx;
  tvy      = g3->tvy;
  tvz      = g3->tvz;
  tpx      = g3->tpx;
  tpy      = g3->tpy;
  tpz      = g3->tpz;
  tptype   = g3->tptype;
  tgen     = g3->tgen;
  tgptype  = g3->tgptype;
  tgppx    = g3->tgppx;
  tgppy    = g3->tgppy;
  tgppz    = g3->tgppz;
  tprivx   = g3->tprivx;
  tprivy   = g3->tprivy;
  tprivz   = g3->tprivz;
  beamx    = g3->beamx;
  beamy    = g3->beamy;
  beamz    = g3->beamz;
  beampx   = g3->beampx;
  beampy   = g3->beampy;
  beampz   = g3->beampz;

}

//___________________________________________________________________________
void GNuMIFluxPassThroughInfo::MakeCopy(const g4numi* g4 )
{

  const int kNearIndx = 0;
  const int kFarIndx  = 0;

  run      = g4->run;
  evtno    = g4->evtno;
  ndxdz    = g4->Ndxdz;
  ndydz    = g4->Ndydz;
  npz      = g4->Npz;
  nenergy  = g4->Nenergy;
  ndxdznea = g4->NdxdzNear[kNearIndx];
  ndydznea = g4->NdydzNear[kNearIndx];
  nenergyn = g4->NenergyN[kNearIndx];
  nwtnear  = g4->NWtNear[kNearIndx];
  ndxdzfar = g4->NdxdzFar[kFarIndx];
  ndydzfar = g4->NdydzFar[kFarIndx];
  nenergyf = g4->NenergyF[kFarIndx];
  nwtfar   = g4->NWtFar[kFarIndx];
  norig    = g4->Norig;
  ndecay   = g4->Ndecay;
  ntype    = g4->Ntype;
  vx       = g4->Vx;
  vy       = g4->Vy;
  vz       = g4->Vz;
  pdpx     = g4->pdPx;
  pdpy     = g4->pdPy;
  pdpz     = g4->pdPz;
  ppdxdz   = g4->ppdxdz;
  ppdydz   = g4->ppdydz;
  pppz     = g4->pppz;
  ppenergy = g4->ppenergy;
  ppmedium = (Int_t)g4->ppmedium;  // int in g3, double in g4!
  ptype    = g4->ptype;
  ppvx     = g4->ppvx;
  ppvy     = g4->ppvy;
  ppvz     = g4->ppvz;
  muparpx  = g4->muparpx;
  muparpy  = g4->muparpy;
  muparpz  = g4->muparpz;
  mupare   = g4->mupare;

  necm     = g4->Necm;
  nimpwt   = g4->Nimpwt;
  xpoint   = g4->xpoint;
  ypoint   = g4->ypoint;
  zpoint   = g4->zpoint;

  tvx      = g4->tvx;
  tvy      = g4->tvy;
  tvz      = g4->tvz;
  tpx      = g4->tpx;
  tpy      = g4->tpy;
  tpz      = g4->tpz;
  tptype   = g4->tptype;
  tgen     = g4->tgen;
  tgptype  = 0 ;  // no equivalent in g4
  tgppx    = 0.;
  tgppy    = 0.;
  tgppz    = 0.;
  tprivx   = 0.;
  tprivy   = 0.;
  tprivz   = 0.;
  beamx    = 0.; // g4->protonX;
  beamy    = 0.; // g4->protonY;
  beamz    = 0.; // g4->protonZ;
  beampx   = 0.; // g4->protonPx;
  beampy   = 0.; // g4->protonPy;
  beampz   = 0.; // g4->protonPz;

#ifndef SKIP_MINERVA_MODS
  //=========================================
  // The following was inserted by MINERvA
  //=========================================
  ntrajectory = g4->ntrajectory;
  overflow    = g4->overflow;
  // Limit number of trajectories to MAX_N_TRAJ
  Int_t ntraj = ntrajectory;
  if ( overflow )
    ntraj = MAX_N_TRAJ;

  for ( Int_t ic = 0; ic < ntraj; ++ic ) {
    pdgcode[ic]  = g4->pdg[ic];
    trackId[ic]  = g4->trackId[ic];
    parentId[ic] = g4->parentId[ic];

    startx[ic]   = g4->startx[ic];
    starty[ic]   = g4->starty[ic];
    startz[ic]   = g4->startz[ic];
    stopx[ic]    = g4->stopx[ic];
    stopy[ic]    = g4->stopy[ic];
    stopz[ic]    = g4->stopz[ic];
    startpx[ic]  = g4->startpx[ic];
    startpy[ic]  = g4->startpy[ic];
    startpz[ic]  = g4->startpz[ic];
    stoppx[ic]   = g4->stoppx[ic];
    stoppy[ic]   = g4->stoppy[ic];
    stoppz[ic]   = g4->stoppz[ic];
    pprodpx[ic]  = g4->pprodpx[ic];
    pprodpy[ic]  = g4->pprodpy[ic];
    pprodpz[ic]  = g4->pprodpz[ic];

    proc[ic] =  getProcessID(g4->proc[ic]);
    ivol[ic] =  getVolID(g4->ivol[ic]);
    fvol[ic] =  getVolID(g4->fvol[ic]);
  }
  //END of minerva additions
#endif

}

//___________________________________________________________________________
void GNuMIFluxPassThroughInfo::MakeCopy(const flugg* f )
{
  run      = f->run;
  evtno    = f->evtno;
  ndxdz    = f->Ndxdz;
  ndydz    = f->Ndydz;
  npz      = f->Npz;
  nenergy  = f->Nenergy;
  ndxdznea = f->Ndxdznea;
  ndydznea = f->Ndydznea;
  nenergyn = f->Nenergyn;
  nwtnear  = f->Nwtnear;
  ndxdzfar = f->Ndxdzfar;
  ndydzfar = f->Ndydzfar;
  nenergyf = f->Nenergyf;
  nwtfar   = f->Nwtfar;
  norig    = f->Norig;
  ndecay   = f->Ndecay;
  ntype    = f->Ntype;
  vx       = f->Vx;
  vy       = f->Vy;
  vz       = f->Vz;
  pdpx     = f->pdPx;
  pdpy     = f->pdPy;
  pdpz     = f->pdPz;
  ppdxdz   = f->ppdxdz;
  ppdydz   = f->ppdydz;
  pppz     = f->pppz;
  ppenergy = f->ppenergy;
  ppmedium = f->ppmedium;
  ptype    = f->ptype;
  ppvx     = f->ppvx;
  ppvy     = f->ppvy;
  ppvz     = f->ppvz;
  muparpx  = f->muparpx;
  muparpy  = f->muparpy;
  muparpz  = f->muparpz;
  mupare   = f->mupare;

  necm     = f->Necm;
  nimpwt   = f->Nimpwt;
  xpoint   = f->xpoint;
  ypoint   = f->ypoint;
  zpoint   = f->zpoint;

  tvx      = f->tvx;
  tvy      = f->tvy;
  tvz      = f->tvz;
  tpx      = f->tpx;
  tpy      = f->tpy;
  tpz      = f->tpz;
  tptype   = f->tptype;
  tgen     = f->tgen;
  tgptype  = f->tgptype;
  tgppx    = f->tgppx;
  tgppy    = f->tgppy;
  tgppz    = f->tgppz;
  tprivx   = f->tprivx;
  tprivy   = f->tprivy;
  tprivz   = f->tprivz;
  beamx    = f->beamx;
  beamy    = f->beamy;
  beamz    = f->beamz;
  beampx   = f->beampx;
  beampy   = f->beampy;
  beampz   = f->beampz;

}

//___________________________________________________________________________
int GNuMIFluxPassThroughInfo::CalcEnuWgt(const TLorentzVector& xyz,
                                         double& enu, double& wgt_xy) const
{

  // Neutrino Energy and Weigth at arbitrary point
  // based on:
  //   NuMI-NOTE-BEAM-0109 (MINOS DocDB # 109)
  //   Title:   Neutrino Beam Simulation using PAW with Weighted Monte Carlos
  //   Author:  Rick Milburn
  //   Date:    1995-10-01

  // history:
  // jzh  3/21/96 grab R.H.Milburn's weighing routine
  // jzh  5/ 9/96 substantially modify the weighting function use dot product 
  //              instead of rotation vecs to get theta get all info except 
  //              det from ADAMO banks neutrino parent is in Particle.inc
  //              Add weighting factor for polarized muon decay
  // jzh  4/17/97 convert more code to double precision because of problems 
  //              with Enu>30 GeV
  // rwh 10/ 9/08 transliterate function from f77 to C++

  // original function description:
  //   Real function for use with PAW Ntuple To transform from destination
  //   detector geometry to the unit sphere moving with decaying hadron with
  //   velocity v, BETA=v/c, etc..  For (pseudo)scalar hadrons the decays will
  //   be isotropic in this  sphere so the fractional area (out of 4-pi) is the
  //   fraction of decays that hit the target.  For a given target point and 
  //   area, and given x-y components of decay transverse location and slope,
  //   and given decay distance from target ans given decay GAMMA and 
  //   rest-frame neutrino energy, the lab energy at the target and the 
  //   fractional solid angle in the rest-frame are determined.
  //   For muon decays, correction for non-isotropic nature of decay is done.

  // Arguments:
  //    double x, y, z :: position to evaluate for (enu,wgt_xy) 
  //        in *beam* frame coordinates  (cm units)
  //    double enu, wgt_xy :: resulting energy and weight
  // Return:
  //    int :: error code
  // Assumptions:
  //    Energies given in GeV
  //    Particle codes have been translated from GEANT into PDG codes

  // for now ... these _should_ come from DB
  // but use these hard-coded values to "exactly" reproduce old code
  //
  const double kPIMASS = 0.13957;
  const double kKMASS  = 0.49368;
  const double kK0MASS = 0.49767;
  const double kMUMASS = 0.105658389;
  const double kOMEGAMASS = 1.67245;

  const int kpdg_nue       =   12;  // extended Geant 53
  const int kpdg_nuebar    =  -12;  // extended Geant 52
  const int kpdg_numu      =   14;  // extended Geant 56
  const int kpdg_numubar   =  -14;  // extended Geant 55

  const int kpdg_muplus     =   -13;  // Geant  5
  const int kpdg_muminus    =    13;  // Geant  6
  const int kpdg_pionplus   =   211;  // Geant  8
  const int kpdg_pionminus  =  -211;  // Geant  9
  const int kpdg_k0long     =   130;  // Geant 10  ( K0=311, K0S=310 )
  const int kpdg_k0short    =   310;  // Geant 16
  const int kpdg_k0mix      =   311;  
  const int kpdg_kaonplus   =   321;  // Geant 11
  const int kpdg_kaonminus  =  -321;  // Geant 12
  const int kpdg_omegaminus =  3334;  // Geant 24
  const int kpdg_omegaplus  = -3334;  // Geant 32

  const double kRDET = 100.0;   // set to flux per 100 cm radius

  double xpos = xyz.X();
  double ypos = xyz.Y();
  double zpos = xyz.Z();

  enu    = 0.0;  // don't know what the final value is
  wgt_xy = 0.0;  // but set these in case we return early due to error


  // in principle we should get these from the particle DB
  // but for consistency testing use the hardcoded values
  double parent_mass = kPIMASS;
  switch ( this->ptype ) {
  case kpdg_pionplus:
  case kpdg_pionminus:
    parent_mass = kPIMASS;
    break;
  case kpdg_kaonplus:
  case kpdg_kaonminus:
    parent_mass = kKMASS;
    break;
  case kpdg_k0long:
  case kpdg_k0short:
  case kpdg_k0mix:
    parent_mass = kK0MASS;
    break;
  case kpdg_muplus:
  case kpdg_muminus:
    parent_mass = kMUMASS;
    break;
  case kpdg_omegaminus:
  case kpdg_omegaplus:
    parent_mass = kOMEGAMASS;
    break;
  default:
    std::cerr << "NU_REWGT unknown particle type " << this->ptype
              << std::endl << std::flush;
    LOG("Flux",pFATAL) << "NU_REWGT unknown particle type " << this->ptype;
    assert(0);
    return 1;
  }

  double parentp2 = ( this->pdpx*this->pdpx +
                      this->pdpy*this->pdpy +
                      this->pdpz*this->pdpz );
  double parent_energy = TMath::Sqrt( parentp2 +
                                     parent_mass*parent_mass);
  double parentp = TMath::Sqrt( parentp2 );

  double gamma     = parent_energy / parent_mass;
  double gamma_sqr = gamma * gamma;
  double beta_mag  = TMath::Sqrt( ( gamma_sqr - 1.0 )/gamma_sqr );

  // Get the neutrino energy in the parent decay CM
  double enuzr = this->necm;
  // Get angle from parent line of flight to chosen point in beam frame
  double rad = TMath::Sqrt( (xpos-this->vx)*(xpos-this->vx) +
                            (ypos-this->vy)*(ypos-this->vy) +
                            (zpos-this->vz)*(zpos-this->vz) );

  double emrat = 1.0;
  double costh_pardet = -999., theta_pardet = -999.;

  // boost correction, but only if parent hasn't stopped
  if ( parentp > 0. ) {
    costh_pardet = ( this->pdpx*(xpos-this->vx) +
                     this->pdpy*(ypos-this->vy) +
                     this->pdpz*(zpos-this->vz) ) 
                     / ( parentp * rad);
    if ( costh_pardet >  1.0 ) costh_pardet =  1.0;
    if ( costh_pardet < -1.0 ) costh_pardet = -1.0;
    theta_pardet = TMath::ACos(costh_pardet);

    // Weighted neutrino energy in beam, approx, good for small theta
    emrat = 1.0 / ( gamma * ( 1.0 - beta_mag * costh_pardet ));
  }

  enu = emrat * enuzr;  // the energy ... normally

  // RWH-debug
  bool debug = false;
  if (debug) {
    std::cout << std::setprecision(15);
    std::cout << "xyz (" << xpos << "," << ypos << "," << zpos << ")" << std::endl;
    std::cout << "ptype " << this->ptype << " m " << parent_mass 
              << " p " << parentp << " e " << parent_energy << " gamma " << gamma
              << " beta " << beta_mag << std::endl;

    std::cout << " enuzr " << enuzr << " rad " << rad << " costh " << costh_pardet
              << " theta " << theta_pardet << " emrat " << emrat 
              << " enu " << enu << std::endl;
  }

#ifdef  GNUMI_TEST_XY_WGT
  gpartials.xdet          = xpos;
  gpartials.ydet          = ypos;
  gpartials.zdet          = zpos;
  gpartials.parent_mass   = parent_mass;
  gpartials.parentp       = parentp;
  gpartials.parent_energy = parent_energy;
  gpartials.gamma         = gamma;
  gpartials.beta_mag      = beta_mag;
  gpartials.enuzr         = enuzr;
  gpartials.rad           = rad;
  gpartials.costh_pardet  = costh_pardet;
  gpartials.theta_pardet  = theta_pardet;
  gpartials.emrat         = emrat;
  gpartials.eneu          = enu;
#endif

  // Get solid angle/4pi for detector element
  // small angle approximation, fixed by Alex Radovic
  //SAA//double sangdet = ( kRDET*kRDET / ( (zpos-this->vz)*(zpos-this->vz)))/4.0;
  
  double sanddetcomp = TMath::Sqrt(( (xpos-this->vx)*(xpos-this->vx) ) +
                                   ( (ypos-this->vy)*(ypos-this->vy) ) +
                                   ( (zpos-this->vz)*(zpos-this->vz) )   );
  double sangdet = ( 1.0 - TMath::Cos(TMath::ATan( kRDET / sanddetcomp)))/2.0;

  // Weight for solid angle and lorentz boost
  wgt_xy = sangdet * ( emrat * emrat );  // ! the weight ... normally

#ifdef  GNUMI_TEST_XY_WGT
  gpartials.sangdet       = sangdet;
  gpartials.wgt           = wgt_xy;
  gpartials.ptype         = this->ptype; // assume already PDG
#endif

  // Done for all except polarized muon decay
  // in which case need to modify weight 
  // (must be done in double precision)
  if ( this->ptype == kpdg_muplus || this->ptype == kpdg_muminus) {
    double beta[3], p_dcm_nu[4], p_nu[3], p_pcm_mp[3], partial;

    // Boost neu neutrino to mu decay CM
    beta[0] = this->pdpx / parent_energy;
    beta[1] = this->pdpy / parent_energy;
    beta[2] = this->pdpz / parent_energy;
    p_nu[0] = (xpos-this->vx)*enu/rad;
    p_nu[1] = (ypos-this->vy)*enu/rad;
    p_nu[2] = (zpos-this->vz)*enu/rad;
    partial = gamma * 
      (beta[0]*p_nu[0] + beta[1]*p_nu[1] + beta[2]*p_nu[2] );
    partial = enu - partial/(gamma+1.0);
    // the following calculation is numerically imprecise
    // especially p_dcm_nu[2] leads to taking the difference of numbers of order ~10's
    // and getting results of order ~0.02's
    // for g3numi we're starting with floats (ie. good to ~1 part in 10^7)
    p_dcm_nu[0] = p_nu[0] - beta[0]*gamma*partial;
    p_dcm_nu[1] = p_nu[1] - beta[1]*gamma*partial;
    p_dcm_nu[2] = p_nu[2] - beta[2]*gamma*partial;
    p_dcm_nu[3] = TMath::Sqrt( p_dcm_nu[0]*p_dcm_nu[0] +
                               p_dcm_nu[1]*p_dcm_nu[1] +
                               p_dcm_nu[2]*p_dcm_nu[2] );

#ifdef  GNUMI_TEST_XY_WGT
    gpartials.betanu[0]     = beta[0];
    gpartials.betanu[1]     = beta[1];
    gpartials.betanu[2]     = beta[2];
    gpartials.p_nu[0]       = p_nu[0];
    gpartials.p_nu[1]       = p_nu[1];
    gpartials.p_nu[2]       = p_nu[2];
    gpartials.partial1      = partial;
    gpartials.p_dcm_nu[0]   = p_dcm_nu[0];
    gpartials.p_dcm_nu[1]   = p_dcm_nu[1];
    gpartials.p_dcm_nu[2]   = p_dcm_nu[2];
    gpartials.p_dcm_nu[3]   = p_dcm_nu[3];
#endif

    // Boost parent of mu to mu production CM
    double particle_energy = this->ppenergy;
    gamma = particle_energy/parent_mass;
    beta[0] = this->ppdxdz * this->pppz / particle_energy;
    beta[1] = this->ppdydz * this->pppz / particle_energy;
    beta[2] =                    this->pppz / particle_energy;
    partial = gamma * ( beta[0]*this->muparpx + 
                        beta[1]*this->muparpy + 
                        beta[2]*this->muparpz );
    partial = this->mupare - partial/(gamma+1.0);
    p_pcm_mp[0] = this->muparpx - beta[0]*gamma*partial;
    p_pcm_mp[1] = this->muparpy - beta[1]*gamma*partial;
    p_pcm_mp[2] = this->muparpz - beta[2]*gamma*partial;
    double p_pcm = TMath::Sqrt ( p_pcm_mp[0]*p_pcm_mp[0] +
                                 p_pcm_mp[1]*p_pcm_mp[1] +
                                 p_pcm_mp[2]*p_pcm_mp[2] );

    //std::cout << " muparpxyz " << this->muparpx << " "
    //          << this->muparpy << " " << this->muparpz << std::endl;
    //std::cout << " beta " << beta[0] << " " << beta[1] << " " << beta[2] << std::endl;
    //std::cout << " gamma " << gamma << " partial " << partial << std::endl;
    //std::cout << " p_pcm_mp " << p_pcm_mp[0] << " " << p_pcm_mp[1] << " " 
    //          << p_pcm_mp[2] << " " << p_pcm << std::endl;

#ifdef  GNUMI_TEST_XY_WGT
    gpartials.muparent_px   = this->muparpx;
    gpartials.muparent_py   = this->muparpy;
    gpartials.muparent_pz   = this->muparpz;
    gpartials.gammamp       = gamma;
    gpartials.betamp[0]     = beta[0];
    gpartials.betamp[1]     = beta[1];
    gpartials.betamp[2]     = beta[2];
    gpartials.partial2      = partial;
    gpartials.p_pcm_mp[0]   = p_pcm_mp[0];
    gpartials.p_pcm_mp[1]   = p_pcm_mp[1];
    gpartials.p_pcm_mp[2]   = p_pcm_mp[2];
    gpartials.p_pcm         = p_pcm;
#endif

    const double eps = 1.0e-30;  // ? what value to use
    if ( p_pcm < eps || p_dcm_nu[3] < eps ) {
      return 3; // mu missing parent info?
    }
    // Calc new decay angle w.r.t. (anti)spin direction
    double costh = ( p_dcm_nu[0]*p_pcm_mp[0] +
                     p_dcm_nu[1]*p_pcm_mp[1] +
                     p_dcm_nu[2]*p_pcm_mp[2] ) /
                   ( p_dcm_nu[3]*p_pcm );
    if ( costh >  1.0 ) costh =  1.0;
    if ( costh < -1.0 ) costh = -1.0;
    // Calc relative weight due to angle difference
    double wgt_ratio = 0.0;
    switch ( this->ntype ) {
    case kpdg_nue:
    case kpdg_nuebar:
      wgt_ratio = 1.0 - costh;
      break;
    case kpdg_numu:
    case kpdg_numubar:
    {
      double xnu = 2.0 * enuzr / kMUMASS;
      wgt_ratio = ( (3.0-2.0*xnu )  - (1.0-2.0*xnu)*costh ) / (3.0-2.0*xnu);
      break;
    }
    default:
      return 2; // bad neutrino type
    }
    wgt_xy = wgt_xy * wgt_ratio;

#ifdef  GNUMI_TEST_XY_WGT
    gpartials.ntype     = this->ntype; // assume converted to PDG
    gpartials.costhmu   = costh;
    gpartials.wgt_ratio = wgt_ratio;
#endif

  } // ptype is muon

  return 0;
}

//___________________________________________________________________________


namespace genie {
namespace flux  {
  ostream & operator << (
    ostream & stream, const genie::flux::GNuMIFluxPassThroughInfo & info)
    {
      // stream << "\n ndecay   = " << info.ndecay << std::endl;
      stream << "\nGNuMIFlux run " << info.run << " evtno " << info.evtno
             << " (pcodes " << info.pcodes << " units " << info.units << ")"
             << "\n random dk: dx/dz " << info.ndxdz 
             << " dy/dz " << info.ndydz
             << " pz " <<  info.npz << " E " << info.nenergy
             << "\n near00 dk: dx/dz " << info.ndxdznea 
             << " dy/dz " << info.ndydznea
             << " E " <<  info.nenergyn << " wgt " << info.nwtnear
             << "\n far00  dk: dx/dz " << info.ndxdzfar
             << " dy/dz " << info.ndydzfar
             << " E " <<  info.nenergyf << " wgt " << info.nwtfar
             << "\n norig " << info.norig << " ndecay " << info.ndecay
             << " ntype " << info.ntype
             << "\n had vtx " << info.vx << " " << info.vy << " " << info.vz
             << "\n parent p3 @ dk " << info.pdpx << " " << info.pdpy << " " << info.pdpz
             << "\n parent prod: dx/dz " << info.ppdxdz 
             << " dy/dz " << info.ppdydz
             << " pz " << info.pppz << " E " << info.ppenergy
             << "\n ppmedium " << info.ppmedium << " ptype " << info.ptype
             << " ppvtx " << info.ppvx << " " << info.ppvy << " " << info.ppvz
             << "\n mu parent p4 " << info.muparpx << " " << info.muparpy
             << " " << info.muparpz << " " << info.mupare
             << "\n necm " << info.necm << " nimpwt " << info.nimpwt
             << "\n point x,y,z " << info.xpoint << " " << info.ypoint
             << " " << info.zpoint
             << "\n tv x,y,z " << info.tvx << " " << info.tvy << " " << info.tvz
             << "\n tptype " << info.tptype << " tgen " << info.tgen
             << " tgptype " << info.tgptype 
             << "\n tgp px,py,pz " << info.tgppx << " " << info.tgppy
             << " " << info.tgppz
             << "\n tpriv x,y,z " << info.tprivx << " " << info.tprivy
             << " " << info.tprivz
             << "\n beam x,y,z " << info.beamx << " " << info.beamy
             << " " << info.beamz
             << "\n beam px,py,pz " << info.beampx << " " << info.beampy
             << " " << info.beampz
        ;

#ifndef SKIP_MINERVA_MODS
      //=========================================
      // The following was inserted by MINERvA
      //=========================================
      stream << "\nNeutrino History : ntrajectories " << info.ntrajectory
             << "\n (trkID, pdg) of nu parents: [";
      Int_t ntraj = info.ntrajectory;
      if ( info.overflow ) ntraj = GNuMIFluxPassThroughInfo::MAX_N_TRAJ;

      for ( Int_t itt = 0; itt < ntraj; ++itt )
        stream << "(" << info.trackId[itt-1] << ", " << info.pdgcode[itt] << ") ";
      stream << "]\n";
      //END of minerva additions
#endif

      stream << "\nCurrent: pdg " << info.fgPdgC 
             << " xywgt " << info.fgXYWgt
             << "\n p4 (beam): " << utils::print::P4AsShortString(&info.fgP4)
             << "\n x4 (beam): " << utils::print::X4AsString(&info.fgX4)
             << "\n p4 (user): " << utils::print::P4AsShortString(&info.fgP4User)
             << "\n x4 (user): " << utils::print::X4AsString(&info.fgX4User);
        ;
#ifdef GNUMI_TEST_XY_WGT
      stream << "\n" << xypartials::GetStaticInstance();
#endif
        

  /*
  //std::cout << "GNuMIFlux::PrintCurrent ....." << std::endl;
  //LOG("Flux", pINFO) 
  LOG("Flux", pNOTICE) 
    << "Current Leaf Values: "
    << " run " << fLf_run << " evtno " << fLf_evtno << "\n"
    << " NenergyN " << fLf_NenergyN << " NWtNear " << fLf_NWtNear
    << " NenergyF " << fLf_NenergyF << " NWtFar  " << fLf_NWtFar << "\n"
    << " Norig " << fLf_Norig << " Ndecay " << fLf_Ndecay << " Ntype " << fLf_Ntype << "\n"
    << " Vxyz " << fLf_Vx << " " << fLf_Vy << " " << fLf_Vz 
    << " pdPxyz " << fLf_pdPx << " " << fLf_pdPy << " " << fLf_pdPz << "\n"
    << " pp dxdz " << fLf_ppdxdz << " dydz " << fLf_ppdydz << " pz " << fLf_pppz << "\n"
    << " pp energy " << fLf_ppenergy << " medium " << fLf_ppmedium 
    << " ptype " << fLf_ptype 
    << " ppvxyz " << fLf_ppvx << " " << fLf_ppvy << " " << fLf_ppvz << "\n"
    << " muparpxyz " << fLf_muparpx << " " << fLf_muparpy << " " << fLf_muparpz
    << " mupare " << fLf_mupare << "\n"
    << ;
  */

     return stream;
  }
}//flux
}//genie

//___________________________________________________________________________

#ifdef  GNUMI_TEST_XY_WGT

double fabserr(double a, double b) 
{ return TMath::Abs(a-b)/TMath::Max(TMath::Abs(b),1.0e-30); }

int outdiff(double a, double b, double eps, const char* label)
{
  double err = fabserr(a,b);
  if ( err > eps ) {
    std::cout << std::setw(15) << label << " err " << err 
         << " vals " << a << " " << b << std::endl;
    return 1;
  }
  return 0;
}

int gnumi2pdg(int igeant) 
{
  switch ( igeant ) {
  case 52: return -12;   // nuebar
  case 53: return  12;   // nue
  case 55: return -14;   // numubar
  case 56: return  14;   // numu
  case 58: return -16;   // nutaubar
  case 59: return  16;   // nutau
  default:
    return genie::pdg::GeantToPdg(igeant);
  }
}

void xypartials::ReadStream(std::ifstream& myfile)
{
  myfile >> parent_mass >> parentp >> parent_energy;
  myfile >> gamma >> beta_mag >> enuzr >> rad;
  myfile >> costh_pardet >> theta_pardet >> emrat >> eneu;
  myfile >> sangdet >> wgt;
  int ptype_g;
  myfile >> ptype_g;
  ptype = gnumi2pdg(ptype_g);
  if ( ptype == 13 || ptype == -13 ) {
    //std::cout << "ReadStream saw ptype " << ptype << std::endl;
    myfile >> betanu[0] >> betanu[1] >> betanu[2]
           >> p_nu[0] >> p_nu[1] >> p_nu[2];
    myfile >> partial1
           >> p_dcm_nu[0] >> p_dcm_nu[1] >> p_dcm_nu[2] >> p_dcm_nu[3];

    myfile >> muparent_px >> muparent_py >> muparent_pz;
    myfile >> gammamp >> betamp[0] >> betamp[1] >> betamp[2];
    myfile >> partial2
           >> p_pcm_mp[0] >> p_pcm_mp[1] >> p_pcm_mp[2] >> p_pcm;

    if ( p_pcm != 0.0 && p_dcm_nu[3] != 0.0 ) {
      int ntype_g;
      myfile >> costhmu >> ntype_g >> wgt_ratio;
      ntype = gnumi2pdg(ntype_g);
    }
  }
}

int xypartials::Compare(const xypartials& other) const
{
  double eps1 = 2.5e-5; // 0.003; // 6.0e-4; // 2.5e-4;
  double eps2 = 2.5e-5; // 6.0e-4; // 2.5e-4;
  double epsX = 2.5e-5; // 2.5e-4;
  int np = 0;
  np += outdiff(parent_mass  ,other.parent_mass  ,eps1,"parent_mass");
  np += outdiff(parentp      ,other.parentp      ,eps1,"parentp");
  np += outdiff(parent_energy,other.parent_energy,eps1,"parent_energy");
  np += outdiff(gamma        ,other.gamma        ,eps1,"gamma");
  np += outdiff(beta_mag     ,other.beta_mag     ,eps1,"beta_mag");
  np += outdiff(enuzr        ,other.enuzr        ,eps1,"enuzr");
  np += outdiff(rad          ,other.rad          ,eps1,"rad");
  np += outdiff(costh_pardet ,other.costh_pardet ,eps1,"costh_pardet");
  //np += outdiff(theta_pardet ,other.theta_pardet ,eps1,"theta_pardet");
  np += outdiff(emrat        ,other.emrat        ,eps1,"emrat");
  np += outdiff(eneu         ,other.eneu         ,epsX,"eneu");
  np += outdiff(sangdet      ,other.sangdet      ,eps1,"sangdet");
  np += outdiff(wgt          ,other.wgt          ,epsX,"wgt");
  if ( ptype != other.ptype ) {
    std::cout << "ptype mismatch " << ptype << " " << other.ptype << std::endl;
    np++;
  }
  if ( TMath::Abs(ptype)==13 || TMath::Abs(other.ptype)==13 ) {
    //std::cout << "== ismu " << std::endl;
    np += outdiff(betanu[0]        ,other.betanu[0]        ,eps2,"betanu[0]");
    np += outdiff(betanu[1]        ,other.betanu[1]        ,eps2,"betanu[1]");
    np += outdiff(betanu[2]        ,other.betanu[2]        ,eps2,"betanu[2]");
    np += outdiff(p_nu[0]          ,other.p_nu[0]          ,eps2,"p_nu[0]");
    np += outdiff(p_nu[1]          ,other.p_nu[1]          ,eps2,"p_nu[1]");
    np += outdiff(p_nu[2]          ,other.p_nu[2]          ,eps2,"p_nu[2]");
    np += outdiff(partial1         ,other.partial1         ,eps2,"partial1");
    np += outdiff(p_dcm_nu[0]      ,other.p_dcm_nu[0]      ,eps2,"p_dcm_nu[0]");
    np += outdiff(p_dcm_nu[1]      ,other.p_dcm_nu[1]      ,eps2,"p_dcm_nu[1]");
    np += outdiff(p_dcm_nu[2]      ,other.p_dcm_nu[2]      ,eps2,"p_dcm_nu[2]");
    np += outdiff(p_dcm_nu[3]      ,other.p_dcm_nu[3]      ,eps2,"p_dcm_nu[3]");

    np += outdiff(muparent_px      ,other.muparent_px      ,eps2,"muparent_px");
    np += outdiff(muparent_py      ,other.muparent_py      ,eps2,"muparent_py");
    np += outdiff(muparent_pz      ,other.muparent_pz      ,eps2,"muparent_pz");
    np += outdiff(gammamp          ,other.gammamp          ,eps1,"gammamp");
    np += outdiff(betamp[0]        ,other.betamp[0]        ,eps1,"betamp[0]");
    np += outdiff(betamp[1]        ,other.betamp[1]        ,eps1,"betamp[1]");
    np += outdiff(betamp[2]        ,other.betamp[2]        ,eps1,"betamp[2]");
    np += outdiff(partial2         ,other.partial2         ,eps1,"partial2");
    np += outdiff(p_pcm_mp[0]      ,other.p_pcm_mp[0]      ,eps1,"p_pcm_mp[0]");
    np += outdiff(p_pcm_mp[1]      ,other.p_pcm_mp[1]      ,eps1,"p_pcm_mp[1]");
    np += outdiff(p_pcm_mp[2]      ,other.p_pcm_mp[2]      ,eps1,"p_pcm_mp[2]");
    np += outdiff(p_pcm            ,other.p_pcm            ,eps1,"p_pcm");

    if ( ntype != other.ntype ) {
      std::cout << "ntype mismatch " << ntype << " " << other.ntype << std::endl;
      np++;
    }
    np += outdiff(costhmu          ,other.costhmu          ,eps1,"costhmu");
    np += outdiff(wgt_ratio        ,other.wgt_ratio        ,eps1,"wgt_ratio");    

  }
  return np;
}

void xypartials::Print(const Option_t* /* opt */) const
{
  std::cout << *this << std::endl;
}

namespace genie {
namespace flux {
  ostream & operator << (ostream & stream, const genie::flux::xypartials & xyp )
  {
    stream << "GNuMIFlux xypartials " << std::endl;
    stream << "  xyzdet (" << xyp.xdet << "," << xyp.ydet << "," 
           << xyp.zdet << ")" << std::endl;
    stream << "  parent: mass=" << xyp.parent_mass << " p=" << xyp.parentp 
           << " e=" << xyp.parent_energy << " gamma=" << xyp.gamma 
           << " beta_mag=" << xyp.beta_mag << std::endl;
    stream << "  enuzr=" << xyp.enuzr << " rad=" << xyp.rad 
           << " costh_pardet=" << xyp.costh_pardet << std::endl;
    stream << "  emrat=" << xyp.emrat << " sangdet=" << xyp.sangdet 
           << " wgt=" << xyp.wgt << std::endl;
    stream << "  ptype=" << xyp.ptype << " "
           << ((TMath::Abs(xyp.ptype) == 13)?"is-muon":"not-muon")
           << std::endl;

    if ( TMath::Abs(xyp.ptype)==13 ) {
      stream << "  betanu: [" << xyp.betanu[0] << "," << xyp.betanu[1] 
             << "," << xyp.betanu[2] << "]" << std::endl;
      stream << "  p_nu: [" << xyp.p_nu[0] << "," << xyp.p_nu[1] 
             << "," << xyp.p_nu[2] << "]" << std::endl;
      stream << "  partial1=" << xyp.partial1 << std::endl;
      stream << "  p_dcm_nu: [" << xyp.p_dcm_nu[0] << "," << xyp.p_dcm_nu[1] 
             << "," << xyp.p_dcm_nu[2] 
             << "," << xyp.p_dcm_nu[3] << "]" << std::endl;
      stream << "  muparent_p: [" << xyp.muparent_px << "," << xyp.muparent_py 
             << "," << xyp.muparent_pz << "]" << std::endl;
      stream << "  gammamp=" << xyp.gammamp << std::endl;
      stream << "  betamp: [" << xyp.betamp[0] << "," << xyp.betamp[1] << "," 
             << xyp.betamp[2] << "]" << std::endl;
      stream << "  partial2=" << xyp.partial2 << std::endl;
      stream << "  p_pcm_mp: [" << xyp.p_pcm_mp[0] << "," << xyp.p_pcm_mp[1] 
             << "," << xyp.p_pcm_mp[2] << "]  p_pcm=" 
             << xyp.p_pcm << std::endl;
      stream << "  ntype=" << xyp.ntype 
             << " costhmu=" << xyp.costhmu 
             << " wgt_ratio=" << xyp.wgt_ratio << std::endl;
    }
    return stream;
  }
} // flux
} // genie

xypartials& xypartials::GetStaticInstance()
{ return gpartials; }
#endif
//___________________________________________________________________________

bool GNuMIFlux::LoadConfig(string cfg)
{
  const char* altxml = gSystem->Getenv("GNUMIFLUXXML");
  if ( altxml ) {
    SetXMLFileBase(altxml);
  }
  genie::flux::GNuMIFluxXMLHelper helper(this);
  return helper.LoadConfig(cfg);
}

//___________________________________________________________________________

void GNuMIFlux::PrintConfig()
{
  
  std::ostringstream s;
  PDGCodeList::const_iterator itr = fPdgCList->begin();
  for ( ; itr != fPdgCList->end(); ++itr) s << (*itr) << " ";
  s << "[rejected: ";
  itr = fPdgCListRej->begin();
  for ( ; itr != fPdgCListRej->end(); ++itr) s << (*itr) << " ";
  s << " ] ";

  std::ostringstream fpattout;
  for (size_t i = 0; i < fNuFluxFilePatterns.size(); ++i)
    fpattout << "\n [" << std::setw(3) << i << "] " << fNuFluxFilePatterns[i];

  std::ostringstream flistout;
  std::vector<std::string> flist = GetFileList();
  for (size_t i = 0; i < flist.size(); ++i)
    flistout << "\n [" << std::setw(3) << i << "] " << flist[i];

  TLorentzVector usr0(0,0,0,0);
  TLorentzVector usr0asbeam;
  User2BeamPos(usr0,usr0asbeam);

  const int w=10, p=6;
  std::ostringstream beamrot_str, beamrotinv_str;
  beamrot_str 
    << "fBeamRot: " << std::setprecision(p) << "\n"
    << "  [ " 
    << std::setw(w) << fBeamRot.XX() << " "
    << std::setw(w) << fBeamRot.XY() << " "
    << std::setw(w) << fBeamRot.XZ() << " ]\n"
    << "  [ " 
    << std::setw(w) << fBeamRot.YX() << " "
    << std::setw(w) << fBeamRot.YY() << " "
    << std::setw(w) << fBeamRot.YZ() << " ]\n"
    << "  [ " 
    << std::setw(w) << fBeamRot.ZX() << " "
    << std::setw(w) << fBeamRot.ZY() << " "
    << std::setw(w) << fBeamRot.ZZ() << " ]";
  beamrotinv_str 
    << "fBeamRotInv: " << std::setprecision(p) << "\n"
    << "  [ " 
    << std::setw(w) << fBeamRotInv.XX() << " "
    << std::setw(w) << fBeamRotInv.XY() << " "
    << std::setw(w) << fBeamRotInv.XZ() << " ]\n"
    << "  [ " 
    << std::setw(w) << fBeamRotInv.YX() << " "
    << std::setw(w) << fBeamRotInv.YY() << " "
    << std::setw(w) << fBeamRotInv.YZ() << " ]\n"
    << "  [ " 
    << std::setw(w) << fBeamRotInv.ZX() << " "
    << std::setw(w) << fBeamRotInv.ZY() << " "
    << std::setw(w) << fBeamRotInv.ZZ() << " ]";

  LOG("Flux", pNOTICE)
    << "GNuMIFlux Config:"
    << "\n Enu_max " << fMaxEv 
    << "\n pdg-codes: " << s.str() << "\n "
    << (fG3NuMI?"g3numi":"") 
    << (fG4NuMI?"g4numi":"") 
    << (fFlugg?"flugg":"")
    << "/" << fNuFluxGen << " "
    << "(" << fNuFluxTreeName << "), " << fNEntries << " entries" 
    << " (FilePOTs " << fFilePOTs << ") "
    <<  "in " << fNFiles << " files: "
    << flistout.str()
    << "\n from file patterns:"
    << fpattout.str()
    << "\n wgt max=" << fMaxWeight << " fudge=" << fMaxWgtFudge << " using "
    << fMaxWgtEntries << " entries"
    << "\n Z0 pushback " << fZ0
    << "\n used entry " << fIEntry << " " << fIUse << "/" << fNUse
    << " times, in " << fICycle << "/" << fNCycles << " cycles"
    << "\n SumWeight " << fSumWeight << " for " << fNNeutrinos << " neutrinos"
    << "\n EffPOTsPerNu " << fEffPOTsPerNu << " AccumPOTs " << fAccumPOTs
    << "\n GenWeighted \"" << (fGenWeighted?"true":"false") << ", "
    << "\", Detector location set \"" << (fDetLocIsSet?"true":"false") << "\""
    << "\n LengthUnits " << fLengthUnits << ", scale b2u " << fLengthScaleB2U
    << ", scale u2b " << fLengthScaleU2B
    << "\n User Flux Window: "
    << "\n       " << utils::print::Vec3AsString(&(fFluxWindowPtUser[0]))
    << "\n       " << utils::print::Vec3AsString(&(fFluxWindowPtUser[1]))
    << "\n       " << utils::print::Vec3AsString(&(fFluxWindowPtUser[2]))
    << "\n Flux Window (cm, beam coord): "
    << "\n  base " << utils::print::X4AsString(&fFluxWindowBase)
    << "\n  dir1 " << utils::print::X4AsString(&fFluxWindowDir1) << " len " << fFluxWindowLen1
    << "\n  dir2 " << utils::print::X4AsString(&fFluxWindowDir2) << " len " << fFluxWindowLen2
    << "\n  normal " << utils::print::Vec3AsString(&(fWindowNormal))
    << "\n User Beam Origin: "
    << "\n  base " << utils::print::X4AsString(&fBeamZero)
    << "\n " << beamrot_str.str() << " "
    << "\n Detector Origin (beam coord): "
    << "\n  base " << utils::print::X4AsString(&usr0asbeam)
    << "\n " << beamrotinv_str.str() << " "
    << "\n UseFluxAtDetCenter " << fUseFluxAtDetCenter;

}

//___________________________________________________________________________
std::vector<std::string> GNuMIFlux::GetFileList() 
{
  std::vector<std::string> flist;
  TObjArray *fileElements=fNuFluxTree->GetListOfFiles();
  TIter next(fileElements);
  TChainElement *chEl=0;
  while (( chEl=(TChainElement*)next() )) {
    flist.push_back(chEl->GetTitle());
  }
  return flist;
}

//___________________________________________________________________________

std::vector<double> GNuMIFluxXMLHelper::GetDoubleVector(std::string str)
{
  // turn string into vector<double>
  // be liberal about separators, users might punctuate for clarity
  std::vector<std::string> strtokens = genie::utils::str::Split(str," ,;:()[]=");
  std::vector<double> vect;
  size_t ntok = strtokens.size();

  if ( fVerbose > 2 ) 
    std::cout << "GetDoubleVector \"" << str << "\"" << std::endl;

  for (size_t i=0; i < ntok; ++i) {
    std::string trimmed = utils::str::TrimSpaces(strtokens[i]);
    if ( " " == trimmed || "" == trimmed ) continue;  // skip empty strings
    double val = strtod(trimmed.c_str(), (char**)NULL);
    if ( fVerbose > 2 ) 
      std::cout << "(" << vect.size() << ") = " << val << std::endl;
    vect.push_back(val);
  }

  return vect;
}

std::vector<long int> GNuMIFluxXMLHelper::GetIntVector(std::string str)
{
  // turn string into vector<long int>
  // be liberal about separators, users might punctuate for clarity
  std::vector<std::string> strtokens = genie::utils::str::Split(str," ,;:()[]=");
  std::vector<long int> vect;
  size_t ntok = strtokens.size();

  if ( fVerbose > 2 ) 
    std::cout << "GetIntVector \"" << str << "\"" << std::endl;

  for (size_t i=0; i < ntok; ++i) {
    std::string trimmed = utils::str::TrimSpaces(strtokens[i]);
    if ( " " == trimmed || "" == trimmed ) continue;  // skip empty strings
    long int val = strtol(trimmed.c_str(),(char**)NULL,10);
    if ( fVerbose > 2 ) 
      std::cout << "(" << vect.size() << ") = " << val << std::endl;
    vect.push_back(val);
  }

  return vect;
}

bool GNuMIFluxXMLHelper::LoadConfig(string cfg)
{
  string fname = utils::xml::GetXMLFilePath(fGNuMI->GetXMLFileBase());

  bool is_accessible = ! (gSystem->AccessPathName(fname.c_str()));
  if (!is_accessible) {
    SLOG("GNuMIFlux", pERROR)
      << "The XML doc doesn't exist! (filename: " << fname << ")";
    return false;
  }

  xmlDocPtr xml_doc = xmlParseFile( fname.c_str() );
  if ( xml_doc == NULL) {
    SLOG("GNuMIFlux", pERROR)
      << "The XML doc can't be parsed! (filename: " << fname << ")";
    return false;
  }

  xmlNodePtr xml_root = xmlDocGetRootElement( xml_doc );
  if ( xml_root == NULL ) {
    SLOG("GNuMIFlux", pERROR)
      << "The XML doc is empty! (filename: " << fname << ")";
    return false;
  }

  string rootele = "gnumi_config";
  if ( xmlStrcmp(xml_root->name, (const xmlChar*)rootele.c_str() ) ) {
    SLOG("GNuMIFlux", pERROR)
      << "The XML doc has invalid root element! (filename: " << fname << ")"
      << " expected \"" << rootele << "\", saw \"" << xml_root->name << "\"";
    return false;
  }

  SLOG("GNuMIFlux", pINFO) << "Attempt to load config \"" << cfg 
                           << "\" from file: " << fname;

  bool found = this->LoadParamSet(xml_doc,cfg);

  xmlFree(xml_doc);
  return found;

}

bool GNuMIFluxXMLHelper::LoadParamSet(xmlDocPtr& xml_doc, string cfg)
{

  xmlNodePtr xml_root = xmlDocGetRootElement( xml_doc );

  // loop over all xml tree nodes that are children of the root node
  // read the entries looking for "param_set" of the right name

  // loop looking for particular config
  bool found = false;
  xmlNodePtr xml_pset = xml_root->xmlChildrenNode;
  for ( ; xml_pset != NULL ; xml_pset = xml_pset->next ) {
    if ( ! xmlStrEqual(xml_pset->name, (const xmlChar*)"param_set") ) continue;
    // every time there is a 'param_set' tag
    string param_set_name = 
      utils::str::TrimSpaces(utils::xml::GetAttribute(xml_pset,"name"));
    
    if ( param_set_name != cfg ) continue;
      
    SLOG("GNuMIFlux", pINFO) << "Found config \"" << cfg;

    this->ParseParamSet(xml_doc,xml_pset);
    found = true;

  } // loop over elements of root
  xmlFree(xml_pset);

  return found;
}

void GNuMIFluxXMLHelper::ParseParamSet(xmlDocPtr& xml_doc, xmlNodePtr& xml_pset)
{
  xmlNodePtr xml_child = xml_pset->xmlChildrenNode;
  for ( ; xml_child != NULL ; xml_child = xml_child->next ) {
    // handle basic gnumi_config/param_set
    // bad cast away const on next line, but function sig requires it
    string pname = 
      utils::xml::TrimSpaces(const_cast<xmlChar*>(xml_child->name));
    if ( pname == "text" || pname == "comment" ) continue;
    string pval  = 
      utils::xml::TrimSpaces(
              xmlNodeListGetString(xml_doc, xml_child->xmlChildrenNode, 1));

    if ( fVerbose > 1 ) 
      SLOG("GNuMIFlux", pINFO)
        << "   pname \"" << pname << "\", string value \"" << pval << "\"";

    if        ( pname == "verbose" ) {
      fVerbose = atoi(pval.c_str());

    } else if ( pname == "using_param_set" ) {
      SLOG("GNuMIFlux", pWARN) << "start using_param_set: \"" << pval << "\"";
      bool found = this->LoadParamSet(xml_doc,pval); // recurse
      if ( ! found ) {
        SLOG("GNuMIFlux", pFATAL) << "using_param_set: \"" << pval << "\" NOT FOUND";
        assert(found);
      }
      SLOG("GNuMIFlux", pWARN) << "done using_param_set: \"" << pval << "\"";
    } else if ( pname == "units" ) {
      double scale = genie::utils::units::UnitFromString(pval);
      fGNuMI->SetLengthUnits(scale);
      SLOG("GNuMIFlux", pINFO) << "set user units to \"" << pval << "\"";

    } else if ( pname == "beamdir" ) {
      ParseBeamDir(xml_doc,xml_child);
      fGNuMI->SetBeamRotation(fBeamRotXML);

    } else if ( pname == "beampos" ) {
      ParseBeamPos(pval);
      fGNuMI->SetBeamCenter(fBeamPosXML);

    } else if ( pname == "window" ) {
      ParseWindowSeries(xml_doc,xml_child);
      // RWH  !!!! MEMORY LEAK!!!!
      //std::cout << " flux window " << std::endl
      //          << " [0] " << utils::print::X4AsString(new TLorentzVector(fFluxWindowPt[0],0)) << std::endl
      //          << " [1] " << utils::print::X4AsString(new TLorentzVector(fFluxWindowPt[1],0)) << std::endl
      //          << " [2] " << utils::print::X4AsString(new TLorentzVector(fFluxWindowPt[2],0)) << std::endl;

      fGNuMI->SetFluxWindow(fFluxWindowPtXML[0],
                            fFluxWindowPtXML[1],
                            fFluxWindowPtXML[2]);

    } else if ( pname == "enumax" ) {
      ParseEnuMax(pval);

    } else if ( pname == "upstreamz" ) {
      double z0usr = -3.4e38;
      std::vector<double> v = GetDoubleVector(pval);
      if ( v.size() > 0 ) z0usr = v[0];
      fGNuMI->SetUpstreamZ(z0usr);
      SLOG("GNuMIFlux", pINFO) << "set upstreamz = " << z0usr;

    } else if ( pname == "reuse" ) {
      long int nreuse = 1;
      std::vector<long int> v = GetIntVector(pval);
      if ( v.size() > 0 ) nreuse = v[0];
      fGNuMI->SetEntryReuse(nreuse);
      SLOG("GNuMIFlux", pINFO) << "set entry reuse = " << nreuse;

    } else {
      SLOG("GNuMIFlux", pWARN)
        << "  NOT HANDLED: pname \"" << pname 
        << "\", string value \"" << pval << "\"";
      
    }

  } // loop over param_set contents
  xmlFree(xml_child);  
}

void GNuMIFluxXMLHelper::ParseBeamDir(xmlDocPtr& xml_doc, xmlNodePtr& xml_beamdir)
{
  fBeamRotXML.SetToIdentity(); // start fresh

  string dirtype = 
    utils::str::TrimSpaces(
      utils::xml::GetAttribute(xml_beamdir,"type"));

  string pval  = 
    utils::xml::TrimSpaces(
      xmlNodeListGetString(xml_doc, xml_beamdir->xmlChildrenNode, 1));

  if        ( dirtype == "series" ) {
    // series of rotations around an axis
    ParseRotSeries(xml_doc,xml_beamdir);

  } else if ( dirtype == "thetaphi3") {
    // G3 style triplet of (theta,phi) pairs
    std::vector<double> thetaphi3 = GetDoubleVector(pval);
    string units = 
      utils::str::TrimSpaces(utils::xml::GetAttribute(xml_beamdir,"units"));
    if ( thetaphi3.size() == 6 ) {
      TRotation fTempRot;
      TVector3 newX = AnglesToAxis(thetaphi3[0],thetaphi3[1],units);
      TVector3 newY = AnglesToAxis(thetaphi3[2],thetaphi3[3],units);
      TVector3 newZ = AnglesToAxis(thetaphi3[4],thetaphi3[5],units);
      fTempRot.RotateAxes(newX,newY,newZ);
      fBeamRotXML = fTempRot;  //.Inverse();
    } else {
      SLOG("GNuMIFlux", pWARN)
        << " type=\"" << dirtype << "\" within <beamdir> needs 6 values";
    }

  } else if ( dirtype == "newxyz" ) {
    // G4 style new axis values
    std::vector<double> newdir = GetDoubleVector(pval);
    if ( newdir.size() == 9 ) {
      TRotation fTempRot;
      TVector3 newX = TVector3(newdir[0],newdir[1],newdir[2]).Unit();
      TVector3 newY = TVector3(newdir[3],newdir[4],newdir[5]).Unit();
      TVector3 newZ = TVector3(newdir[6],newdir[7],newdir[8]).Unit();
      fTempRot.RotateAxes(newX,newY,newZ);
      fBeamRotXML = fTempRot.Inverse(); // weirdly necessary: frame vs. obj rot
    } else {
      SLOG("GNuMIFlux", pWARN)
        << " type=\"" << dirtype << "\" within <beamdir> needs 9 values";
    }

  } else {
    // yet something else ... what? 3 choices weren't sufficient?
    SLOG("GNuMIFlux", pWARN)
      << " UNHANDLED type=\"" << dirtype << "\" within <beamdir>";
  }

  if ( fVerbose > 1 ) {
    int w=10, p=6;
    std::cout << " fBeamRotXML: " << std::setprecision(p) << std::endl;
    std::cout << " [ " 
              << std::setw(w) << fBeamRotXML.XX() << " "
              << std::setw(w) << fBeamRotXML.XY() << " "
              << std::setw(w) << fBeamRotXML.XZ() << endl
              << "   " 
              << std::setw(w) << fBeamRotXML.YX() << " "
              << std::setw(w) << fBeamRotXML.YY() << " "
              << std::setw(w) << fBeamRotXML.YZ() << endl
              << "   " 
              << std::setw(w) << fBeamRotXML.ZX() << " "
              << std::setw(w) << fBeamRotXML.ZY() << " "
              << std::setw(w) << fBeamRotXML.ZZ() << " ] " << std::endl;
    std::cout << std::endl;
  }

}

void GNuMIFluxXMLHelper::ParseBeamPos(std::string str)
{
  std::vector<double> xyz = GetDoubleVector(str);
  if ( xyz.size() == 3 ) {
    fBeamPosXML = TVector3(xyz[0],xyz[1],xyz[2]);
  } else if ( xyz.size() == 6 ) {
    // should check for '=' between triplets but we won't be so pedantic
    // ( userx, usery, userz ) = ( beamx, beamy, beamz )
    TVector3 userpos(xyz[0],xyz[1],xyz[2]);
    TVector3 beampos(xyz[3],xyz[4],xyz[5]);
    fBeamPosXML = userpos - fBeamRotXML*beampos;
  } else {
    SLOG("GNuMIFlux", pWARN)
      << "Unable to parse " << xyz.size() << " values in <beampos>";
    return;
   }
  if ( fVerbose > 1 ) {
    int w=16, p=10;
    std::cout << " fBeamPosXML: [ " << std::setprecision(p) 
              << std::setw(w) << fBeamPosXML.X() << " , "
              << std::setw(w) << fBeamPosXML.Y() << " , "
              << std::setw(w) << fBeamPosXML.Z() << " ] "
              << std::endl;
  }
}

void GNuMIFluxXMLHelper::ParseEnuMax(std::string str)
{
  std::vector<double> v = GetDoubleVector(str);
  size_t n = v.size();
  if ( n > 0 ) {
    fGNuMI->SetMaxEnergy(v[0]);
    if ( fVerbose > 1 ) 
      std::cout << "ParseEnuMax SetMaxEnergy(" << v[0] << ") " << std::endl;
  }
  if ( n > 1 ) {
    fGNuMI->SetMaxEFudge(v[1]);
    if ( fVerbose > 1 ) 
      std::cout << "ParseEnuMax SetMaxEFudge(" << v[1] << ")" << std::endl;
  }
  if ( n > 2 ) {
    if ( n == 3 ) {
      fGNuMI->SetMaxWgtScan(v[2]);
      if ( fVerbose > 1 ) 
        std::cout << "ParseEnuMax SetMaxWgtScan(" << v[2] << ")" << std::endl;
    } else {
      long int nentries = (long int)v[3];
      fGNuMI->SetMaxWgtScan(v[2],nentries);
      if ( fVerbose > 1 ) 
        std::cout << "ParseEnuMax SetMaxWgtScan(" << v[2] << "," << nentries << ")" << std::endl;
    }
  }
}

void GNuMIFluxXMLHelper::ParseRotSeries(xmlDocPtr& xml_doc, xmlNodePtr& xml_pset)
{
  TRotation fTempRot; // reset matrix

  xmlNodePtr xml_child = xml_pset->xmlChildrenNode;
  for ( ; xml_child != NULL ; xml_child = xml_child->next ) {
    // in a <beamdir> of type "series"
    // should be a sequence of <rotation> entries
    string name = 
      utils::xml::TrimSpaces(const_cast<xmlChar*>(xml_child->name));
    if ( name == "text" || name == "comment" ) continue;

    if ( name == "rotation" ) {
      string val = utils::xml::TrimSpaces(
          xmlNodeListGetString(xml_doc, xml_child->xmlChildrenNode, 1));
      string axis = 
        utils::str::TrimSpaces(utils::xml::GetAttribute(xml_child,"axis"));

      string units = 
        utils::str::TrimSpaces(utils::xml::GetAttribute(xml_child,"units"));

      double rot = atof(val.c_str());
      // assume radians unless given a hint that it's degrees
      if ( 'd' == units[0] || 'D' == units[0] ) rot *= TMath::DegToRad();

      if ( fVerbose > 0 )
        SLOG("GNuMIFlux", pINFO)
          << " rotate " << rot << " radians around " << axis << " axis";

      if      ( axis[0] == 'x' || axis[0] == 'X' ) fTempRot.RotateX(rot);
      else if ( axis[0] == 'y' || axis[0] == 'Y' ) fTempRot.RotateY(rot);
      else if ( axis[0] == 'z' || axis[0] == 'Z' ) fTempRot.RotateZ(rot);
      else {
        SLOG("GNuMIFlux", pINFO)
          << " no " << axis << " to rotate around";
      }

    } else {
      SLOG("GNuMIFlux", pWARN)
        << " found <" << name << "> within <beamdir type=\"series\">";
    }
  }
  // TRotation rotates objects not frames, so we want the inverse
  fBeamRotXML = fTempRot.Inverse();
  xmlFree(xml_child);  
}

void GNuMIFluxXMLHelper::ParseWindowSeries(xmlDocPtr& xml_doc, xmlNodePtr& xml_pset)
{
  int ientry = -1;

  xmlNodePtr xml_child = xml_pset->xmlChildrenNode;
  for ( ; xml_child != NULL ; xml_child = xml_child->next ) {
    // in a <windowr> element
    // should be a sequence of <point> entries
    string name = 
      utils::xml::TrimSpaces(const_cast<xmlChar*>(xml_child->name));
    if ( name == "text" || name == "comment" ) continue;

    if ( name == "point" ) {
      string val  = 
        utils::xml::TrimSpaces(
          xmlNodeListGetString(xml_doc, xml_child->xmlChildrenNode, 1));
      string coord = 
        utils::str::TrimSpaces(utils::xml::GetAttribute(xml_child,"coord"));

      std::vector<double> xyz = GetDoubleVector(val);
      if ( xyz.size() != 3 || coord != "det" ) {
        SLOG("GNuMIFlux", pWARN)
          << "parsing <window> found <point> but size=" << xyz.size()
          << " (expect 3) and coord=\"" << coord << "\" (expect \"det\")"
          << " IGNORE problem";
      }
      ++ientry;
      if ( ientry < 3 && ientry >= 0 ) {
        TVector3 pt(xyz[0],xyz[1],xyz[2]);
        if ( fVerbose > 0 ) {
          int w=16, p=10;
          std::cout << " point[" << ientry <<"] = [ " << std::setprecision(p) 
                    << std::setw(w) << pt.X() << " , "
                    << std::setw(w) << pt.Y() << " , "
                    << std::setw(w) << pt.Z() << " ] "
                    << std::endl;
        }
        fFluxWindowPtXML[ientry] = pt;  // save the point
      } else {
        SLOG("GNuMIFlux", pWARN)
          << " <window><point> ientry " << ientry << " out of range (0-2)";
      }

    } else {
      SLOG("GNuMIFlux", pWARN)
        << " found <" << name << "> within <window>";
    }
  }
  xmlFree(xml_child);  
}

TVector3 GNuMIFluxXMLHelper::AnglesToAxis(double theta, double phi, std::string units)
{
  double xyz[3];
  // assume radians unless given a hint that it's degrees
  double scale = ('d'==units[0]||'D'==units[0]) ? TMath::DegToRad() : 1.0 ;

  xyz[0] = TMath::Cos(scale*phi)*TMath::Sin(scale*theta);
  xyz[1] = TMath::Sin(scale*phi)*TMath::Sin(scale*theta);
  xyz[2] = TMath::Cos(scale*theta);
  // condition vector to eliminate most floating point errors
  for (int i=0; i<3; ++i) {
    const double eps = 1.0e-15;
    if (TMath::Abs(xyz[i])   < eps ) xyz[i] =  0;
    if (TMath::Abs(xyz[i]-1) < eps ) xyz[i] =  1;
    if (TMath::Abs(xyz[i]+1) < eps ) xyz[i] = -1;
  }
  return TVector3(xyz[0],xyz[1],xyz[2]);                    
}

TVector3 GNuMIFluxXMLHelper::ParseTV3(const string& str)
{
  std::vector<double> xyz = GetDoubleVector(str);
  if ( xyz.size() != 3 ) {
    return TVector3();
    SLOG("GNuMIFlux", pWARN)
      << " ParseTV3 \"" << str << "\" had " << xyz.size() << " elements ";
  }
  return TVector3(xyz[0],xyz[1],xyz[2]);

}
//___________________________________________________________________________

#ifndef SKIP_MINERVA_MODS
//=========================================
// The following was inserted by MINERvA
//=========================================
int GNuMIFluxPassThroughInfo::getVolID(TString sval){
  int ival=0;
  if(sval=="AddedLV")ival=1;
  else if(sval=="AlTube1LV")ival=2;
  else if(sval=="AlTube2LV")ival=3;
  else if(sval=="Al_BLK1")ival=4;
  else if(sval=="Al_BLK2")ival=5;
  else if(sval=="Al_BLK3")ival=6;
  else if(sval=="Al_BLK4")ival=7;
  else if(sval=="Al_BLK5")ival=8;
  else if(sval=="Al_BLK6")ival=9;
  else if(sval=="Al_BLK7")ival=10;
  else if(sval=="Al_BLK8")ival=11;
  else if(sval=="AlholeL")ival=12;
  else if(sval=="AlholeR")ival=13;
  else if(sval=="BEndLV")ival=14;
  else if(sval=="BFrontLV")ival=15;
  else if(sval=="BeDWLV")ival=16;
  else if(sval=="BeUp1LV")ival=17;
  else if(sval=="BeUp2LV")ival=18;
  else if(sval=="BeUp3LV")ival=19;
  else if(sval=="BodyLV")ival=20;
  else if(sval=="BudalMonitor")ival=21;
  else if(sval=="CLid1LV")ival=22;
  else if(sval=="CLid2LV")ival=23;
  else if(sval=="CShld_BLK11")ival=24;
  else if(sval=="CShld_BLK12")ival=25;
  else if(sval=="CShld_BLK2")ival=26;
  else if(sval=="CShld_BLK3")ival=27;
  else if(sval=="CShld_BLK4")ival=28;
  else if(sval=="CShld_BLK7")ival=29;
  else if(sval=="CShld_BLK8")ival=30;
  else if(sval=="CShld_stl,BLK")ival=31;
  else if(sval=="CerTubeLV")ival=32;
  else if(sval=="CeramicRod")ival=33;
  else if(sval=="ConcShield")ival=34;
  else if(sval=="Concrete Chase Section")ival=35;
  else if(sval=="Conn1LV")ival=36;
  else if(sval=="Conn2LV")ival=37;
  else if(sval=="Conn3LV")ival=38;
  else if(sval=="DNWN")ival=39;
  else if(sval=="DPIP")ival=40;
  else if(sval=="DVOL")ival=41;
  else if(sval=="DuratekBlock")ival=42;
  else if(sval=="DuratekBlockCovering")ival=43;
  else if(sval=="HadCell")ival=44;
  else if(sval=="HadronAbsorber")ival=45;
  else if(sval=="MuMonAlcvFill_0")ival=46;
  else if(sval=="MuMonAlcv_0")ival=47;
  else if(sval=="MuMonAlcv_1")ival=48;
  else if(sval=="MuMonAlcv_2")ival=49;
  else if(sval=="MuMon_0")ival=50;
  else if(sval=="MuMon_1")ival=51;
  else if(sval=="MuMon_2")ival=52;
  else if(sval=="PHorn1CPB1slv")ival=53;
  else if(sval=="PHorn1CPB2slv")ival=54;
  else if(sval=="PHorn1F")ival=55;
  else if(sval=="PHorn1Front")ival=56;
  else if(sval=="PHorn1IC")ival=57;
  else if(sval=="PHorn1InsRingslv")ival=58;
  else if(sval=="PHorn1OC")ival=59;
  else if(sval=="PHorn2CPB1slv")ival=60;
  else if(sval=="PHorn2CPB2slv")ival=61;
  else if(sval=="PHorn2F")ival=62;
  else if(sval=="PHorn2Front")ival=63;
  else if(sval=="PHorn2IC")ival=64;
  else if(sval=="PHorn2InsRingslv")ival=65;
  else if(sval=="PHorn2OC")ival=66;
  else if(sval=="PVHadMon")ival=67;
  else if(sval=="Pipe1")ival=68;
  else if(sval=="Pipe1_water")ival=69;
  else if(sval=="Pipe2")ival=70;
  else if(sval=="Pipe2_water")ival=71;
  else if(sval=="Pipe3")ival=72;
  else if(sval=="Pipe3_water")ival=73;
  else if(sval=="Pipe4")ival=74;
  else if(sval=="Pipe4_water")ival=75;
  else if(sval=="Pipe5")ival=76;
  else if(sval=="Pipe5_water")ival=77;
  else if(sval=="Pipe6")ival=78;
  else if(sval=="Pipe6_water")ival=79;
  else if(sval=="Pipe7")ival=80;
  else if(sval=="Pipe8")ival=81;
  else if(sval=="Pipe8_water")ival=82;
  else if(sval=="Pipe9")ival=83;
  else if(sval=="PipeAdapter1")ival=84;
  else if(sval=="PipeAdapter1_water")ival=85;
  else if(sval=="PipeAdapter2")ival=86;
  else if(sval=="PipeAdapter2_water")ival=87;
  else if(sval=="PipeBellowB")ival=88;
  else if(sval=="PipeBellowB_water")ival=89;
  else if(sval=="PipeBellowT")ival=90;
  else if(sval=="PipeBellowT_water")ival=91;
  else if(sval=="PipeC1")ival=92;
  else if(sval=="PipeC1_water")ival=93;
  else if(sval=="PipeC2")ival=94;
  else if(sval=="PipeC2_water")ival=95;
  else if(sval=="ROCK")ival=96;
  else if(sval=="Ring1LV")ival=97;
  else if(sval=="Ring2LV")ival=98;
  else if(sval=="Ring3LV")ival=99;
  else if(sval=="Ring4LV")ival=100;
  else if(sval=="Ring5LV")ival=101;
  else if(sval=="SC01")ival=102;
  else if(sval=="SpiderSupport")ival=103;
  else if(sval=="Stl_BLK1")ival=104;
  else if(sval=="Stl_BLK10")ival=105;
  else if(sval=="Stl_BLK2")ival=106;
  else if(sval=="Stl_BLK3")ival=107;
  else if(sval=="Stl_BLK4")ival=108;
  else if(sval=="Stl_BLK5")ival=109;
  else if(sval=="Stl_BLK6")ival=110;
  else if(sval=="Stl_BLK7")ival=111;
  else if(sval=="Stl_BLK8")ival=112;
  else if(sval=="Stl_BLK9")ival=113;
  else if(sval=="Stlhole")ival=114;
  else if(sval=="TGAR")ival=115;
  else if(sval=="TGT1")ival=116;
  else if(sval=="TGTExitCyl2LV")ival=117;
  else if(sval=="TUNE")ival=118;
  else if(sval=="Tube1aLV")ival=119;
  else if(sval=="Tube1bLV")ival=120;
  else if(sval=="UpWn1")ival=121;
  else if(sval=="UpWn2")ival=122;
  else if(sval=="UpWnAl1SLV")ival=123;
  else if(sval=="UpWnAl2SLV")ival=124;
  else if(sval=="UpWnAl3SLV")ival=125;
  else if(sval=="UpWnFe1SLV")ival=126;
  else if(sval=="UpWnFe2SLV")ival=127;
  else if(sval=="UpWnPolyCone")ival=128;
  else if(sval=="blu_BLK25")ival=129;
  else if(sval=="blu_BLK26")ival=130;
  else if(sval=="blu_BLK27")ival=131;
  else if(sval=="blu_BLK28")ival=132;
  else if(sval=="blu_BLK29")ival=133;
  else if(sval=="blu_BLK32")ival=134;
  else if(sval=="blu_BLK37")ival=135;
  else if(sval=="blu_BLK38")ival=136;
  else if(sval=="blu_BLK39")ival=137;
  else if(sval=="blu_BLK40")ival=138;
  else if(sval=="blu_BLK45")ival=139;
  else if(sval=="blu_BLK46")ival=140;
  else if(sval=="blu_BLK47")ival=141;
  else if(sval=="blu_BLK48")ival=142;
  else if(sval=="blu_BLK49")ival=143;
  else if(sval=="blu_BLK50")ival=144;
  else if(sval=="blu_BLK51")ival=145;
  else if(sval=="blu_BLK53")ival=146;
  else if(sval=="blu_BLK55")ival=147;
  else if(sval=="blu_BLK57")ival=148;
  else if(sval=="blu_BLK59")ival=149;
  else if(sval=="blu_BLK61")ival=150;
  else if(sval=="blu_BLK63")ival=151;
  else if(sval=="blu_BLK64")ival=152;
  else if(sval=="blu_BLK65")ival=153;
  else if(sval=="blu_BLK66")ival=154;
  else if(sval=="blu_BLK67")ival=155;
  else if(sval=="blu_BLK68")ival=156;
  else if(sval=="blu_BLK69")ival=157;
  else if(sval=="blu_BLK70")ival=158;
  else if(sval=="blu_BLK72")ival=159;
  else if(sval=="blu_BLK73")ival=160;
  else if(sval=="blu_BLK75")ival=161;
  else if(sval=="blu_BLK77")ival=162;
  else if(sval=="blu_BLK78")ival=163;
  else if(sval=="blu_BLK79")ival=164;
  else if(sval=="blu_BLK81")ival=165;
  else if(sval=="conc_BLK")ival=166;
  else if(sval=="pvBaffleMother")ival=167;
  else if(sval=="pvDPInnerTrackerTube")ival=168;
  else if(sval=="pvMHorn1Mother")ival=169;
  else if(sval=="pvMHorn2Mother")ival=170;
  else if(sval=="pvTargetMother")ival=171;
  else if(sval=="stl_slab1")ival=172;
  else if(sval=="stl_slab4")ival=173;
  else if(sval=="stl_slab5")ival=174;
  else if(sval=="stl_slabL")ival=175;
  else if(sval=="stl_slabR")ival=176;
  return ival;
}

int GNuMIFluxPassThroughInfo::getProcessID(TString sval){
  int ival=0;
  if(sval=="AntiLambdaInelastic")ival=1;
  else if(sval=="AntiNeutronInelastic")ival=2;
  else if(sval=="AntiOmegaMinusInelastic")ival=3;
  else if(sval=="AntiProtonInelastic")ival=4;
  else if(sval=="AntiSigmaMinusInelastic")ival=5;
  else if(sval=="AntiSigmaPlusInelastic")ival=6;
  else if(sval=="AntiXiMinusInelastic")ival=7;
  else if(sval=="AntiXiZeroInelastic")ival=8;
  else if(sval=="Decay")ival=9;
  else if(sval=="KaonMinusInelastic")ival=10;
  else if(sval=="KaonPlusInelastic")ival=11;
  else if(sval=="KaonZeroLInelastic")ival=12;
  else if(sval=="KaonZeroSInelastic")ival=13;
  else if(sval=="LambdaInelastic")ival=14;
  else if(sval=="NeutronInelastic")ival=15;
  else if(sval=="OmegaMinusInelastic")ival=16;
  else if(sval=="PionMinusInelastic")ival=17;
  else if(sval=="PionPlusInelastic")ival=18;
  else if(sval=="Primary")ival=19;
  else if(sval=="ProtonInelastic")ival=20;
  else if(sval=="SigmaMinusInelastic")ival=21;
  else if(sval=="SigmaPlusInelastic")ival=22;
  else if(sval=="XiMinusInelastic")ival=23;
  else if(sval=="XiZeroInelastic")ival=24;
  else if(sval=="hElastic")ival=25;
  return ival;
}
//END of minerva additions
#endif

