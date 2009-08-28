//____________________________________________________________________________
/*
 Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Robert Hatcher <rhatcher@fnal.gov>
         Fermi National Accelerator Laboratory

         Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory

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
*/
//____________________________________________________________________________

#include <cstdlib>
#include <fstream>
#include <vector>
#include <sstream>
#include <cassert>

#include "libxml/xmlmemory.h"
#include "libxml/parser.h"

#include "Utils/XmlParserUtils.h"
#include "Utils/StringUtils.h"

#include <TFile.h>
#include <TChain.h>
#include <TSystem.h>
#include <TStopwatch.h>

#include "Conventions/Units.h"
#include "Conventions/GBuild.h"

#include "FluxDrivers/GNuMIFlux.h"
#include "FluxDrivers/GNuMINtuple/g3numi.h"
#include "FluxDrivers/GNuMINtuple/g3numi.C"
#include "FluxDrivers/GNuMINtuple/g4numi.h"
#include "FluxDrivers/GNuMINtuple/g4numi.C"

#include "Messenger/Messenger.h"
#include "Numerical/RandomGen.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGCodeList.h"
#include "Utils/MathUtils.h"
#include "Utils/PrintUtils.h"
#include "Utils/UnitUtils.h"

using std::endl;

#include <vector>
#include <algorithm>
#include <iomanip>
#include "TRegexp.h"
#include "TString.h"

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
      TVector3    fBeamPos;
      TRotation   fBeamRot;
      TVector3    fFluxWindowPt[3];
    };
  }
}

ClassImp(GNuMIFluxPassThroughInfo)

// some nominal positions used in the original g3 ntuple
const TLorentzVector kPosCenterNearBeam(0.,0.,  1039.35,0.);
const TLorentzVector kPosCenterFarBeam (0.,0.,735340.00,0.);

//____________________________________________________________________________
GNuMIFlux::GNuMIFlux()
{
  this->Initialize();
}
//___________________________________________________________________________
GNuMIFlux::~GNuMIFlux()
{
  this->CleanUp();
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
       LOG("Flux", pERROR)
         << "** Fractional weight = " << f 
         << " > 1 !! Bump fMaxWeight estimate."
         << PassThroughInfo();
       fMaxWeight = this->Weight() * fMaxWgtFudge; // bump the weight
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
  if ( ! fG3NuMI && ! fG4NuMI ) {
     LOG("Flux", pERROR)
          << "The flux driver has not been properly configured";
     return false;	
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
      if (fICycle < fNCycles || fNCycles == 0 ) {
        fICycle++;
        fIEntry=0;
      } else {
        LOG("Flux", pWARN)
          << "No more entries in input flux neutrino ntuple, cycle "
          << fICycle << " of " << fNCycles;
        fEnd = true;
        //assert(0);
        return false;	
      }
    }
    
    if ( fG3NuMI ) { 
      fG3NuMI->GetEntry(fIEntry); 
      fCurEntry->MakeCopy(fG3NuMI); 
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
      LOG("Flux",pDEBUG) 
        << "got " << fNNeutrinos << " new fIEntry " << fIEntry 
        << " evtno " << fCurEntry->evtno;
#endif
    } else { 
      fG4NuMI->GetEntry(fIEntry); 
      fCurEntry->MakeCopy(fG4NuMI); 
    }

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

  if ( ! fPdgCList->ExistsInPDGCodeList(fCurEntry->fgPdgC) ) {
     LOG("Flux", pWARN)
          << "Unknown decay mode or decay mode producing an undeclared"
          << " neutrino species: "
          << fCurEntry->ntype 
          << " (pcodes=" << fCurEntry->pcodes << ")"
          << "\nDeclared list of neutrino species: " << *fPdgCList;
     //exit(123);
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
  fWeight = fCurEntry->nimpwt * fCurEntry->fgXYWgt;  // full weight

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
  TVector3 dirNu = fCurEntry->fgX4.Vect() - fgX4dkvtx.Vect();
  double dirnorm = 1.0 / dirNu.Mag();
  fCurEntry->fgP4.SetPxPyPzE( Ev*dirnorm*dirNu.X(), 
                              Ev*dirnorm*dirNu.Y(),
                              Ev*dirnorm*dirNu.Z(), Ev);

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
    << "\n p4: " << utils::print::P4AsShortString(&fCurEntry->fgP4)
    << "\n x4: " << utils::print::X4AsString(&fCurEntry->fgX4);
#endif
  if ( Ev > fMaxEv ) {
    LOG("Flux", pFATAL)
      << "Generated neutrino had E_nu = " << Ev << " > " << fMaxEv 
      << " maximum ";
    assert(0);
  }


  // update the # POTs, sum of weights & number of neutrinos 
  fAccumPOTs += fEffPOTsPerNu / fMaxWeight;
  fSumWeight += this->Weight();
  fNNeutrinos++;

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
  fEffPOTsPerNu = area_ratio * ( fFilePOTs / fNEntries );
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
void GNuMIFlux::LoadBeamSimData(string filename, string det_loc)
{
// Loads in a gnumi beam simulation root file (converted from hbook format)
// into the GNuMIFlux driver.

  bool found_cfg = this->LoadConfig(det_loc);
  if ( ! found_cfg ) {
    LOG("Flux", pFATAL) 
      << "LoadBeamSimData could not find XML config \"" << det_loc << "\"\n";
    exit(1);
  }

  fNuFluxFilePattern = filename;
  LOG("Flux", pNOTICE)
        << "Loading gnumi flux tree from ROOT file(s): " << filename;

  // !WILDCARD only works for file name ... NOT directory
  string dirname = gSystem->UnixPathName(gSystem->WorkingDirectory());
  size_t slashpos = filename.find_last_of("/");
  size_t fbegin;
  if ( slashpos != std::string::npos ) {
    dirname = filename.substr(0,slashpos);
    LOG("Flux", pINFO) << "Look for flux using directory " << dirname;
    fbegin = slashpos + 1;
  } else { fbegin = 0; }

  void* dirp = gSystem->OpenDirectory(gSystem->ExpandPathName(dirname.c_str()));
    // create a (sortable) vector of file names
  std::vector<std::string> fnames;
  if ( dirp ) {
    std::string basename = 
      filename.substr(fbegin,filename.size()-fbegin);
    TRegexp re(basename.c_str(),kTRUE);
    const char* onefile;
    while ( ( onefile = gSystem->GetDirEntry(dirp) ) ) {
      TString afile = onefile;
      if ( afile=="." || afile==".." ) continue;
      if ( basename!=afile && afile.Index(re) == kNPOS ) continue;
      std::string fullname = dirname + "/" + afile.Data();
      fnames.push_back(fullname);
    }
    gSystem->FreeDirectory(dirp);
    // sort the list
    std::sort(fnames.begin(),fnames.end());
    for (unsigned int indx=0; indx < fnames.size(); ++indx) {
      //std::cout << "  [" << std::setw(3) << indx << "]  \"" 
      //          << fnames[indx] << "\"" << std::endl;
      bool isok = ! (gSystem->AccessPathName(fnames[indx].c_str()));
      if ( isok ) {
        // open the file to see what it contains
        TFile tf(fnames[indx].c_str());
        const std::string tnames[] = { "h10", "nudata" };
        for (int j = 0; j < 2 ; ++j ) { 
          TTree* atree = (TTree*)tf.Get(tnames[j].c_str());
          if ( atree ) {
            // create chain if none exists
            if ( ! fNuFluxTree ) {
              this->SetTreeName(tnames[j]);
              fNuFluxTree = new TChain(fNuFluxTreeName.c_str());
              // here we should scan for estimated POTs/file
              // also maybe the check on max weight
            }
            // sanity check for mixing g3/g4 files
            if ( fNuFluxTreeName !=  tnames[j] ) {
              LOG("Flux", pFATAL) 
                << "Inconsistent flux file types\n"
                << "The input gnumi flux file \"" << fnames[indx] 
                << "\"\ncontains a '" << tnames[j] << "' g" << int(j+3)
                << "numi ntuple " 
                // 0->1  1->0 with 1-j,  0->4 1->3 with 4-j 
                << "but a '" << tnames[1-j] << "' g" << int(4-j)
                << "numi ntuple has alread been seen in the chain";
              exit(1);
            } // sanity mix/match g3/g4
            // add the file to the chain
            this->AddFile(atree,fnames[indx]);
          } // found a tree
        } // loop over either g3 or g4
        tf.Close();
      }
    } // loop over sorted file names
  } // legal directory
  if ( fNuFluxTreeName == "" ) {
    LOG("Flux", pFATAL)
     << "The input gnumi flux file doesn't exist! Initialization failed!";
    exit(1);
  }
  if ( fNuFluxTreeName == "h10"    ) fG3NuMI = new g3numi(fNuFluxTree);
  if ( fNuFluxTreeName == "nudata" ) fG4NuMI = new g4numi(fNuFluxTree);

#ifdef OLD_STUFF
  bool is_accessible = ! (gSystem->AccessPathName( filename.c_str() ));
  if (!is_accessible) {
    LOG("Flux", pFATAL)
     << "The input gnumi flux file doesn't exist! Initialization failed!";
    exit(1);
  }

  fNuFluxFile = new TFile(filename.c_str(), "read");
  if (fNuFluxFile) {
      LOG("Flux", pINFO) << "Getting flux tree: " << fNuFluxTreeName;
      fNuFluxTree = (TTree*) fNuFluxFile->Get(fNuFluxTreeName.c_str());
      if (!fNuFluxTree) {
          LOG("Flux", pERROR)
             << "** Couldn't get flux tree: " << fNuFluxTreeName;
          return;
      }
  } else {
      LOG("Flux", pERROR) << "** Couldn't open: " << filename;
      return;
  }
#endif

  // this will open all files and read header!!
  fNEntries = fNuFluxTree->GetEntries();

  LOG("Flux", pNOTICE)
      << "Loaded flux tree contains " <<  fNEntries << " entries";

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
    if ( ! fG3NuMI ) strindx += 2;
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
                       << wgtgenmx << ", energy = " << enumx;

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
void GNuMIFlux::SetFluxParticles(const PDGCodeList & particles)
{
  if (!fPdgCList) {
     fPdgCList = new PDGCodeList;
  }
  fPdgCList->Copy(particles);

  LOG("Flux", pINFO)
    << "Declared list of neutrino species: " << *fPdgCList;
}
//___________________________________________________________________________
void GNuMIFlux::SetMaxEnergy(double Ev)
{
  fMaxEv = TMath::Max(0.,Ev);

  LOG("Flux", pINFO)
    << "Declared maximum flux neutrino energy: " << fMaxEv;
}
//___________________________________________________________________________
void GNuMIFlux::SetUpstreamZ(double z0)
{
// The flux neutrino position (x,y) is given on the user specified flux window.
// This method sets the preferred user coord starting z position upstream of
// detector face. Each flux neutrino will be backtracked from the initial
// flux window to the input z0.  If the value is unreasonable (> 10^30) 
// then the ray is left on the flux window.

  fZ0 = z0;
}
//___________________________________________________________________________
void GNuMIFlux::SetNumOfCycles(long int ncycle)
{
// The flux ntuples can be recycled for a number of times to boost generated
// event statistics without requiring enormous beam simulation statistics.
// That option determines how many times the driver is going to cycle through
// the input flux ntuple.
// With ncycle=0 the flux ntuple will be recycled an infinite amount of times so
// that the event generation loop can exit only on a POT or event num check.

  fNCycles = TMath::Max(0L, ncycle);
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
void GNuMIFlux::Initialize(void)
{
  LOG("Flux", pNOTICE) << "Initializing GNuMIFlux driver";

  fMaxEv           =  0;
  fEnd             =  false;
  fPdgCList        = new PDGCodeList;
  fCurEntry    = new GNuMIFluxPassThroughInfo;

  fNuFluxTree      =  0;
  fG3NuMI          =  0;
  fG4NuMI          =  0;
  fNuFluxTreeName  = "";
  fNFiles          =  0;

  fNEntries        =  0;
  fIEntry          = -1;
  fNCycles         =  0;
  fICycle          =  0;
  fNUse            =  1;
  fIUse            =  999999;

  fNuTot           = 0;
  fFilePOTs        = 0;

  fMaxWeight       = -1;
  fMaxWgtFudge     =  1.05;
  fMaxWgtEntries   = 2500000;
  fMaxEFudge       =  0;

  fZ0              =  -3.4e38;
  fSumWeight       =  0;
  fNNeutrinos      =  0;
  fEffPOTsPerNu    =  0;
  fAccumPOTs       =  0;

  fGenWeighted     = false;
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

  if (fPdgCList) delete fPdgCList;
  if (fCurEntry) delete fCurEntry;

  if ( fG3NuMI ) delete fG3NuMI;
  if ( fG4NuMI ) delete fG4NuMI;

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
  Int_t evtno = 0;
  TBranch* br_evtno = 0;
  thetree->SetBranchAddress("evtno",&evtno, &br_evtno);
  thetree->GetEntry(0);
  Int_t evt_1 = evtno;
  Int_t est_1 = (TMath::FloorNint(evt_1/100.))*100 + 1;
  thetree->GetEntry(nentries-1);
  Int_t evt_N = evtno;
  Int_t est_N = (TMath::FloorNint((evt_N-1)/100.)+1)*100;
  ULong64_t npots = est_N - est_1 + 1;

  LOG("Flux",pNOTICE) //INFO)
    << fNuFluxTreeName << "->AddFile() of " << nentries << " entries ["
    << evt_1 << ":" << evt_N << "(" <<  est_1 << ":" << est_N << ")=" 
    << npots <<" POTs] in file: " << fname;
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
    std::cout << "ptype " << this->ptype << " m " << parent_mass 
              << " p " << parentp << " e " << parent_energy << " gamma " << gamma
              << " beta " << beta_mag << std::endl;

    std::cout << " enuzr " << enuzr << " rad " << rad << " costh " << costh_pardet
              << " theta " << theta_pardet << " emrat " << emrat 
              << " enu " << enu << std::endl;
  }

#ifdef  GNUMI_TEST_XY_WGT
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
  double sangdet = ( kRDET*kRDET / ( (zpos-this->vz)*(zpos-this->vz)))/4.0;

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
             << "\n norg " << info.norig << " ndecay " << info.ndecay
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
      stream << "\nCurrent: pdg " << info.fgPdgC 
             << " xywgt " << info.fgXYWgt
             << "\n p4 (beam): " << utils::print::P4AsShortString(&info.fgP4)
             << "\n x4 (beam): " << utils::print::X4AsString(&info.fgX4)
             << "\n p4 (user): " << utils::print::P4AsShortString(&info.fgP4User)
             << "\n x4 (user): " << utils::print::X4AsString(&info.fgX4User);
        ;
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

void xypartials::ReadStream(ifstream& myfile)
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
  std::cout << "GNuMIFlux xypartials " << std::endl;
  std::cout << "  parent: mass=" << parent_mass << " p=" << parentp 
            << " e=" << parent_energy << " gamma=" << gamma << " beta_mag=" << beta_mag << std::endl;
  std::cout << "  enuzr=" << enuzr << " rad=" << rad 
            << " costh_pardet=" << costh_pardet << std::endl;
  std::cout << "  emrat=" << emrat << " sangdet=" << sangdet << " wgt=" << wgt << std::endl;
  std::cout << "  ptype=" << ptype
            << ((TMath::Abs(ptype) == 13)?"is-muon":"not-muon")
            << std::endl;

  if ( TMath::Abs(ptype)==13 ) {
    std::cout << "  betanu: [" << betanu[0] << "," << betanu[1] << "," << betanu[2] << "]" << std::endl;
    std::cout << "  p_nu: [" << p_nu[0] << "," << p_nu[1] << "," << p_nu[2] << "]" << std::endl;
    std::cout << "  partial1=" << partial1 << std::endl;
    std::cout << "  p_dcm_nu: [" << p_dcm_nu[0] << "," << p_dcm_nu[1] << "," << p_dcm_nu[2] << "," << p_dcm_nu[3] << "]" << std::endl;
    std::cout << "  muparent_p: [" << muparent_px << "," << muparent_py << "," << muparent_pz << "]" << std::endl;
    std::cout << "  gammamp=" << gammamp << std::endl;
    std::cout << "  betamp: [" << betamp[0] << "," << betamp[1] << "," << betamp[2] << "]" << std::endl;
    std::cout << "  partial2=" << partial2 << std::endl;
    std::cout << "  p_pcm_mp: [" << p_pcm_mp[0] << "," << p_pcm_mp[1] << "," << p_pcm_mp[2] << "]  p_pcm=" << p_pcm << std::endl;
    std::cout << "  ntype=" << ntype 
              << " costhmu=" << costhmu << " wgt_ratio=" << wgt_ratio << std::endl;
  }
}

xypartials& xypartials::GetStaticInstance()
{ return gpartials; }
#endif
//___________________________________________________________________________

bool GNuMIFlux::LoadConfig(string cfg)
{
  genie::flux::GNuMIFluxXMLHelper helper(this);
  return helper.LoadConfig(cfg);
}

//___________________________________________________________________________

void GNuMIFlux::PrintConfig()
{
  
  std::ostringstream s;
  PDGCodeList::const_iterator itr = fPdgCList->begin();
  for ( ; itr != fPdgCList->end(); ++itr)
    s << (*itr) << " ";

  LOG("Flux", pNOTICE)
    << "GNuMIFlux Config:"
    << "\n Enu_max " << fMaxEv 
    << "\n pdg-codes: " << s.str()
    << "\n " << (fG3NuMI?"g3numi":"") << (fG4NuMI?"g4numi":"") << "("
    << fNuFluxTreeName << "), " << fNEntries << " entries " 
    << " (FilePOTs " << fFilePOTs << ") "
    <<  "in " << fNFiles << " files like: "
    << "\n " << fNuFluxFilePattern
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
    << "\n User Beam Origin: "
    << "\n  base " << utils::print::X4AsString(&fBeamZero)
    << "\n BeamRot/BeamRotInv ... not yet implemented"
    << "\n UseFluxAtDetCenter " << fUseFluxAtDetCenter;

}

//___________________________________________________________________________
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
  string fname = utils::xml::GetXMLFilePath("GNuMIFlux.xml");

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
      this->LoadParamSet(xml_doc,pval); // recurse
      SLOG("GNuMIFlux", pWARN) << "done using_param_set: \"" << pval << "\"";
    } else if ( pname == "units" ) {
      double scale = genie::utils::units::UnitFromString(pval);
      fGNuMI->SetLengthUnits(scale);
      SLOG("GNuMIFlux", pINFO) << "set user units to \"" << pval << "\"";

    } else if ( pname == "beamdir" ) {
      ParseBeamDir(xml_doc,xml_child);
      fGNuMI->SetBeamRotation(fBeamRot);

    } else if ( pname == "beampos" ) {
      ParseBeamPos(pval);
      fGNuMI->SetBeamCenter(fBeamPos);

    } else if ( pname == "window" ) {
      ParseWindowSeries(xml_doc,xml_child);
      // RWH  !!!! MEMORY LEAK!!!!
      //std::cout << " flux window " << std::endl
      //          << " [0] " << utils::print::X4AsString(new TLorentzVector(fFluxWindowPt[0],0)) << std::endl
      //          << " [1] " << utils::print::X4AsString(new TLorentzVector(fFluxWindowPt[1],0)) << std::endl
      //          << " [2] " << utils::print::X4AsString(new TLorentzVector(fFluxWindowPt[2],0)) << std::endl;

      fGNuMI->SetFluxWindow(fFluxWindowPt[0],fFluxWindowPt[1],fFluxWindowPt[2]);

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
  fBeamRot.SetToIdentity(); // start fresh

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
      TVector3 newX = AnglesToAxis(thetaphi3[0],thetaphi3[1],units);
      TVector3 newY = AnglesToAxis(thetaphi3[2],thetaphi3[3],units);
      TVector3 newZ = AnglesToAxis(thetaphi3[4],thetaphi3[5],units);
      fBeamRot.RotateAxes(newX,newY,newZ);
    } else {
      SLOG("GNuMIFlux", pWARN)
        << " type=\"" << dirtype << "\" within <beamdir> needs 6 values";
    }

  } else if ( dirtype == "newxyz" ) {
    // G4 style new axis values
    std::vector<double> newdir = GetDoubleVector(pval);
    if ( newdir.size() == 9 ) {
      TVector3 newX = TVector3(newdir[0],newdir[1],newdir[2]).Unit();
      TVector3 newY = TVector3(newdir[3],newdir[4],newdir[5]).Unit();
      TVector3 newZ = TVector3(newdir[6],newdir[7],newdir[8]).Unit();
      fBeamRot.RotateAxes(newX,newY,newZ);
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
    std::cout << " fBeamRot: " << std::setprecision(p) << std::endl;
    std::cout << " [ " 
              << std::setw(w) << fBeamRot.XX() << " "
              << std::setw(w) << fBeamRot.XY() << " "
              << std::setw(w) << fBeamRot.XZ() << endl
              << "   " 
              << std::setw(w) << fBeamRot.YX() << " "
              << std::setw(w) << fBeamRot.YY() << " "
              << std::setw(w) << fBeamRot.YZ() << endl
              << "   " 
              << std::setw(w) << fBeamRot.ZX() << " "
              << std::setw(w) << fBeamRot.ZY() << " "
              << std::setw(w) << fBeamRot.ZZ() << " ] " << std::endl;
    std::cout << std::endl;
  }

}

void GNuMIFluxXMLHelper::ParseBeamPos(std::string str)
{
  std::vector<double> xyz = GetDoubleVector(str);
  if ( xyz.size() == 3 ) {
    fBeamPos = TVector3(xyz[0],xyz[1],xyz[2]);
  } else if ( xyz.size() == 6 ) {
    // should check for '=' between triplets but we won't be so pedantic
    // ( userx, usery, userz ) = ( beamx, beamy, beamz )
    TVector3 userpos(xyz[0],xyz[1],xyz[2]);
    TVector3 beampos(xyz[3],xyz[4],xyz[5]);
    fBeamPos = userpos - fBeamRot*beampos;
  } else {
    SLOG("GNuMIFlux", pWARN)
      << "Unable to parse " << xyz.size() << " values in <beampos>";
    return;
   }
  if ( fVerbose > 1 ) {
    int w=16, p=10;
    std::cout << " fBeamPos: [ " << std::setprecision(p) 
              << std::setw(w) << fBeamPos.X() << " , "
              << std::setw(w) << fBeamPos.Y() << " , "
              << std::setw(w) << fBeamPos.Z() << " ] "
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
  fBeamRot = TRotation(); // reset matrix

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

      if      ( axis[0] == 'x' || axis[0] == 'X' ) fBeamRot.RotateX(rot);
      else if ( axis[0] == 'y' || axis[0] == 'Y' ) fBeamRot.RotateY(rot);
      else if ( axis[0] == 'z' || axis[0] == 'Z' ) fBeamRot.RotateZ(rot);
      else {
        SLOG("GNuMIFlux", pINFO)
          << " no " << axis << " to rotate around";
      }

    } else {
      SLOG("GNuMIFlux", pWARN)
        << " found <" << name << "> within <beamdir type=\"series\">";
    }
  }
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
        fFluxWindowPt[ientry] = pt;  // save the point
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
