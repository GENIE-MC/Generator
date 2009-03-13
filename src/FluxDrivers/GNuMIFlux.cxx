//____________________________________________________________________________
/*
 Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory

         Robert Hatcher <rhatcher@fnal.gov>
         Fermi National Accelerator Laboratory

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ June 27, 2008 - CA
   The first implementation of this concrete flux driver was first added in
   the development version 2.5.1
*/
//____________________________________________________________________________

#include <cstdlib>
#include <fstream>
#include <vector>

#include "libxml/xmlmemory.h"
#include "libxml/parser.h"

#include "Utils/XmlParserUtils.h"

#include <TFile.h>
#include <TChain.h>
#include <TSystem.h>

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

       // these should go in more general package
       std::vector<std::string> TokenizeString(std::string str, std::string sep);
       std::vector<double> GetDoubleVector(std::string str);
       std::string GetXMLPathList();
       std::string GetXMLFilePath(std::string basename);

    private:
      void     ParseParamSet(xmlDocPtr&, xmlNodePtr&);
      void     ParseBeamDir(xmlDocPtr&, xmlNodePtr&);
      void     ParseBeamPos(std::string);
      void     ParseRotSeries(xmlDocPtr&, xmlNodePtr&);
      void     ParseWindowSeries(xmlDocPtr&, xmlNodePtr&);
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
const TLorentzVector kPosCenterNearBeam(0.,0.,103935.,0.);
const TLorentzVector kPosCenterFarBeam(0.,0.,73534000.,0.);

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
  if (!fDetLocIsSet) {
     LOG("Flux", pERROR)
       << "Specify a detector location before generating flux neutrinos";
     return false;	
  }
  if (fMaxWeight<=0) {
     LOG("Flux", pERROR)
       << "Run ScanForMaxWeight() before generating unweighted flux neurtinos";
     return false;	
  }

  RandomGen* rnd = RandomGen::Instance();
  while ( true ) {
     // Check for end of flux ntuple
     bool end = this->End();
     if ( end ) return false;

     // Get next weighted flux ntuple entry
     bool nextok = this->GenerateNext_weighted();
     if ( fGenWeighted ) return nextok;
     if ( ! nextok ) continue;

     if ( fNCycles == 0 ) {
       LOG("Flux", pNOTICE)
          << "Got flux entry: " << fIEntry
          << " - Cycle: "<< fICycle << "/ infinite";
     } else {
       LOG("Flux", pNOTICE)
          << "Got flux entry: "<< fIEntry
          << " - Cycle: "<< fICycle << "/"<< fNCycles;
     }

     // Get fractional weight & decide whether to accept curr flux neutrino
     double f = this->Weight() / fMaxWeight;
     LOG("Flux", pNOTICE)
        << "Curr flux neutrino fractional weight = " << f;
     if (f > 1.) {
       LOG("Flux", pERROR)
           << "** Fractional weight = " << f << " > 1 !!";
       fMaxWeight = this->Weight(); // bump the weight
     }
     double r = (f < 1.) ? rnd->RndFlux().Rndm() : 0;
     bool accept = ( r < f );
     if ( accept ) {
       fWeight = 1.;
       return true;
     }

     LOG("Flux", pNOTICE)
       << "** Rejecting current flux neutrino based on the flux weight only";
  }
  return false;
}
//___________________________________________________________________________
bool GNuMIFlux::GenerateNext_weighted(void)
{
// Get next (weighted) flux ntuple entry on the specified detector location
//
  // Reset previously generated neutrino code / 4-p / 4-x
  this->ResetCurrent();

  // Check whether a flux ntuple has been loaded
  if ( ! fG3NuMI && ! fG4NuMI ) {
     LOG("Flux", pWARN)
          << "The flux driver has not been properly configured";
     return false;	
  }

  // Read next flux ntuple entry
  if (fIEntry >= fNEntries) {
     // Run out of entries @ the current cycle.
     // Check whether more (or infinite) number of cycles is requested
     if (fICycle < fNCycles || fNCycles == 0 ) {
        fICycle++;
        fIEntry=0;
     } else {
        LOG("Flux", pWARN)
            << "No more entries in input flux neutrino ntuple";
        return false;	
     }
  }

  //rwh//old: fNuFluxTree->GetEntry(fIEntry);
  if ( fG3NuMI ) { fG3NuMI->GetEntry(fIEntry); fCurrentEntry->Copy(fG3NuMI); }
  else           { fG4NuMI->GetEntry(fIEntry); fCurrentEntry->Copy(fG4NuMI); }

  fIEntry++;
  fCurrentEntry->pcodes = 0;  // fetched entry has geant codes
  fCurrentEntry->units  = 0;  // fetched entry has original units

  // Convert the current gnumi neutrino flavor mode into a neutrino pdg code
  // Also convert other particle codes in GNuMIFluxPassThroughInfo to PDG
  fCurrentEntry->ConvertPartCodes();
  fgPdgC = fCurrentEntry->ntype;

  // Check neutrino pdg against declared list of neutrino species declared
  // by the current instance of the NuMI neutrino flux driver.
  // No undeclared neutrino species will be accepted at this point as GENIE
  // has already been configured to handle the specified list.
  // Make sure that the appropriate list of flux neutrino species was set at
  // initialization via GNuMIFlux::SetFluxParticles(const PDGCodeList &)

  if ( ! fPdgCList->ExistsInPDGCodeList(fgPdgC) ) {
     LOG("Flux", pWARN)
          << "Unknown decay mode or decay mode producing an undeclared neutrino species";
     LOG("Flux", pWARN)
          << "Declared list of neutrino species: " << *fPdgCList;
     LOG("Flux", pWARN)
       << "gnumi neutrino flavor: " << fCurrentEntry->ntype 
       << "(pcodes=" << fCurrentEntry->pcodes << ")" ;
     return false;	
  }

  // Update the curr neutrino weight and energy

  // Check current neutrino energy against the maximum flux neutrino energy declared
  // by the current instance of the NuMI neutrino flux driver.
  // No flux neutrino exceeding that maximum energy will be accepted at this point as
  // that maximum energy has already been used for normalizing the interaction probabilities.
  // Make sure that the appropriate maximum flux neutrino energy was set at
  // initialization via GNuMIFlux::SetMaxEnergy(double Ev)
  fWeight = fCurrentEntry->nimpwt;   // start with importance weight
  double Ev = 0;
  fgP4 = fFluxWindowBase;
  switch ( fUseFluxAtDetCenter ) {
  case -1:  // near detector
    fWeight *= fCurrentEntry->nwtnear;
    Ev       = fCurrentEntry->nenergyn;
    break;
  case +1:  // far detector
    fWeight *= fCurrentEntry->nwtfar;
    Ev       = fCurrentEntry->nenergyf;
    break;
  default:  // recalculate on x-y window
    double wgt_xy = 0;
    RandomGen * rnd = RandomGen::Instance();
    fgP4 += ( rnd->RndFlux().Rndm()*fFluxWindowDir1 +
              rnd->RndFlux().Rndm()*fFluxWindowDir2   );
#ifdef  GNUMI_TEST_XY_WGT
    xypartials partials;
    fCurrentEntry->CalcEnuWgt(fgX4.X(),fgX4.Y(),fgX4.Z(),Ev,wgt_xy,partials);
#else
    fCurrentEntry->CalcEnuWgt(fgX4.X(),fgX4.Y(),fgX4.Z(),Ev,wgt_xy);
#endif
    fWeight *= wgt_xy;
    break;
  }

  if (Ev > fMaxEv) {
     LOG("Flux", pWARN)
          << "Flux neutrino energy exceeds declared maximum neutrino energy";
     LOG("Flux", pWARN)
          << "Ev = " << Ev << "(> Ev{max} = " << fMaxEv << ")";
  }

  // Set the current flux neutrino 4-momentum
  // RWH!!! currently in *beam* coordinates
  TLorentzVector dkVtx(fCurrentEntry->vx,fCurrentEntry->vy,fCurrentEntry->vz,0.);
  TLorentzVector dirNu = fgX4 - dkVtx;
  double dirnorm = 1.0 / dirNu.Vect().Mag();  // magnitude of position offset
  fgP4.SetPxPyPzE (Ev*dirnorm*dirNu.X(), Ev*dirnorm*dirNu.Y(), Ev*dirnorm*dirNu.Z(), Ev);

  // Set the current flux neutrino 4-position, direction in user coord
  Beam2UserDir(fgP4,fgP4User);
  Beam2UserPos(fgX4,fgX4User);

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("Flux", pINFO)
	<< "Generated neutrino: "
	<< "\n pdg-code: " << fgPdgC
        << "\n p4: " << utils::print::P4AsShortString(&fgP4)
        << "\n x4: " << utils::print::X4AsString(&fgX4);
#endif

  // update the sum of weights & number of neutrinos
  fSumWeight += this->Weight();
  fNNeutrinos++;

  return true;
}
//___________________________________________________________________________
double GNuMIFlux::POT_curr(void)
{
// Compute current number of flux POTs

  if (!fNuFluxTree) {
     LOG("Flux", pWARN)
          << "The flux driver has not been properly configured";
     return 0;	
  }

  return 0;
}
//___________________________________________________________________________
void GNuMIFlux::LoadBeamSimData(string filename)
{
// Loads in a gnumi beam simulation root file (converted from hbook format)
// into the GNuMIFlux driver.

  fNuFluxFilePattern = filename;
  LOG("Flux", pNOTICE)
        << "Loading gnumi flux tree from ROOT file(s): " << filename;

  // !WILDCARD only works for file name ... NOT directory
  string dirname = gSystem->UnixPathName(gSystem->WorkingDirectory());
  size_t slashpos = filename.find_last_of("/");
  if ( slashpos != std::string::npos ) {
    dirname = filename.substr(0,slashpos);
    LOG("Flux", pINFO) << "Look for flux using directory " << dirname;
  } else { slashpos = -1; }

  void* dirp = gSystem->OpenDirectory(gSystem->ExpandPathName(dirname.c_str()));
    // create a (sortable) vector of file names
  std::vector<std::string> fnames;
  if ( dirp ) {
    std::string basename = 
      filename.substr(slashpos+1,filename.size()-slashpos-1);
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
            LOG("Flux",pNOTICE) //INFO)
              << fNuFluxTreeName << "->AddFile() of "
              << fnames[indx];
            fNuFluxTree->AddFile(fnames[indx].c_str());
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

  // this will open all files and head header!!
  fNEntries = fNuFluxTree->GetEntries();

  LOG("Flux", pNOTICE)
      << "Loaded flux tree contains " <<  fNEntries << " entries";

  // current ntuple cycle # (flux ntuples may be recycled)
  fICycle = 1;
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
    if ( ipos_estimator == 0 )
      LOG("Flux", pERROR)
        << "Don't know how to estimate max weight for the flux window base Z=" 
        << zbase << " position";
  }
  if ( ipos_estimator > 0 ) {
    if ( fG3NuMI ) fNuFluxTree->Draw("Nimpwt*Nwtfar","","goff");
    else           fNuFluxTree->Draw("Nimpwt*NWtFar[0]","","goff");
  }
  if ( ipos_estimator < 0 ) {
    if ( fG3NuMI ) fNuFluxTree->Draw("Nimpwt*Nwtnear","","goff");
    else           fNuFluxTree->Draw("Nimpwt*Nwtnear[0]","","goff");
  }
  Long64_t idx = TMath::LocMax(
                     fNuFluxTree->GetSelectedRows(),
                     fNuFluxTree->GetV1());
  fMaxWeight = fNuFluxTree->GetV1()[idx];
  LOG("Flux", pNOTICE) << "Maximum flux weight = " << fMaxWeight;
  if ( fMaxWeight <= 0 ) {
      LOG("Flux", pFATAL) << "Non-positive maximum flux weight!";
      exit(1);
  }
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
void GNuMIFlux::SetFilePOT(double pot)
{
// POTs in input flux file

  fFilePOT = pot;
}
//___________________________________________________________________________
void GNuMIFlux::SetUpstreamZ(double z0)
{
// The flux neutrino position (x,y) is given at the detector coord system
// at z=0. This method sets the preferred starting z position upstream of
// the upstream detector face. Each flux neutrino will be backtracked from
// z=0 to the input z0.

  fZ0 = z0;
}
//___________________________________________________________________________
void GNuMIFlux::SetNumOfCycles(int n)
{
// The flux ntuples can be recycled for a number of times to boost generated
// event statistics without requiring enormous beam simulation statistics.
// That option determines how many times the driver is going to cycle through
// the input flux ntuple.
// With n=0 the flux ntuple will be recycled an infinite amount of times so
// that the event generation loop can exit only on a POT or event num check.

  fNCycles = TMath::Max(0, n);
}
//___________________________________________________________________________
void GNuMIFlux::SetTreeName(string name)
{
  fNuFluxTreeName = name;
}
//___________________________________________________________________________
void GNuMIFlux::UseFluxAtNearDetCenter(void)
{
  fFluxWindowBase      = kPosCenterNearBeam;
  fFluxWindowDir1      = TLorentzVector();   // no extent
  fFluxWindowDir2      = TLorentzVector();
  fUseFluxAtDetCenter  = -1;
  fDetLocIsSet         = true;
}
//___________________________________________________________________________
void GNuMIFlux::UseFluxAtFarDetCenter(void)
{
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
  beam0 *= (1./fLengthScaleB2U);
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
  // NOTE: internally these are in "cm", but user might have set a preference
  TVector3 beam0 = fBeamZero.Vect();
  beam0 *= fLengthScaleB2U;
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
  usrxyz = fBeamRot*beamxyz + fBeamZero;
}
void GNuMIFlux::Beam2UserDir(const TLorentzVector& beamdir, 
                                   TLorentzVector& usrdir) const
{
  usrdir = fBeamRot*beamdir;
}
void GNuMIFlux::User2BeamPos(const TLorentzVector& usrxyz,
                                   TLorentzVector& beamxyz) const
{
  beamxyz = fBeamRotInv*(usrxyz-fBeamZero);
}
void GNuMIFlux::User2BeamDir(const TLorentzVector& usrdir,
                                   TLorentzVector& beamdir) const
{
  beamdir = fBeamRotInv*usrdir;
}

//___________________________________________________________________________
void GNuMIFlux::PrintCurrent(void)
{
  LOG("Flux", pNOTICE) << "CurrentEntry:" << *fCurrentEntry;
}
//___________________________________________________________________________
void GNuMIFlux::Initialize(void)
{
  LOG("Flux", pNOTICE) << "Initializing GNuMIFlux driver";

  fMaxEv           = 0;
  fPdgCList        = new PDGCodeList;
  fCurrentEntry    = new GNuMIFluxPassThroughInfo;

  fNuFluxTree      = 0;
  fG3NuMI          = 0;
  fG4NuMI          = 0;
  fNuFluxTreeName  = "";

  fNEntries        = 0;
  fIEntry          = 0;
  fMaxWeight       =-1;
  fFilePOT         = 0;
  fZ0              = 0;
  fNCycles         = 0;
  fICycle          = 0;
  fSumWeight       = 0;
  fNNeutrinos      = 0;
  fGenWeighted     = false;
  fUseFluxAtDetCenter = 0;
  fDetLocIsSet        = false;
  // by default assume user length is cm
  SetLengthUnits(genie::utils::units::UnitFromString("cm"));

  this->SetDefaults();
  this->ResetCurrent();
}
//___________________________________________________________________________
void GNuMIFlux::SetDefaults(void)
{
// - Set default neutrino species list (nue, nuebar, numu, numubar) and
//   maximum energy (125 GeV).
//   These defaults can be overwritten by user calls (at the driver init) to
//   GNuMIlux::SetMaxEnergy(double Ev) and
//   GNuMIFlux::SetFluxParticles(const PDGCodeList & particles)
// - Set the default file normalization to 1E+21 POT
// - Set the default flux neutrino start z position at -5m (z=0 is the
//   detector centre).
// - Set number of cycles to 1

  LOG("Flux", pNOTICE) << "Setting default GJPARCNuFlux driver options";

  PDGCodeList particles;
  particles.push_back(kPdgNuMu);
  particles.push_back(kPdgAntiNuMu);
  particles.push_back(kPdgNuE);
  particles.push_back(kPdgAntiNuE);

  this->SetFluxParticles (particles);
  this->SetMaxEnergy     (125./*GeV*/);  // was 200, but that would be wasteful
  this->SetFilePOT       (1E+21);
  this->SetUpstreamZ     (-5.0);
  this->SetNumOfCycles   (1);
}
//___________________________________________________________________________
void GNuMIFlux::ResetCurrent(void)
{
// reset running values of neutrino pdg-code, 4-position & 4-momentum
// and the input ntuple leaves

  fgPdgC  = 0;
  fWeight = 0;
  fgP4.SetPxPyPzE (0.,0.,0.,0.);
  fgX4.SetXYZT    (0.,0.,0.,0.);

  fCurrentEntry->Reset();
}
//___________________________________________________________________________
void GNuMIFlux::CleanUp(void)
{
  LOG("Flux", pNOTICE) << "Cleaning up...";

  if (fPdgCList)        delete fPdgCList;
  if (fCurrentEntry)    delete fCurrentEntry;

  if ( fG3NuMI ) delete fG3NuMI;
  if ( fG4NuMI ) delete fG4NuMI;

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

  double cm = genie::utils::units::UnitFromString("cm");
  fLengthScaleB2U = user_units / cm ;

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
{ Reset(); }
//___________________________________________________________________________
void GNuMIFluxPassThroughInfo::Reset()
{
  pcodes = -1;
  units  = -1;
  
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
void GNuMIFluxPassThroughInfo::Copy(const g3numi* g3 )
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
void GNuMIFluxPassThroughInfo::Copy(const g4numi* g4 )
{

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
int GNuMIFluxPassThroughInfo::CalcEnuWgt(double xpos, double ypos, double zpos,
                                         double& enu, double& wgt_xy
#ifdef  GNUMI_TEST_XY_WGT
                                         , xypartials& partials
#endif
                                         ) const
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

  const int kpdg_nue       =   12;  // extended Geant 53
  const int kpdg_nuebar    =  -12;  // extended Geant 52
  const int kpdg_numu      =   14;  // extended Geant 56
  const int kpdg_numubar   =  -14;  // extended Geant 55

  const int kpdg_muplus    =  -13;  // Geant  5
  const int kpdg_muminus   =   13;  // Geant  6
  const int kpdg_pionplus  =  211;  // Geant  8
  const int kpdg_pionminus = -211;  // Geant  9
  const int kpdg_k0long    =  130;  // Geant 10  ( K0=311, K0S=310 )
  const int kpdg_k0short   =  310;  // Geant 16
  const int kpdg_k0mix     =  311;  
  const int kpdg_kaonplus  =  321;  // Geant 11
  const int kpdg_kaonminus = -321;  // Geant 12

  const double kRDET = 100.0;   // set to flux per 100 cm radius

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
  default:
    std::cerr << "NU_REWGT unknown particle type " << this->ptype
              << std::endl;
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

  double costh_pardet = ( this->pdpx*(xpos-this->vx) +
                          this->pdpy*(ypos-this->vy) +
                          this->pdpz*(zpos-this->vz) ) 
                        / ( parentp * rad);
  if ( costh_pardet >  1.0 ) costh_pardet =  1.0;
  if ( costh_pardet < -1.0 ) costh_pardet = -1.0;
  double theta_pardet = TMath::ACos(costh_pardet);

  // Weighted neutrino energy in beam, approx, good for small theta
  //RWH//double emrat = 1.0 / ( gamma * ( 1.0 - beta_mag * TMath::Cos(theta_pardet)));
  double emrat = 1.0 / ( gamma * ( 1.0 - beta_mag * costh_pardet ));
  enu = emrat * enuzr;  // ! the energy ... normally

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
  partials.parent_mass   = parent_mass;
  partials.parentp       = parentp;
  partials.parent_energy = parent_energy;
  partials.gamma         = gamma;
  partials.beta_mag      = beta_mag;
  partials.enuzr         = enuzr;
  partials.rad           = rad;
  partials.costh_pardet  = costh_pardet;
  partials.theta_pardet  = theta_pardet;
  partials.emrat         = emrat;
  partials.eneu          = enu;
#endif

  // Get solid angle/4pi for detector element
  double sangdet = ( kRDET*kRDET / ( (zpos-this->vz)*(zpos-this->vz)))/4.0;

  // Weight for solid angle and lorentz boost
  wgt_xy = sangdet * ( emrat * emrat );  // ! the weight ... normally

#ifdef  GNUMI_TEST_XY_WGT
  partials.sangdet       = sangdet;
  partials.wgt           = wgt_xy;
  partials.ptype         = this->ptype; // assume already PDG
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
    partials.betanu[0]     = beta[0];
    partials.betanu[1]     = beta[1];
    partials.betanu[2]     = beta[2];
    partials.p_nu[0]       = p_nu[0];
    partials.p_nu[1]       = p_nu[1];
    partials.p_nu[2]       = p_nu[2];
    partials.partial1      = partial;
    partials.p_dcm_nu[0]   = p_dcm_nu[0];
    partials.p_dcm_nu[1]   = p_dcm_nu[1];
    partials.p_dcm_nu[2]   = p_dcm_nu[2];
    partials.p_dcm_nu[3]   = p_dcm_nu[3];
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
    partials.muparent_px   = this->muparpx;
    partials.muparent_py   = this->muparpy;
    partials.muparent_pz   = this->muparpz;
    partials.gammamp       = gamma;
    partials.betamp[0]     = beta[0];
    partials.betamp[1]     = beta[1];
    partials.betamp[2]     = beta[2];
    partials.partial2      = partial;
    partials.p_pcm_mp[0]   = p_pcm_mp[0];
    partials.p_pcm_mp[1]   = p_pcm_mp[1];
    partials.p_pcm_mp[2]   = p_pcm_mp[2];
    partials.p_pcm         = p_pcm;
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
    partials.ntype     = this->ntype; // assume converted to PDG
    partials.costhmu   = costh;
    partials.wgt_ratio = wgt_ratio;
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

void xypartials::Print() const
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
#endif
//___________________________________________________________________________

bool GNuMIFlux::LoadConfig(string cfg)
{
  genie::flux::GNuMIFluxXMLHelper helper(this);
  return helper.LoadConfig(cfg);
}

//___________________________________________________________________________
//___________________________________________________________________________

std::vector<std::string> 
GNuMIFluxXMLHelper::TokenizeString(std::string str, std::string sep)
{
  // Separate "str" string into elements under the assumption
  // that they are separated by any of the characters in "sep"
  std::vector<std::string> rlist;
  size_t pos_beg = 0;
  size_t str_end = str.size();
  while ( pos_beg != string::npos && pos_beg < str_end ) {
    size_t pos_end = str.find_first_of(sep.c_str(),pos_beg);
    std::string oneval = str.substr(pos_beg,pos_end-pos_beg);
    pos_beg = pos_end+1;
    if ( pos_end == string::npos ) pos_beg = str_end;
    if ( oneval != "" ) rlist.push_back(oneval);
  }
  return rlist;
}

std::string GNuMIFluxXMLHelper::GetXMLPathList()
{
  // get a separated list of potential locations for xml files
  // e.g. ".:$MYSITEXML:/path/to/exp/version:$GALGCONF:$GENIE"
  string pathlist; 
  const char* p = gSystem->Getenv("GXMLPATHS");
  if ( p ) { pathlist = std::string(p) + ":"; }
  // alternative path current supported
  p = gSystem->Getenv("GALGCONF");
  if ( p ) { pathlist = std::string(p) + ":"; }
  pathlist += "$GENIE/config";  // standard path in case no env
  return pathlist;
}

std::string GNuMIFluxXMLHelper::GetXMLFilePath(std::string basename)
{
  // return a full path to a real XML file
  // e.g. passing in "GNuMIFlux.xml"
  //   will return   "/Users/rhatcher/Software/GENIE/HEAD/config/GNuMIFlux.xml"
  // allow ::colon:: ::semicolon:: and ::comma:: as separators
  std::string pathlist = GetXMLPathList();
  std::vector<std::string> paths = TokenizeString(pathlist,":;,");
  // expand any wildcards, etc.
  size_t np = paths.size();
  for ( size_t i=0; i< np; ++i ) {
    const char* tmppath = paths[i].c_str();
    std::string onepath = gSystem->ExpandPathName(tmppath);
    onepath += "/";
    onepath += basename;
    bool noAccess = gSystem->AccessPathName(onepath.c_str());
    if ( ! noAccess ) return onepath;  // found one
  }
  // didn't find it, return basename in case it is in "." and that
  // wasn't listed in the XML path list.   If you want "." to take
  // precedence then it needs to be explicitly listed in GXMLPATHS.
  return basename;  

}

std::vector<double> GNuMIFluxXMLHelper::GetDoubleVector(std::string str)
{
  // turn string into vector<double>
  // be liberal about separators, users might punctuate for clarity
  std::vector<std::string> strtokens = TokenizeString(str," ,;:()[]=");
  std::vector<double> vect;
  size_t ntok = strtokens.size();

  if ( fVerbose > 2 ) 
    std::cout << "GetDoubleVector \"" << str << "\"" << std::endl;

  for (size_t i=0; i < ntok; ++i) {
    std::string trimmed = utils::str::TrimSpaces(strtokens[i]);
    if ( " " == trimmed || "" == trimmed ) continue;  // skip empty strings
    double val = atof(trimmed.c_str());
    if ( fVerbose > 2 ) 
      std::cout << "(" << vect.size() << ") = " << val << std::endl;
    vect.push_back(val);
  }

  return vect;
}

bool GNuMIFluxXMLHelper::LoadConfig(string cfg)
{
  //   (search at $GALGCONF or use the default at: $GENIE/config)
  //string fname = (gSystem->Getenv("GALGCONF")) ?
  //  string(gSystem->Getenv("GALGCONF")) + string("/") : 
  //  string(gSystem->Getenv("GENIE")) + string("/config/");
  //fname += "GNuMIFlux.xml";

  string fname = GetXMLFilePath("GNuMIFlux.xml");

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


  // loop over all xml tree nodes that are children of the root node
  // read the entries looking for "param_set" of the right name

  SLOG("GNuMIFlux", pINFO) << "Attempt to load config \"" << cfg 
                           << "\" from file: " << fname;

  // loop
  xmlNodePtr xml_pset = xml_root->xmlChildrenNode;
  for ( ; xml_pset != NULL ; xml_pset = xml_pset->next ) {
    if ( ! xmlStrEqual(xml_pset->name, (const xmlChar*)"param_set") ) continue;
    // every time there is a 'param_set' tag
    string param_set_name = 
      utils::str::TrimSpaces(XmlParserUtils::GetAttribute(xml_pset,"name"));
    
    if ( param_set_name != cfg ) continue;
      
    SLOG("GNuMIFlux", pINFO) << "Found config \"" << cfg << "\" in file:"
                               << fname;
    this->ParseParamSet(xml_doc,xml_pset);
  } // loop over elements of root
  xmlFree(xml_pset);
  xmlFree(xml_doc);
  return true;

}

void GNuMIFluxXMLHelper::ParseParamSet(xmlDocPtr& xml_doc, xmlNodePtr& xml_pset)
{
  xmlNodePtr xml_child = xml_pset->xmlChildrenNode;
  for ( ; xml_child != NULL ; xml_child = xml_child->next ) {
    // handle basic gnumi_config/param_set
    // bad cast away const on next line, but function sig requires it
    string pname = 
      XmlParserUtils::TrimSpaces(const_cast<xmlChar*>(xml_child->name));
    if ( pname == "text" || pname == "comment" ) continue;
    string pval  = 
      XmlParserUtils::TrimSpaces(
              xmlNodeListGetString(xml_doc, xml_child->xmlChildrenNode, 1));

    if ( fVerbose > 1 ) 
      SLOG("GNuMIFlux", pINFO)
        << "   pname \"" << pname << "\", string value \"" << pval << "\"";

    if        ( pname == "verbose" ) {
      fVerbose = atoi(pval.c_str());
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
      fGNuMI->SetFluxWindow(fFluxWindowPt[0],fFluxWindowPt[1],fFluxWindowPt[2]);

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
      XmlParserUtils::GetAttribute(xml_beamdir,"type"));

  string pval  = 
    XmlParserUtils::TrimSpaces(
      xmlNodeListGetString(xml_doc, xml_beamdir->xmlChildrenNode, 1));

  if        ( dirtype == "series" ) {
    // series of rotations around an axis
    ParseRotSeries(xml_doc,xml_beamdir);

  } else if ( dirtype == "thetaphi3") {
    // G3 style triplet of (theta,phi) pairs
    std::vector<double> thetaphi3 = GetDoubleVector(pval);
    string units = 
      utils::str::TrimSpaces(XmlParserUtils::GetAttribute(xml_beamdir,"units"));
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

void GNuMIFluxXMLHelper::ParseRotSeries(xmlDocPtr& xml_doc, xmlNodePtr& xml_pset)
{
  fBeamRot = TRotation(); // reset matrix

  xmlNodePtr xml_child = xml_pset->xmlChildrenNode;
  for ( ; xml_child != NULL ; xml_child = xml_child->next ) {
    // in a <beamdir> of type "series"
    // should be a sequence of <rotation> entries
    string name = 
      XmlParserUtils::TrimSpaces(const_cast<xmlChar*>(xml_child->name));
    if ( name == "text" || name == "comment" ) continue;

    if ( name == "rotation" ) {
      string val  = 
        XmlParserUtils::TrimSpaces(
          xmlNodeListGetString(xml_doc, xml_child->xmlChildrenNode, 1));
      string axis = 
        utils::str::TrimSpaces(XmlParserUtils::GetAttribute(xml_child,"axis"));

      string units = 
        utils::str::TrimSpaces(XmlParserUtils::GetAttribute(xml_child,"units"));

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
      XmlParserUtils::TrimSpaces(const_cast<xmlChar*>(xml_child->name));
    if ( name == "text" || name == "comment" ) continue;

    if ( name == "point" ) {
      string val  = 
        XmlParserUtils::TrimSpaces(
          xmlNodeListGetString(xml_doc, xml_child->xmlChildrenNode, 1));
      string coord = 
        utils::str::TrimSpaces(XmlParserUtils::GetAttribute(xml_child,"coord"));

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
