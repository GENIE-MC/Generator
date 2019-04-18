//________________________________________________________________________________________
/*!

\program gevgen_upmu

\brief   A GENIE upgoing muon event generation application.

         *** Synopsis :

           gevgen_upmu [-h]
                       [-r run#]
                        -f flux
                        -n n_of_events
                        -d detector_bounding_box_size
                        -g rock_composition
                       [--seed random_number_seed]
                        --cross-sections xml_file
                       [--message-thresholds xml_file]

         *** Options :

           [] Denotes an optional argument.

           -h
              Prints out the syntax and exits.
           -r
              Specifies the MC run number
              [default: 1000]
           -f
              Specifies the input flux files
              The general syntax is: `-f simulation:/path/file.data[neutrino_code],...'
              [Notes]
               - The `simulation' string can be either `FLUKA' or `BGLRS' (so that
                 input data are binned using the correct FLUKA and BGLRS energy and
                 costheta binning). See comments in
                 - $GENIE/src/Flux/GFLUKAAtmoFlux.h
                 - $GENIE/src/Flux/GBGLRSAtmoFlux.h
                 and follow the links to the FLUKA and BGLRS atmo. flux web pages.
               - The neutrino codes are the PDG ones, for numu and anumu only
               - The /path/file.data,neutrino_code part of the option can be
                 repeated multiple times (separated by commas), once for each
                 flux neutrino species you want to consider,
                 eg. '-f FLUKA:~/data/sdave_numu07.dat[14],~/data/sdave_nue07.dat[12]'
                 eg. '-f BGLRS:~/data/flux10_271003_z.kam_nue[12]'
           -n
              Specifies how many events to generate.
           -d
              Specifies side length (in mm) of the detector bounding box.
              [default 100m (100000mm)]
           -g
              rock composition << NOT IMPLEMENTED YET >>
              Rock composition should be implemented in terms of materials defined in
              src/MuELoss/MuELMaterial.h and their weight fraction
              eg like -g 'silicon[0.30],calcium[0.29],iron[0.02],...'
           --seed
              Random number seed.
           --cross-sections
              Name (incl. full path) of an XML file with pre-computed
              cross-section values used for constructing splines.
           --message-thresholds
              Allows users to customize the message stream thresholds.
              The thresholds are specified using an XML file.
              See $GENIE/config/Messenger.xml for the XML schema.

         *** Examples:

           (1) Generate 100k events (run number 999210) for nu_mu only, using the
               sdave_numu07.dat FLUKA flux file (files in /data/flx/), and a detector of
               side length 50m.

               % gevgen_upmu -r 999210 -n 100000 -d 50000
                       -f FLUKA:/data/flx/sdave_numu07.dat[14]


         You can further control the GENIE behaviour by setting its standard
         environmental variables.
         Please read the GENIE User Manual for more information.

\created April 15, 2011

\author  Jen Truby <jen.truby \at sky.com>
         Oxford University - University of Liverpool & STFC Rutherford Appleton Lab summer student

         Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

\cpright Copyright (c) 2003-2019, The GENIE Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//_________________________________________________________________________________________

#include <cassert>
#include <cstdlib>
#include <cctype>
#include <string>
#include <vector>
#include <sstream>
#include <map>

#include <TRotation.h>
#include <TH1D.h>
#include <TH3D.h>
#include <TTree.h>
#include <TLorentzVector.h>

#include "Framework/Algorithm/Algorithm.h"
#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Units.h"
#include "Framework/EventGen/GFluxI.h"
#include "Framework/EventGen/GMCJDriver.h"
#include "Framework/EventGen/GMCJMonitor.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/Utils/XSecSplineList.h"
#include "Framework/Utils/StringUtils.h"
#include "Framework/Utils/SystemUtils.h"
#include "Framework/Utils/PrintUtils.h"
#include "Framework/Utils/UnitUtils.h"
#include "Framework/Utils/KineUtils.h"
#include "Framework/Utils/AppInit.h"
#include "Framework/Utils/RunOpt.h"
#include "Framework/Utils/CmdLnArgParser.h"

#ifdef __GENIE_FLUX_DRIVERS_ENABLED__
#include "Tools/Flux/GFLUKAAtmoFlux.h"
#include "Tools/Flux/GBGLRSAtmoFlux.h"
#endif

using std::string;
using std::vector;
using std::map;
using std::ostringstream;

using namespace genie;
using namespace genie::flux;
using namespace genie::constants;

void       GetCommandLineArgs     (int argc, char ** argv);
void       PrintSyntax            (void);
GFluxI *   GetFlux                (void);
void       GenerateUpNu           (GFluxI * flux_driver);
TH3D *     BuildEmuEnuCosThetaPdf (int nu_code);
TH1D *     GetEmuPdf              (double Enu, double costheta, const TH3D* pdf3d);
TVector3   GetDetectorVertex      (double CosTheta, double Enu);
double     GetCrossSection        (int nu_code, double Enu, double Emu);
double     ProbabilityEmu         (int nu_code, double Enu, double Emu);

// User-specified options:
//
Long_t          gOptRunNu;                     // run number
string          gOptFluxSim;                   // flux simulation (FLUKA or BGLRS)
map<int,string> gOptFluxFiles;                 // neutrino pdg code -> flux file map
int             gOptNev = -1;                  // exposure - in terms of number of events
double          gOptDetectorSide;              // detector side length, in mm.
long int        gOptRanSeed;                   // random number seed
string          gOptInpXSecFile;               // cross-section splines

// Constants
const double a = 2e+6;           // a = 2 MeV / (g cm-2)
const double e = 500e+9;         // e = 500 GeV

// Defaults:
//
double kDefOptDetectorSide = 1e+5;     // side length of detector, 100m (in mm)

//________________________________________________________________________________________
int main(int argc, char** argv)
{
  // Parse command line arguments
  GetCommandLineArgs(argc,argv);

  if ( ! RunOpt::Instance()->Tune() ) {
    LOG("gmkspl", pFATAL) << " No TuneId in RunOption";
    exit(-1);
  }
  RunOpt::Instance()->BuildTune();

  // Init
  utils::app_init::MesgThresholds(RunOpt::Instance()->MesgThresholdFiles());
  utils::app_init::RandGen(gOptRanSeed);
  utils::app_init::XSecTable(gOptInpXSecFile, true);


  // Get requested flux driver
  GFluxI * flux_driver = GetFlux();

  // Create output tree to store generated up-going muons
  TTree * ntupmuflux = new TTree("ntupmuflux","GENIE Upgoing Muon Event Tree");
  // Tree branches
  int    brIev        = 0;      // Event number
  int    brNuCode     = 0;      // Neutrino PDG code
  double brEmu        = 0;      // Muon energy (GeV)
  double brEnu        = 0;      // Neutrino energy (GeV)
  double brCosTheta   = 0;      // Neutrino cos(zenith angle), muon assumed to be collinear
  double brWghtFlxNu  = 0;      // Weight associated with the current flux neutrino
  double brWghtEmuPdf = 0;      // Weight associated with the Emu pdf for current Enu, costheta bins
  double brVx         = 0;      // Muon x (mm) - intersection with box surrounding the detector volume
  double brVy         = 0;      // Muon y (mm) - intersection with box surrounding the detector volume
  double brVz         = 0;      // Muon z (mm) - intersection with box surrounding the detector volume
  double brXSec       = 0;      //
  ntupmuflux->Branch("iev",          &brIev,         "iev/I"         );
  ntupmuflux->Branch("nu_code",      &brNuCode,      "nu_code/I"     );
  ntupmuflux->Branch("Emu",          &brEmu,         "Emu/D"         );
  ntupmuflux->Branch("Enu",          &brEnu,         "Enu/D"         );
  ntupmuflux->Branch("costheta",     &brCosTheta,    "costheta/D"    );
  ntupmuflux->Branch("wght_fluxnu",  &brWghtFlxNu,   "wght_fluxnu/D" );
  ntupmuflux->Branch("wght_emupdf",  &brWghtEmuPdf,  "wght_emupdf/D" );
  ntupmuflux->Branch("vx",           &brVx,          "vx/D"          );
  ntupmuflux->Branch("vy",           &brVy,          "vy/D"          );
  ntupmuflux->Branch("vz",           &brVz,          "vz/D"          );
  ntupmuflux->Branch("xsec",         &brXSec,        "xsec/D"        );

  // Build 3-D pdfs describing the the probability of a muon neutrino (or anti-neutrino)
  // of energy Enu and zenith angle costheta producing a mu- (or mu+) of energy E_mu
  TH3D * pdf3d_numu    = BuildEmuEnuCosThetaPdf (kPdgNuMu    );
  TH3D * pdf3d_numubar = BuildEmuEnuCosThetaPdf (kPdgAntiNuMu);

  // Up-going muon event loop
  for(brIev = 0; brIev < gOptNev; brIev++) {

    // Loop until upgoing nu is generated.
    GenerateUpNu(flux_driver);

    // Get neutrino code, Enu, costheta and weight for the generated neutrino
    brNuCode    = flux_driver->PdgCode();
    brEnu       = flux_driver->Momentum().E();
    brCosTheta  = -1. * flux_driver->Momentum().Pz() / flux_driver->Momentum().Vect().Mag();
    brWghtFlxNu = flux_driver->Weight();

    LOG("gevgen_upmu", pNOTICE)
        << "Generated flux neutrino: code = " << brNuCode
        << ", Ev = " << brEnu << " GeV"
        << ", cos(theta) = " << brCosTheta
        << ", weight = " << brWghtFlxNu;

    // Get Emu pdf my slicing the 3-D Enu,Emu,costheta pdf.
    TH3D * pdf3d  = (brNuCode == kPdgNuMu) ? pdf3d_numu : pdf3d_numubar;
    TH1D * pdfEmu = GetEmuPdf(brEnu,brCosTheta,pdf3d);

    // Get a random Emu from the pdf, and get the weight for that Emu.
    brEmu        = pdfEmu->GetRandom();
    brWghtEmuPdf = pdfEmu->Integral("width");
    LOG("gevgen_upmu", pNOTICE)
        << "Selected muon has energy Emu = " << brEmu
        << " and Emu pdg weight = " << brWghtEmuPdf;

    // Randomly select the point at which the neutrino crosses the detector.
    // assumes a simple cube detector of side length gOptDetectorSide.
    TVector3 Vertex = GetDetectorVertex(brCosTheta,brEnu);
    brVx = Vertex.X();
    brVy = Vertex.Y();
    brVz = Vertex.Z();

    LOG("gevgen_upmu", pNOTICE)
      << "Generated muon position: (" << brVx << ", " << brVy << ", " << brVz << ") m";

    // Save the cross section for the interaction generated.
    brXSec = GetCrossSection(brNuCode, brEnu, brEmu);

    // save all relevant values to the ntuple
    ntupmuflux->Fill();

    // Clean-up
    delete pdfEmu;
  }

  // Save the muon ntuple and calculate 3-D pdfs

  ostringstream outfname;
  outfname << "./genie." << gOptRunNu << ".upmu.root";
  TFile outf(outfname.str().c_str(), "recreate");
  ntupmuflux    -> Write("ntupmu");
  pdf3d_numu    -> Write("pdf3d_numu");
  pdf3d_numubar -> Write("pdf3d_numubar");
  outf.Close();

  // Clean-up
  delete flux_driver;

  return 0;
}
//________________________________________________________________________________________
TVector3 GetDetectorVertex(double costheta, double Enu)
{
// Find the point at which the muon crosses the detector. Returns 0,0,0 if the muon misses.
// Detector is a cube of side length l (=gOptDetectorSide).

  // Get a RandomGen instance
  RandomGen * rnd = RandomGen::Instance();

  // Generate random phi.
  double phi = 2.*kPi* rnd->RndFlux().Rndm();

  // Set distance at which incoming muon is considered
  double rad = 0.87*gOptDetectorSide;

  // Set transverse radius of a circle
  double rad_trans = 0.87*gOptDetectorSide;

  // Get necessary trig
  double sintheta = TMath::Sqrt(1-costheta*costheta);
  double cosphi   = TMath::Cos(phi);
  double sinphi   = TMath::Sin(phi);

  // Compute the muon position at distance Rad.
  double z = rad * costheta;
  double y = rad * sintheta * cosphi;
  double x = rad * sintheta * sinphi;

  // Displace muon randomly on a circular surface of radius rad_trans,
  // perpendicular to a sphere radius rad at that position [x,y,z].
  TVector3 vec(x,y,z);                // vector towards selected point
  TVector3 dvec1 = vec.Orthogonal();  // orthogonal vector
  TVector3 dvec2 = dvec1;             // second orthogonal vector
  dvec2.Rotate(-kPi/2.0,vec);         // rotate second vector by 90deg -> Cartesian coords

  // Select a random point on the transverse surface, within radius rad_trans
  double psi = 2 * kPi * rnd->RndFlux().Rndm();      // rndm angle [0,2pi]
  double random = rnd->RndFlux().Rndm();             // rndm number [0,1]
  dvec1.SetMag(TMath::Sqrt(random)*rad_trans*TMath::Cos(psi));
  dvec2.SetMag(TMath::Sqrt(random)*rad_trans*TMath::Sin(psi));
  x += dvec1.X() + dvec2.X();    // x-coord of a point the muon passes through
  y += dvec1.Y() + dvec2.Y();    // y-coord of a point the muon passes through
  z += dvec1.Z() + dvec2.Z();    // z-coord of a point the muon passes through

  // Find out if the muon passes through any side of the detector.

  // Get momentum vector
  double pz = Enu * costheta;
  double py = Enu * sintheta * sinphi;
  double px = Enu * sintheta * cosphi;
  TVector3 p3(-px,-py,-pz);

  // Set up other vectors needed for line-box intersection.
  TVector3 x3(x,y,z);
  TVector3 temp3(x3);
  TVector3 Hit3(0,0,0);

  // Find out if the line of the muon intersects the detector.
  bool HitFound = false;
  double l = gOptDetectorSide;
  TVector3 unit(0,0,0);
  unit = p3.Unit();
  double unitx = unit.X();
  double unity = unit.Y();
  double unitz = unit.Z();

  // Check to see if muon intersects with z=-l/2 surface.
  if (x3.Z() < -l/2 && unitz > 0 && !HitFound){
    while (temp3.Z() < -l/2){
      temp3 += unit;
    }
    if ( -l/2<temp3.X() && l/2>temp3.X() && -l/2<temp3.Y() && l/2>temp3.Y() ){
      Hit3.SetXYZ(temp3.X(),temp3.Y(),-l/2);
      HitFound = true;
    }
    else temp3 = x3;
  }

  // Check to see if muon intersects with x=l/2 surface.
  if (x3.X() > l/2 && unitx < 0 && !HitFound){
    while (temp3.X() > l/2){
      temp3 += p3.Unit();
    }
    if ( -l/2<temp3.Y() && l/2>temp3.Y() && -l/2<temp3.Z() && l/2>temp3.Z() ){
      Hit3.SetXYZ(l/2,temp3.Y(),temp3.Z());
      HitFound = true;
    }
    else temp3 = x3;
  }

  // Check to see if muon intersects with y=l/2 surface.
  if (x3.Y() > l/2 && unity < 0 && !HitFound){
    while (temp3.Y() > l/2){
      temp3 += p3.Unit();
    }
    if ( -l/2<temp3.X() && l/2>temp3.X() && -l/2<temp3.Z() && l/2>temp3.Z() ){
      Hit3.SetXYZ(temp3.X(),l/2,temp3.Z());
      HitFound = true;
    }
    else temp3 = x3;
  }

  // Check to see if muon intersects with x=-l/2 surface.
  if (x3.X() < -l/2 && unitx > 0 && !HitFound){
    while (temp3.X() < -l/2){
      temp3 += p3.Unit();
    }
    if ( -l/2<temp3.Y() && l/2>temp3.Y() && -l/2<temp3.Z() && l/2>temp3.Z() ){
      Hit3.SetXYZ(-l/2,temp3.Y(),temp3.Z());
      HitFound = true;
    }
    else temp3 = x3;
  }

  // Check to see if muon intersects with y=-l/2 surface.
  if (x3.Y() < -l/2 && unity > 0 && !HitFound){
    while (temp3.Y() < -l/2){
      temp3 += p3.Unit();
    }
    if ( -l/2<temp3.X() && l/2>temp3.X() && -l/2<temp3.Z() && l/2>temp3.Z() ){
      Hit3.SetXYZ(temp3.X(),-l/2,temp3.Z());
      HitFound = true;
    }
    else temp3 = x3;
  }

  return Hit3;
}
//________________________________________________________________________________________
void GenerateUpNu(GFluxI * flux_driver)
{
// Generate a Neutrino. Keep generating neutrinos until an upgoing one is generated.

  while (1)
  {
    LOG("gevgen_upmu", pINFO) << "Pulling next neurtino from the flux driver...";
    flux_driver->GenerateNext();

    int nu_code = flux_driver->PdgCode();

    const TLorentzVector & p4 = flux_driver->Momentum();
    double costheta = -p4.Pz() / p4.Vect().Mag();

    bool keep = (costheta<0) && (nu_code==kPdgNuMu || nu_code==kPdgAntiNuMu);
    if(keep) return;
  }
}
//________________________________________________________________________________________
double ProbabilityEmu(int nu_code, double Enu, double Emu)
{
// Calculate the probability of an incoming neutrino of energy Enu
// generating a muon of energy Emu.
  double dxsec_dxdy = GetCrossSection(nu_code,Enu,Emu);
  double Int = e * constants::kNA * dxsec_dxdy / (a * (1 + Emu/e));
  return Int;
}
//________________________________________________________________________________________
TH3D* BuildEmuEnuCosThetaPdf(int nu_code)
{
// Set up a 3D histogram, with axes Emu, Enu, CosTheta.
// Bin convention is defined at the start.
// Content of each bin is given by the probability of getting a muon
// of energy Emu from a neutrino of energy Enu and zenith ange costheta

  // Bin convention for Enu, consistent with BGLRS
  const double Enumin            = 0.1;
  const int    nEnubinsPerDecade = 10;
  const int    nEnubins          = 31;

  // Bin convention for Emu, consistent with BGLRS
  const double Emumin            = 0.1;
  const int    nEmubinsPerDecade = 10;
  const int    nEmubins          = 31;

  // Bin convention for CosTheta, consistent with BGLRS
  const double costheta_min = -1;
  const double costheta_max = +1;
  const int    ncostheta    = 20;

  // Set up an array of CosTheta bins.
  double costhetabinwidth = ((costheta_max-costheta_min)/ncostheta);
  double CosThetaBinEdges[ncostheta+1];
  for (int i=0; i<=ncostheta; i++)
    {
      CosThetaBinEdges[i] = costheta_min + i*costhetabinwidth;
    }

  // Set up an array of Enu bins.
  Double_t MinLogEnu = log(Enumin);
  Double_t MaxLogEnu = log(10*Enumin);
  Double_t LogBinWidthEnu = ((MaxLogEnu-MinLogEnu)/nEnubinsPerDecade);
  Double_t EnuBinEdges[nEnubins+1];
  for (int i=0; i<=nEnubins; i++)
    {
      EnuBinEdges[i] = exp(MinLogEnu + i*LogBinWidthEnu);
    }

  // Set up an array of Emu bins.
  Double_t MinLogEmu = log(Emumin);
  Double_t MaxLogEmu = log(10*Emumin);
  Double_t LogBinWidthEmu = ((MaxLogEmu-MinLogEmu)/nEmubinsPerDecade);
  Double_t EmuBinEdges[nEmubins+2];
  for (int i=0; i<=nEmubins+1; i++)
    {
      EmuBinEdges[i] = exp(MinLogEmu + i*LogBinWidthEmu);
    }

  // Create 3D histogram. X-axis: Emu; Y-axis: Enu; Z-axis: CosTheta.
  TH3D *h3 = new TH3D("h3","",nEmubins+1,EmuBinEdges,nEnubins,EnuBinEdges,ncostheta,CosThetaBinEdges);

  // Draw histogram.
  //h3->Draw();

  // Calculate Probability for each triplet.
  for (int i=1; i< h3->GetNbinsX(); i++)
    {
      double Emu = h3->GetXaxis()->GetBinCenter(i);
      for (int j=1; j<= h3->GetNbinsY(); j++)
        {
          double Enu = h3->GetBinCenter(j);
          double Int = ProbabilityEmu(nu_code,Enu,Emu);
          for (int k=1; k<= h3->GetNbinsZ(); k++)
            {
              h3->SetBinContent(i,j,k,Int); // Bin Content will is Int.
              //h3->SetBinContent(i,j,k,1); // Bin content constant (to test)
            }
        }
    }
  return h3;
}
//________________________________________________________________________________________
TH1D * GetEmuPdf(double Enu, double costheta, const TH3D* pdf3d)
{
  // Build 1-D Emu pdf
  int Emu_nbins = pdf3d->GetXaxis()->GetNbins();
  const double * Emu_bins = pdf3d->GetXaxis()->GetXbins()->GetArray();
  TH1D * pdf_Emu = new TH1D("pdf_Emu","",Emu_nbins,Emu_bins);
  pdf_Emu->SetDirectory(0);

  // Find appropriate bins for Enu and CosTheta.
  int Enu_bin      = pdf3d->GetYaxis()->FindBin(Enu);
  int costheta_bin = pdf3d->GetZaxis()->FindBin(costheta);

  // For those Enu and costheta bins, build Emu pdf
  for (int Emu_bin = 1; Emu_bin < pdf_Emu->GetNbinsX(); Emu_bin++) {
    pdf_Emu->SetBinContent(Emu_bin, pdf3d->GetBinContent(Emu_bin,Enu_bin,costheta_bin));
  }

  return pdf_Emu;
}
//________________________________________________________________________________________
double GetCrossSection(int nu_code, double Enu, double Emu)
{
// Get the interaction cross section for a neutrino of energy Enu
// to produce a muon of energy Emu.

  double dxsec_dxdy = 0;

  if ( Emu >= Enu ) {
    dxsec_dxdy = 0;
  }
  else {
    // get the algorithm factory & config pool
    AlgFactory * algf = AlgFactory::Instance();

    // instantiate some algorithms
    AlgId id("genie::QPMDISPXSec","Default");
    Algorithm * alg = algf->AdoptAlgorithm(id);
    XSecAlgorithmI * xsecalg = dynamic_cast<XSecAlgorithmI*>(alg);

    Interaction * vp = Interaction::DISCC(kPdgTgtFreeP, kPdgProton,  nu_code, Enu);
    Interaction * vn = Interaction::DISCC(kPdgTgtFreeN, kPdgNeutron, nu_code, Enu);

    // Integrate over x [0,1] and y [0,1-Emu/Enu], 100 steps for each.
    double ymax = 1 - Emu/Enu;

    double dy   = ymax/10;
    double dx   = 0.1;
    double x    = 0;
    double y    = 0;

    for (int i = 0; i<=10; i++){
      x = i*dx;
      for (int j = 0; j<=10; j++){
        y = j * dy;
        double W  = 0;
        double Q2 = 0;
        utils::kinematics::XYtoWQ2(Enu,kNucleonMass,W,Q2,x,y);
        vp->KinePtr()->Setx(x);
        vp->KinePtr()->Sety(y);
        vp->KinePtr()->SetQ2(Q2);
        vp->KinePtr()->SetW(W);
        vn->KinePtr()->Setx(x);
        vn->KinePtr()->Sety(y);
        vn->KinePtr()->SetQ2(Q2);
        vn->KinePtr()->SetW(W);

        dxsec_dxdy += 0.5*(xsecalg->XSec(vp,kPSxyfE) + xsecalg->XSec(vn,kPSxyfE)) / units::cm2;
        }
    }

    delete vp;
    delete vn;
  }

  return dxsec_dxdy;
}
//________________________________________________________________________________________
GFluxI* GetFlux(void)
{
  GFluxI * flux_driver = 0;

#ifdef __GENIE_FLUX_DRIVERS_ENABLED__

  // Instantiate appropriate concrete flux driver
  GAtmoFlux * atmo_flux_driver = 0;
  if(gOptFluxSim == "FLUKA") {
     GFLUKAAtmoFlux * fluka_flux = new GFLUKAAtmoFlux;
     atmo_flux_driver = dynamic_cast<GAtmoFlux *>(fluka_flux);
  } else
  if(gOptFluxSim == "BGLRS") {
     GBGLRSAtmoFlux * bartol_flux = new GBGLRSAtmoFlux;
     atmo_flux_driver = dynamic_cast<GAtmoFlux *>(bartol_flux);
  } else {
     LOG("gevgen_upmu", pFATAL) << "Uknonwn flux simulation: " << gOptFluxSim;
     gAbortingInErr = true;
     exit(1);
  }


  atmo_flux_driver->GenerateWeighted(true);

  // Configure GAtmoFlux options (common to all concrete atmospheric flux drivers)
  // set flux files:
  map<int,string>::const_iterator file_iter = gOptFluxFiles.begin();
  for( ; file_iter != gOptFluxFiles.end(); ++file_iter) {
    int neutrino_code = file_iter->first;
    string filename   = file_iter->second;
    atmo_flux_driver->AddFluxFile(neutrino_code, filename);
  }
  atmo_flux_driver->LoadFluxData();
  // configure flux generation surface:
  atmo_flux_driver->SetRadii(1, 1);
  // Cast to GFluxI, the generic flux driver interface
  flux_driver = dynamic_cast<GFluxI *>(atmo_flux_driver);

#else
  LOG("gevgen_upmu", pFATAL) << "You need to enable the GENIE flux drivers first!";
  LOG("gevgen_upmu", pFATAL) << "Use --enable-flux-drivers at the configuration step.";
  gAbortingInErr = true;
  exit(1);
#endif

  return flux_driver;
}
//________________________________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{
// Get the command line arguments

  LOG("gevgen_upmu", pINFO) << "Parsing command line arguments";

  // Common run options. Set defaults and read.
  RunOpt::Instance()->ReadFromCommandLine(argc,argv);

  // Parse run options for this app

  CmdLnArgParser parser(argc,argv);

  // help?
  bool help = parser.OptionExists('h');
  if(help) {
    PrintSyntax();
    exit(0);
  }

  //
  // *** run number
  //
  if( parser.OptionExists('r') ) {
    LOG("gevgen_upmu", pDEBUG) << "Reading MC run number";
    gOptRunNu = parser.ArgAsLong('r');
  } else {
    LOG("gevgen_upmu", pDEBUG) << "Unspecified run number - Using default";
    gOptRunNu = 1000;
  } //-r

  //
  // *** exposure
  //

  // in number of events
  bool have_required_statistics = false;
  if( parser.OptionExists('n') ) {
    LOG("gevgen_upmu", pDEBUG)
        << "Reading number of events to generate";
    gOptNev = parser.ArgAsInt('n');
    have_required_statistics = true;
  }//-n
  if(!have_required_statistics) {
    LOG("gevgen_upmu", pFATAL)
       << "You must request exposure either in terms of number of events"
       << "\nUse the -n option";
    PrintSyntax();
    gAbortingInErr = true;
    exit(1);
  }

  //
  // *** detector side length
  //

  if( parser.OptionExists('d') ) {
    LOG("gevgen_upmu", pDEBUG)
      << "Reading detector side length";
    gOptDetectorSide = parser.ArgAsDouble('d');
  } else {
    LOG("gevgen_upmu", pDEBUG)
      << "Unspecified detector side length - Using default";
    gOptDetectorSide = kDefOptDetectorSide;
  }//-d

  //
  // *** flux files
  //

  // syntax:
  // simulation:/path/file.data[neutrino_code],/path/file.data[neutrino_code],...
  //
  if( parser.OptionExists('f') ) {
    LOG("gevgen_upmu", pDEBUG) << "Getting input flux files";
    string flux = parser.ArgAsString('f');

    // get flux simulation info (FLUKA,BGLRS) so as to instantiate the
    // appropriate flux driver
    string::size_type jsimend = flux.find_first_of(":",0);
    if(jsimend==string::npos) {
       LOG("gevgen_upmu", pFATAL)
           << "You need to specify the flux file source";
       PrintSyntax();
       gAbortingInErr = true;
       exit(1);
    }
    gOptFluxSim = flux.substr(0,jsimend);
    for(string::size_type i=0; i<gOptFluxSim.size(); i++) {
       gOptFluxSim[i] = toupper(gOptFluxSim[i]);
    }
    if((gOptFluxSim != "FLUKA") && (gOptFluxSim != "BGLRS")) {
        LOG("gevgen_upmu", pFATAL)
             << "The flux file source needs to be one of <FLUKA,BGLRS>";
        PrintSyntax();
        gAbortingInErr = true;
        exit(1);
    }
    // now get the list of input files and the corresponding neutrino codes.
    flux.erase(0,jsimend+1);
    vector<string> fluxv = utils::str::Split(flux,",");
    vector<string>::const_iterator fluxiter = fluxv.begin();
    for( ; fluxiter != fluxv.end(); ++fluxiter) {
       string filename_and_pdg = *fluxiter;
       string::size_type open_bracket  = filename_and_pdg.find("[");
       string::size_type close_bracket = filename_and_pdg.find("]");
       if (open_bracket ==string::npos ||
           close_bracket==string::npos)
       {
           LOG("gevgen_upmu", pFATAL)
              << "You made an error in specifying the flux info";
           PrintSyntax();
           gAbortingInErr = true;
           exit(1);
       }
       string::size_type ibeg = 0;
       string::size_type iend = open_bracket;
       string::size_type jbeg = open_bracket+1;
       string::size_type jend = close_bracket;
       string flux_filename   = filename_and_pdg.substr(ibeg,iend-ibeg);
       string neutrino_pdg    = filename_and_pdg.substr(jbeg,jend-jbeg);
       gOptFluxFiles.insert(
          map<int,string>::value_type(atoi(neutrino_pdg.c_str()), flux_filename));
    }
    if(gOptFluxFiles.size() == 0) {
       LOG("gevgen_upmu", pFATAL)
          << "You must specify at least one flux file!";
       PrintSyntax();
       gAbortingInErr = true;
       exit(1);
    }

  } else {
    LOG("gevgen_upmu", pFATAL) << "No flux info was specified! Use the -f option.";
    PrintSyntax();
    gAbortingInErr = true;
    exit(1);
  }

  //
  // *** random number seed
  //
  if( parser.OptionExists("seed") ) {
    LOG("gevgen_upmu", pINFO) << "Reading random number seed";
    gOptRanSeed = parser.ArgAsLong("seed");
  } else {
    LOG("gevgen_upmu", pINFO) << "Unspecified random number seed - Using default";
    gOptRanSeed = -1;
  }

  //
  // *** input cross-section file
  //
  if( parser.OptionExists("cross-sections") ) {
    LOG("gevgen_upmu", pINFO) << "Reading cross-section file";
    gOptInpXSecFile = parser.ArgAsString("cross-sections");
  } else {
    LOG("gevgen_upmu", pINFO) << "Unspecified cross-section file";
    gOptInpXSecFile = "";
  }


  //
  // print-out summary
  //

  PDGLibrary * pdglib = PDGLibrary::Instance();

  ostringstream fluxinfo;
  fluxinfo << "Using " << gOptFluxSim << " flux files: ";
  map<int,string>::const_iterator file_iter = gOptFluxFiles.begin();
  for( ; file_iter != gOptFluxFiles.end(); ++file_iter) {
     int neutrino_code = file_iter->first;
     string filename   = file_iter->second;
     TParticlePDG * p = pdglib->Find(neutrino_code);
     if(p) {
        string name = p->GetName();
        fluxinfo << "(" << name << ") -> " << filename << " / ";
     }
  }

  ostringstream expinfo;
  if(gOptNev > 0) { expinfo << gOptNev << " events";   }

  LOG("gevgen_atmo", pNOTICE)
     << "\n\n"
     << utils::print::PrintFramedMesg("gevgen_upmu job configuration");

  LOG("gevgen_upmu", pNOTICE)
     << "\n"
     << "\n @@ Run number: " << gOptRunNu
     << "\n @@ Random number seed: " << gOptRanSeed
     << "\n @@ Using cross-section file: " << gOptInpXSecFile
     << "\n @@ Flux"
     << "\n\t" << fluxinfo.str()
     << "\n @@ Exposure"
     << "\n\t" << expinfo.str()
     << "\n\n";
}
//________________________________________________________________________________________
void PrintSyntax(void)
{
  LOG("gevgen_upmu", pFATAL)
   << "\n **Syntax**"
   << "\n gevgen_upmu [-h]"
   << "\n             [-r run#]"
   << "\n              -f simulation:flux_file[neutrino_code],..."
   << "\n              -n n_of_events,"
   << "\n             [-d detector side length (mm)]"
   << "\n             [--seed random_number_seed]"
   << "\n              --cross-sections xml_file"
   << "\n            [--message-thresholds xml_file]"
   << "\n"
   << " Please also read the detailed documentation at http://www.genie-mc.org"
   << "\n";
}
//________________________________________________________________________________________
