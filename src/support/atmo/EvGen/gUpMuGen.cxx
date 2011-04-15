//________________________________________________________________________________________
/*!

\program gupmugen

\brief   A GENIE upgoing muon event generation application.

         *** Synopsis :

           gupmugen    [-h] 
                       [-r run#] 
                        -f flux
                        -n n_of_events
			-d detector side length

         *** Options :

           [] Denotes an optional argument.

           -h Prints out the syntax and exits

           -r Specifies the MC run number 
              [default: 100000000]

           -f Specifies the input flux files
              The general syntax is: `-f simulation:/path/file.data[neutrino_code],...'
              [Notes] 
               - The `simulation' string can be either `FLUKA' or `BGLRS' (so that
                 input data are binned using the correct FLUKA and BGLRS energy and
                 costheta binning). See comments in 
                 - $GENIE/src/Flux/GFlukaAtmo3DFlux.h
                 - $GENIE/src/Flux/GBartolAtmoFlux.h
                 and follow the links to the FLUKA and BGLRS atmo. flux web pages.
               - The neutrino codes are the PDG ones, for numu and anumu only
               - The /path/file.data,neutrino_code part of the option can be 
                 repeated multiple times (separated by commas), once for each 
                 flux neutrino species you want to consider, 
                 eg. '-f FLUKA:~/data/sdave_numu07.dat[14],~/data/sdave_nue07.dat[12]'
                 eg. '-f BGLRS:~/data/flux10_271003_z.kam_nue[12]'
                     
           -n Specifies how many events to generate.

	   -d Specifies detector side length in mm - if not set use default 100m (100000mm)


         *** Examples:

           (1) Generate 100k events (run number 999210)for nu_mu only, using the 
	       sdave_numu07.dat FLUKA flux file (files in /data/flx/), and a detector of
	       side length 50m.

               % gupmugen -r 999210 -n 100000 -d 50000
                       -f FLUKA:/data/flx/sdave_numu07.dat[14]

      
         You can further control the GENIE behaviour by setting its standard 
         environmental variables.
         Please read the GENIE User Manual for more information.

\created April 15, 2011

\author  Jen Truby <jen.truby \at sky.com>
         Oxford University

\cpright Copyright (c) 2003-2011, GENIE Neutrino MC Generator Collaboration
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
#include <iostream.h>

#include <TRotation.h>
#include <TH1D.h>
#include <TH3D.h>
#include <TNtupleD.h>
#include <TLorentzVector.h>

#include "Algorithm/Algorithm.h"
#include "Algorithm/AlgFactory.h"
#include "Base/XSecAlgorithmI.h"
#include "Conventions/Constants.h"
#include "Conventions/Units.h"
#include "EVGDrivers/GFluxI.h"
#include "EVGDrivers/GMCJDriver.h"
#include "EVGDrivers/GMCJMonitor.h"
#include "Numerical/RandomGen.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGLibrary.h"
#include "Utils/StringUtils.h"
#include "Utils/SystemUtils.h"
#include "Utils/UnitUtils.h"
#include "Utils/KineUtils.h"
#include "Utils/CmdLnArgParser.h"

#ifdef __GENIE_FLUX_DRIVERS_ENABLED__
#include "FluxDrivers/GFlukaAtmo3DFlux.h"
#include "FluxDrivers/GBartolAtmoFlux.h"
#endif

using std::string;
using std::vector;
using std::map;
using std::ostringstream;

using namespace genie;
using namespace genie::flux;
using namespace genie::constants;

void            GetCommandLineArgs (int argc, char ** argv);
void            PrintSyntax        (void);
GFluxI *        GetFlux            (void);

TVector3        GetDetectorVertex  (double CosTheta, double Enu);
double          GetCrossSection    (int nu_code, double Enu, double Emu);
TH3D *          Make3DProbHist     (int nu_code);
TH1D *          MakeEmuHistEmpty   (void);
TH1D *          GetEmuHist         (double Enu, double CosTheta, const TH3D* hist3D, TH1D* h1);
double          GetRandomEmu       (const TH1D* hist1D);
double          ProbabilityEmu     (int nu_code, double Enu, double Emu);
void            GenerateUpNu       (GFluxI * flux_driver);

// User-specified options:
//
Long_t          gOptRunNu;                     // run number
string          gOptFluxSim;                   // flux simulation (FLUKA or BGLRS)
map<int,string> gOptFluxFiles;                 // neutrino pdg code -> flux file map
int             gOptNev = -1;                  // exposure - in terms of number of events
double          gOptDetectorSide;              // detector side length, in mm.

// constants
const double NA = 6.0221415e+23; // avogadro's number.
const double a = 2e+6; // a = 2 MeV / (g cm-2)
const double e = 500e+9; // e = 500 GeV

// Defaults:
//
double          kDefOptDetectorSide = 1e+5;     // side length of detector, 100m (in mm)

//________________________________________________________________________________________
int main(int argc, char** argv)
{
  // Parse command line arguments
  GetCommandLineArgs(argc,argv); // Changed by JEN!!

  // get flux driver
  GFluxI * flux_driver = GetFlux();

  // make 3D histograms for nu_code = 14 and -14, containing the probability of a
  // neutrino of energy Ev producing a muon of energy Emu.
  TH3D * h3numu = Make3DProbHist(kPdgNuMu);
  TH3D * h3anumu = Make3DProbHist(kPdgAntiNuMu);

  // Save the 3D histogram
  TFile file0("./hist3D.root","recreate");
  h3numu->Write("hist3D");
  file0.Close();

  // make 1D histograms for nu_code =14, -14 - empty histograms to be filled by
  // a slice of Emu probabilities for a given Enu, CosTheta
  TH1D * h1numu = MakeEmuHistEmpty();
  TH1D * h1anumu = MakeEmuHistEmpty();

  // create ntuple to store values in.
  TNtupleD * UpMuFluxEv = new TNtupleD("UpMuFluxEv","","ievent/D:Emu:Enu:CosTheta:nu_code:flxwght:muwght:vx:vy:vz:xsec");

  // event loop
  for(int iev = 0; iev < gOptNev; iev++) {

    // loop until upgoing nu is generated.
    GenerateUpNu(flux_driver);

    // Get Enu, CosTheta and weight for the generated neutrino
    const TLorentzVector & p4 = flux_driver->Momentum();
    int nu_code = flux_driver->PdgCode();
    double Enu = p4.E();
    double CosTheta = -p4.Pz() / p4.Vect().Mag();
    double FluxWeight = flux_driver->Weight();

    LOG("gupmugen", pNOTICE) 
        << "Generated neutrino: code = " << nu_code
        << ", Ev = " << Enu << "GeV"
        << ", cos(theta) = " << CosTheta 
        << ", weight = " << FluxWeight; 


    // fill the 1D histogram with a slice from the 3D histogram for Enu, Costheta.
    // then  get a random Emu from the 1D histogram, and get the weight for that Emu.
    double Emu = 0;
    double MuWeight = 0;

    if(nu_code ==kPdgNuMu){
      h1numu = GetEmuHist(Enu, CosTheta, h3numu, h1numu);
      MuWeight = h1numu->Integral("width");
      Emu = GetRandomEmu(h1numu);
    }
    else {
      h1anumu = GetEmuHist(Enu, CosTheta, h3anumu, h1anumu);
      MuWeight = h1anumu->Integral("width");
      Emu = GetRandomEmu(h1anumu);
    } 
    

    // Randomly select the point at which the neutrino crosses the detector.
    // assumes a simple cube detector of side length gOptDetectorSide.
    TVector3 Vertex = GetDetectorVertex(CosTheta,Enu);
    double vx = Vertex.X();
    double vy = Vertex.Y();
    double vz = Vertex.Z();

    LOG("gupmugen", pNOTICE) 
      << "Generated up-going muon: Emu = " << Emu << " GeV"
      << "/n at vertex: " << vx << ", " << vy << ", " << vz;
    
    // Save the cross section for the interaction generated.
    double xsec = GetCrossSection(nu_code, Enu, Emu);

    // save all relevant values to the ntuple
    UpMuFluxEv->Fill(iev,Emu,Enu,CosTheta,nu_code,FluxWeight,MuWeight,vx,vy,vz,xsec);

  }

  // Save the ntuple!
  TFile file("./UpMuFluxNtuple.root","recreate");
  UpMuFluxEv->Write("UpMuFlux");
  file.Close();

  // clean-up
  delete flux_driver;

  return 0;
}
//________________________________________________________________________________________

TVector3 GetDetectorVertex(double CosTheta, double Enu)
{

  // Find the point at which the muon crosses the detector. Returns 0,0,0 if the muon misses.
  // Detector is a cube of side length l (=gOptDetectorSide).

  // Get a RandomGen instance
  RandomGen * rnd = RandomGen::Instance();

  // Generate random phi.
  double Phi = 2.*kPi* rnd->RndFlux().Rndm();

  // Set distance at which incoming muon is considered
  double Rad = 0.87*gOptDetectorSide;

  // Set transverse radius of a circle
  double RadTrans = 0.87*gOptDetectorSide;

  // Get necessary trig
  double SinTheta = TMath::Sqrt(1-CosTheta*CosTheta);
  double CosPhi = TMath::Cos(Phi);
  double SinPhi = TMath::Sin(Phi);

  // Compute the muon position at distance Rad.
  double z = Rad * CosTheta;
  double y = Rad * SinTheta * CosPhi;
  double x = Rad * SinTheta * SinPhi;

  // Displace muon randomly on a circular surface of radius RadTrans, 
  // perpendicular to a sphere radius Rad at that position [x,y,z].
  TVector3 vec(x,y,z);                // vector towards selected point
  TVector3 dvec1 = vec.Orthogonal();  // orthogonal vector
  TVector3 dvec2 = dvec1;             // second orthogonal vector
  dvec2.Rotate(-kPi/2.0,vec);         // rotate second vector by 90deg -> Cartesian coords

  // Select a random point on the transverse surface, within radius RadTrans
  double Psi = 2 * kPi * rnd->RndFlux().Rndm();      // rndm angle [0,2pi]
  double random = rnd->RndFlux().Rndm();             // rndm number [0,1]
  dvec1.SetMag(TMath::Sqrt(random)*RadTrans*TMath::Cos(Psi));
  dvec2.SetMag(TMath::Sqrt(random)*RadTrans*TMath::Sin(Psi));
  x += dvec1.X() + dvec2.X();    // x-coord of a point the muon passes through
  y += dvec1.Y() + dvec2.Y();    // y-coord of a point the muon passes through
  z += dvec1.Z() + dvec2.Z();    // z-coord of a point the muon passes through

  // Find out if the muon passes through any side of the detector.

  // Get momentum vector
  double pz = Enu * CosTheta;
  double py = Enu * SinTheta * SinPhi;
  double px = Enu * SinTheta * CosPhi;
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

  while (1){

    flux_driver->GenerateNext();
    cout << "generated next" << endl;

    int nu_code = flux_driver->PdgCode();

    const TLorentzVector & p4 = flux_driver->Momentum();
    //double CosTheta = p4.Pz() / p4.Vect().Mag(); // test with Pz the wrong way! i.e. selecting downgoing
    double CosTheta = -p4.Pz() / p4.Vect().Mag();

    bool keep = (CosTheta<0) && (nu_code==kPdgNuMu || nu_code==kPdgAntiNuMu);
    if(keep) return;
  }
}
//________________________________________________________________________________________
double GetRandomEmu(const TH1D* hist1D)
{

  // Get a random value of Emu from the 1D histogram containing a slice of constant
  // Ev and CosTheta from the 3D histogram of probabilities of interaction.

  double Emu = hist1D->GetRandom();
  return Emu;

}
//________________________________________________________________________________________
TH1D* GetEmuHist(double Enu, double CosTheta, const TH3D* hist3D, TH1D* h1)
{
 
  // Get a 1D slice of constant Ev, CosTheta from the 3D histogram of probabilities.

  // Reset the bin, i.e. remove all contents.
  h1->Reset("C");
 
  //Find appropriate bins for Enu and CosTheta.
  int EnuBin = hist3D->GetYaxis()->FindBin(Enu);
  int CosThetaBin = hist3D->GetZaxis()->FindBin(CosTheta);


  // For those Enu and costheta bins, make a new hist from hist3D.
  for (int i=1; i< h1->GetNbinsX(); i++)
    {
      int EmuBin = i;

      h1->SetBinContent(EmuBin, hist3D->GetBinContent(EmuBin,EnuBin,CosThetaBin));
    }
  
  //h1->Draw();
  return h1;

}
//________________________________________________________________________________________
TH1D* MakeEmuHistEmpty(void)
{

  // Make an empty 1D histogram ready to fill with the probability of getting each value
  // of Emu for a given Enu, CosTheta

  // Bin convention for Emu, consistent with BGLRS
  const double Emumin            = 0.1;
  const int    nEmubinsPerDecade = 10;
  const int    nEmubins      =31;

  // Set up an array of Emu bins.
  Double_t MinLogEmu = log(Emumin);
  Double_t MaxLogEmu = log(10*Emumin);
  Double_t LogBinWidthEmu = ((MaxLogEmu-MinLogEmu)/nEmubinsPerDecade);
  Double_t EmuBinEdges[nEmubins+2];
  for (int i=0; i<=nEmubins+1; i++)
    {
      EmuBinEdges[i] = exp(MinLogEmu + i*LogBinWidthEmu);
    }

  // Create empty 1D histogram.
  TH1D *h1 = new TH1D("hist1D","",nEmubins,EmuBinEdges);

  return h1;
}
//________________________________________________________________________________________
double ProbabilityEmu(int nu_code, double Enu, double Emu)
{

  // Calculate the probability of an incoming neutrino of energy Enu
  // generating a muon of energy Emu.

  double dxsec_dxdy = GetCrossSection(nu_code,Enu,Emu);

  double Int = e * NA * dxsec_dxdy / (a * (1 + Emu/e));

  return Int;

}
//________________________________________________________________________________________
TH3D* Make3DProbHist(int nu_code)
{

// Set up a 3D histogram, with axes Emu, Enu, CosTheta.
// Bin convention is defined at the start.
// Content of each bin is given by the probability of getting a muon
// of energy Emu from a neutrino of energy Enu.

  // Bin convention for Enu, consistent with BGLRS
  const double Enumin            = 0.1;
  const int    nEnubinsPerDecade = 10;
  const int    nEnubins      =31;

  // Bin convention for Emu, consistent with BGLRS
  const double Emumin            = 0.1;
  const int    nEmubinsPerDecade = 10;
  const int    nEmubins      =31;

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
double GetCrossSection(int nu_code, double Enu, double Emu)
{

  // Get the interaction cross section for a neutrino of energy Enu
  // to produce a muon of energy Emu.

  //cout << "Enu: " << Enu << " and Emu: " << Emu << endl;
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

  //cout << "cross section " << dxsec_dxdy << endl;
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
     GFlukaAtmo3DFlux * fluka_flux = new GFlukaAtmo3DFlux;
     atmo_flux_driver = dynamic_cast<GAtmoFlux *>(fluka_flux);
  } else
  if(gOptFluxSim == "BGLRS") {
     GBartolAtmoFlux * bartol_flux = new GBartolAtmoFlux;
     atmo_flux_driver = dynamic_cast<GAtmoFlux *>(bartol_flux);
  } else {
     LOG("gupmugen", pFATAL) << "Uknonwn flux simulation: " << gOptFluxSim;
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
    atmo_flux_driver->SetFluxFile(neutrino_code, filename);
  }
  atmo_flux_driver->LoadFluxData();
  // configure flux generation surface:
  atmo_flux_driver->SetRadii(1, 1);
  // Cast to GFluxI, the generic flux driver interface 
  flux_driver = dynamic_cast<GFluxI *>(atmo_flux_driver);

#else
  LOG("gupmugen", pFATAL) << "You need to enable the GENIE flux drivers first!";
  LOG("gupmugen", pFATAL) << "Use --enable-flux-drivers at the configuration step.";
  gAbortingInErr = true;
  exit(1);
#endif

  return flux_driver;
}
//________________________________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv) // JEN!!
{
// Get the command line arguments

  LOG("gupmugen", pNOTICE) << "Parsing command line arguments";

  CmdLnArgParser parser(argc,argv);

  // help?
  bool help = parser.OptionExists('h');
  if(help) {
    PrintSyntax(); // JEN!!
      exit(0);
  }

  //
  // *** run number
  //
  if( parser.OptionExists('r') ) {
    LOG("gupmugen", pDEBUG) << "Reading MC run number";
    gOptRunNu = parser.ArgAsLong('r');
  } else {
    LOG("gupmugen", pDEBUG) << "Unspecified run number - Using default";
    gOptRunNu = 100000000;
  } //-r



  //
  // *** exposure
  // 

  // in number of events
  bool have_required_statistics = false;
  if( parser.OptionExists('n') ) {
    LOG("gupmugen", pDEBUG) 
        << "Reading number of events to generate";
    gOptNev = parser.ArgAsInt('n');
    have_required_statistics = true;
  }//-n
  if(!have_required_statistics) {
    LOG("gupmugen", pFATAL) 
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
    LOG("gupmugen", pDEBUG)
      << "Reading detector side length";
    gOptDetectorSide = parser.ArgAsDouble('d');
  } else {
    LOG("gupmugen", pDEBUG)
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
    LOG("gupmugen", pDEBUG) << "Getting input flux files";
    string flux = parser.ArgAsString('f');

    // get flux simulation info (FLUKA,BGLRS) so as to instantiate the
    // appropriate flux driver
    string::size_type jsimend = flux.find_first_of(":",0);
    if(jsimend==string::npos) {
       LOG("gupmugen", pFATAL) 
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
        LOG("gupmugen", pFATAL) 
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
           LOG("gupmugen", pFATAL) 
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
       LOG("gupmugen", pFATAL) 
          << "You must specify at least one flux file!"; 
       PrintSyntax();
       gAbortingInErr = true;
       exit(1);
    }

  } else {
    LOG("gupmugen", pFATAL) << "No flux info was specified! Use the -f option.";
    PrintSyntax();
    gAbortingInErr = true;
    exit(1);
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
  if(gOptNev > 0)            { expinfo << gOptNev            << " events";   } 

  LOG("gupmugen", pNOTICE) 
     << "\n"
     << "\n ****** MC Job (" << gOptRunNu << ") Settings ****** "
     << "\n @@ Flux"
     << "\n\t" << fluxinfo.str()
     << "\n @@ Exposure" 
     << "\n\t" << expinfo.str()
     << "\n\n";

}
//________________________________________________________________________________________
void PrintSyntax(void)
{
  LOG("gupmugen", pFATAL) 
   << "\n **Syntax**"
   << "\n gupmugen [-h]"
   << "\n           [-r run#]"
   << "\n            -f simulation:flux_file[neutrino_code],..."
   << "\n            -n n_of_events,"
   << "\n           [-d detector side length (mm)]"
   << "\n"
   << " Please also read the detailed documentation at http://www.genie-mc.org"
   << "\n";
}
//________________________________________________________________________________________

