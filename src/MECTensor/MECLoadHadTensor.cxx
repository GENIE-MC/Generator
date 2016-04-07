//_______________________________________________________________
/*

The key bits of code that pulls the prediction from the hadron tensor

Originated for us in an insight by Federico Sanchez, after working
with the Valencial model in FORTRAN in 2013, that the optimal way to 
push this to an accept/reject loop is to separate the model
at the point where the leptonic and hadronic tensors are contracted.
The leptonic tensor codes nu and nubar, energy, and Q-value,
while the hadronic tensor has no dependence on those.

Manolo Vicente Vacas and Juan Nieves modified their FORTRAN to
make this split work.  Sanchez made the first skeleton 
accept-reject loop in C++ that called their FORTRAN directly
used for the Gran, Nieves, Sanchez, Vicente Vacas 2013 paper

This code (the XSec method) is a C++ reimplementation of the
100 line contraction of the hadron and lepton tensor fortran code
and produces output identical within machine precision.

The physics coded in the HadTensor inputs and here is described:

J. Nieves, I. Ruiz Simo, M.J. Vicente Vacas
"Inclusive quasi-elastic neutrino reactions" PRC 83 (2011) 045501

R. Gran, J. Nieves, F. Sanchez, M.J. Vicente Vacas
"Neutrino-nucleus quasi-elastic and 2p2h interactions up to 10 GeV"
PRD 88 (2013) 113007

J. Schwehr, D. Cherdack, R. Gran
"GENIE implementation of IFIC Valencia model for
 QE-like 2p2h neutrino-nucleus cross sections"
arXiv:1601.02038 whitepaper (maybe later to JINST).

2014/09/15 first working version by ~J.Schwehr
2015/01/14 all neutrino flavors, oxygen and carbon ~J.Schwehr
2015/11/xx streamline the codes and methods J.Schwehr, G. Perdue
2015/02/06 added spline generation code ~R. Gran
2016/02/20 general clean-up and fix a problem with PDD ~R. Gran
2016/03/01 structural modifications to cover all nuclei ~R. Gran

*/
//_______________________________________________________________


// --- Includes --- //

// 
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <string>
#include <fstream>
#include <cassert>
#include <string>
#include <iostream>
#include <vector>
#include <algorithm>

#include "Conventions/Constants.h"
#include "HadronTransport/INukeHadroData.h"
#include "Messenger/Messenger.h"
#include "Numerical/Spline.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "PDG/PDGLibrary.h"
#include "Utils/CmdLnArgParser.h"
#include "Numerical/BLI2D.h"
#include "MECTensor/MECLoadHadTensor.h"

#include <TSystem.h>
#include <TFile.h>
#include <TGraph2D.h>
#include <TH2.h>

using std::ostringstream;
using std::istream;
using std::ios;
using std::endl;

using namespace genie;
using namespace genie::constants;

//________________________________________________________________
// 
MECLoadHadTensor * MECLoadHadTensor::fInstance = 0;
//________________________________________________________________
//
MECLoadHadTensor::MECLoadHadTensor()
{

  // This list should be abstracted to a configuration file
  // Along with the known target method.
  
  //  three targets have explicit tensor tables
  fKnownTensors.push_back(kPdgTgtC12); // C12 
  fKnownTensors.push_back(kPdgTgtO16); // O16
  fKnownTensors.push_back(1000140280); // Si28
  fKnownTensors.push_back(kPdgTgtCa40); // Ca40, also for Ar40
  fKnownTensors.push_back(1000280560); // Ni56, actually it is pseudo-Fe56
  fKnownTensors.push_back(1000561120); // Ba112, actually it is pseudo-Cd112
  fKnownTensors.push_back(1001042080); // Rf208, actually it is pseudo-Pb208

  // why pseudo ?  The had tensor calculation requires the nuclear density function
  // simple two-parameter hartree fock Fermi function for heavy nuclei .
  // and in principle modified harmonic oscillator function for light nuclei
  // we are never really going to want isoscalar Ni56, so I used density for Fe56.
  // likewise never want Rf208, I used the density for Pb208
  // likewise never want Ba112, I used the density for Cd112

  // this one loads all known targets when instantiated, regardless of what the user wants.
  // is there a point to being efficient, check if the requested one is loaded ?
  // TODO  there is a fail if the spline file is missing one of the known targets
  // and the user requests a mix of two targets.
  for(std::vector<int>::const_iterator it=fKnownTensors.begin(); it!=fKnownTensors.end(); ++it){

    this->LoadTensorTables(*it);

    // this method will fail silently if the target list here does not match later.
    
  }
  
  fInstance = 0;
}
//________________________________________________________________
// 
MECLoadHadTensor::~MECLoadHadTensor()
{
}
//________________________________________________________________
//
MECLoadHadTensor * MECLoadHadTensor::Instance()
{
  // check if we have a full instance of this class.  If not, make one.
  if(fInstance == 0){
    //LOG("MECLoadHadTensor", pINFO) << "MECLoadHadTensor late initialization";
    static MECLoadHadTensor::Cleaner cleaner;
    cleaner.DummyMethodAndSilentCompiler();
    fInstance = new MECLoadHadTensor();
  }
  
  // now we have one for sure, return one to whomever requested it.
  return fInstance;
}

//_________________________________________________________________
bool MECLoadHadTensor::KnownTarget(int targetpdg)
{

  if(std::count(fKnownTensors.begin(), fKnownTensors.end(), targetpdg)!=0){
    return true;
  }

  // abstract these limits.  They are used also in XSecFour.
  // hard coding things in two locations isn't a good idea.
  // also, future-person might give a had tensor and want to not allow scaling.

  int Arequest = pdg::IonPdgCodeToA(targetpdg);
  //int Zrequest = pdg::IonPdgCodeToZ(targetpdg);
  if(Arequest >= 10){
    return true;
  }

  return false;
  
} 

bool MECLoadHadTensor::KnownTensor(int targetpdg)
{
  return std::count(fKnownTensors.begin(), fKnownTensors.end(), targetpdg)!=0;
  
}

//_________________________________________________________________
// Create the hadron tensor tables
void MECLoadHadTensor::LoadTensorTables(int targetpdg)
{
  
  // define dimensions of data in the hadron tensor files
  int nwpoints = 5;
  const int nq0points = 120;  //240;
  const int nqzpoints = 120;  //240;
  const int nq0qzpoints = nq0points*nqzpoints;
  double arraystep = 0.01;  // GeV
  // if later we use tables that are not 120x120
  // then extract these constants to the config file with the input table specs
  // or find a way for the input tables to be self-descriptive.


  // will start with ${GENIE}/data/evgen/mectensor/nieves/HadTensor...dat
  string modelLocation = string("/data/evgen/mectensor/nieves");
  string modelFileStart = "HadTensor120-Nieves-";
  string modelFileEnd = "-20150210.dat";


  
  if(!KnownTarget(targetpdg)){
    LOG("MEC", pERROR) << "LoadTensorTables asked for table for " << targetpdg << ", for which we don't have a table";
    return;
  }

  // define directory of hadron tensor files
  // Ideally, the xml configuration can override the default location
  string data_dir = string(gSystem->Getenv("GENIE")) + modelLocation;

  // define arrays to fill from data files
  double hadtensor_q0_array[nq0points];
  double hadtensor_qz_array[nqzpoints];
  double hadtensor_w_array[nwpoints][nq0qzpoints];

  // fill q0 array (in GeV)
  // 240 5MeV bins or 120 10 MeV bins
  for (int a = 0; a < nq0points; a++){
    hadtensor_q0_array[a]=double(a+1)*arraystep;
  }

  // fill qz array (in GeV)
  for (int a = 0; a < nqzpoints; a++){
    hadtensor_qz_array[a]=double(a+1)*arraystep;
  }

  
  // possible future feature, allow a model to not deliver Delta tensors.
  std::map<HadTensorTables::EHadTensorType, std::string> tensorTypeNames;
  tensorTypeNames[HadTensorTables::kFullAll]  = "FullAll";
  tensorTypeNames[HadTensorTables::kFullpn]   = "Fullpn";
  tensorTypeNames[HadTensorTables::kDeltaAll] = "DeltaAll";
  tensorTypeNames[HadTensorTables::kDeltapn]  = "Deltapn";

  // iterate over all four hadron tensor types in the map above.
  for(int tensorType=0; tensorType < HadTensorTables::kNHadTensorTypes; ++tensorType){
    // build filenames from the bits of string and twine and gubbins
    ostringstream datafile;

    datafile << data_dir << "/" << modelFileStart << targetpdg << "-" 
	     << tensorTypeNames[(HadTensorTables::EHadTensorType)tensorType]
	     << modelFileEnd;

    // make sure data files are available
    LOG("MEC", pDEBUG) << "Asserting that file " << datafile.str().c_str() << " exists...";
    //std::cout <<  "Asserting that file " << datafile.str().c_str() << " exists..." << std::endl;
      
    assert (! gSystem->AccessPathName(datafile.str().c_str()));
  
    // read data file
    ReadHadTensorqzq0File(datafile.str(), nwpoints, nqzpoints, nq0points, hadtensor_w_array);
  
    //loop over all 5 tensors 
    for (int i=0; i<nwpoints; i++){
   
      // "save" xsecs into non uniform grid
      genie::BLI2DNonUnifGrid *hadTensorGrid = new genie::BLI2DNonUnifGrid(nq0points, nqzpoints, hadtensor_q0_array, hadtensor_qz_array, hadtensor_w_array[i]);

      // we are storing these in a map
      // map key is pdg code of target
      //   std::map<int, HadTensorTables> fTargetTensorTables;
      // seems nested in a way that defies scrutiny, but I guess its okay.
      fTargetTensorTables[targetpdg].fTables[(HadTensorTables::EHadTensorType)tensorType].push_back(hadTensorGrid);

      
    }
  }// end loop over tables (probably need to clear things here...)  
}

//___________________________________________________________________
void MECLoadHadTensor::ReadHadTensorqzq0File( 
					     string filename, int nwpoints, int nqzpoints, int nq0points, double hadtensor_w_array[][14400]){
  // ------ open and check file ---------- //

  // open file
  std::ifstream tensor_stream(filename.c_str(), ios::in);

  // check file exists
  if(!tensor_stream.good()){
    LOG("MEC", pERROR) << "Bad file name: " << filename;
    return;
  }
  
  //--------- read file ---------- //

  double temp;
  
  for (int ij = 0; ij < (nqzpoints*nq0points); ij++){  
    for (int k = 0; k < nwpoints; k++){
      tensor_stream >> temp;
      hadtensor_w_array[k][ij]=temp;
    }}

}

// There is a note here because I'm not sure where to put it.
// The had tensor calculation requires
// parameters that describe the nuclear density.
// Juan takes them as shown in
// C. Garcia-Recio, Nieves, Oset Nucl.Phys.A 547 (1992) 473-487 Table I
// and simplifies that the p and n density are not necessarily the same?
// Standard tables such as deVries et al are similar ~5%
// This is the only nuclear input on the hadron side of the Hadron Tensor.
// and is what makes PseudoPb and PseudoFe different than Ni56 and Rf208


double MECLoadHadTensor::Qvalue(int targetpdg, int nupdg)
{

  int nearpdg = targetpdg;

  if(nupdg > 0){
    // get Z+1 for CC interaction e.g. 1000210400 from 1000200400
    nearpdg = targetpdg + 10000;
  } else {
    // get Z-1 for CC interaction e.g. 1000210400 from 1000200400
    nearpdg = targetpdg - 10000;
  }
  
  TParticlePDG *partf = PDGLibrary::Instance()->Find(nearpdg);
  TParticlePDG *parti = PDGLibrary::Instance()->Find(targetpdg);

  if(NULL == partf || NULL == parti){
    // maybe not every valid nucleus in the table has Z+1 or Z-1
    // for example, Ca40 did not.
    LOG("MEC", pERROR) << "Can't get qvalue, nucleus " << targetpdg << " or " << nearpdg << " is not in the table of nuclei in /data/evgen/pdg ";    
  }

  double massf = partf->Mass();
  double massi = parti->Mass();

  // keep this here as a reminder of what not to do
  //+/- emass;  // subtract electron from Z+1, add it to Z-1
  // the lookup table gives the nucleus mass, not the atomic mass, so don't need this.

  //LOG("MECLoadHadTensor", pINFO) << "Qvalue " << massf - massi << " final " << massf << " initial " << massi << " " << targetpdg;
  // to see that the Qvalue matters, override the value for Ar40 with Ca40
  // which is a relatively big shift.  By itself that boosts the cross section.

  
  return massf - massi;  // in GeV , thats what I want.



  
  // From Table I in Nieves et al. PRC 70 055503 (2004)
  // Also in Table I in Nieves et al. PRC 73 025504 (2006) hep-ph/0511204
  //
  // With the replacement table for GENIE's old nucleus mass
  // G. Audi et al., CPC(HEP & NP) v36 n12 (Dec 2012) 1157-1286
  // don't need the hardcoded version.  Keep them here commented out for posterity.
  
  //if(targetpdg == kPdgTgtC12){
  //  if (nupdg > 0) myQvalue = 16.827/1000.; // GeV neutrino
  //  else myQvalue = 13.880/1000.; // GeV anti-neutrino
  //} else if(targetpdg == kPdgTgtO16){
  //  if (nupdg > 0) myQvalue = 14.906/1000.; // GeV neutrino
  //  else myQvalue = 10.931/1000.; // GeV anti-neutrino
  //} else if(targetpdg == kPdgTgtCa40){
  //  if (nupdg > 0) myQvalue = 13.809/1000.; // GeV neutrino
  //  else myQvalue = 1.822/1000.; // GeV anti-neutrino
  //}
  
}

//____________________________________________________________________
// Xsec

double MECLoadHadTensor::XSec(int targetpdg, int tensorpdg, int nupdg, double Enu, double Tmu, double Costheta,  vector <genie::BLI2DNonUnifGrid *> HadTensor)
{
  // These units should deliver 10^{41} cm^2 / GeV for d2sigma/(dTmu dcos_mu)

  // This implements the lepton tensor contraction from FORTRAN code lectorV42.f
  // with the hadron tensor provided in tabular form
  // the lepton tensor is expressed formally in Nieves PRC 70 (2004) 055503
    
  //--- variables and constants ---//

  // lepton:
  double Mlep = 0.0;
  if( nupdg == kPdgNuE || nupdg == kPdgAntiNuE ){
    Mlep = kElectronMass;
  }
  else if( nupdg == kPdgNuMu || nupdg == kPdgAntiNuMu ){
    Mlep = kMuonMass;
  }
  else if( nupdg == kPdgNuTau || nupdg == kPdgAntiNuTau){
    Mlep = kTauMass;
  }
  double xsec;
  TLorentzVector v4lep;
  TLorentzVector v4Nu = TLorentzVector(0,0,Enu,Enu);   // assuming traveling along z:
  TLorentzVector v4q;
  double q0nucleus;
  double facconv = 0.0389391289e15; // std::pow(0.19733,2)*1e15;
  // These units should deliver 10^{41} cm^2 / GeV for d2sigma/(dTmu dcos_mu)

  double myQvalue = this->Qvalue(targetpdg, nupdg);

  //--- more kinematics to precompute ---//
  //angles
  double Sintheta = 1. - Costheta * Costheta;
  if(Sintheta < 0.0)Sintheta = 0.0;
  else Sintheta = TMath::Sqrt(Sintheta);   
  
  double Cosh = TMath::Cos(TMath::ACos(Costheta)/2.);
  double Sinh = TMath::Sin(TMath::ACos(Costheta)/2.);
  // lepton
  v4lep.SetE( Tmu + Mlep );
  // energy transfer from the lepton
  double q0 = v4Nu.E() - v4lep.E();
  // energy transfer that actually gets to the nucleons
  q0nucleus = q0 - myQvalue;

  // Define some calculation placeholders
  double part1, part2;
  double modkprime ;
  double w1, w2, w3, w4, w5;  
  double wtotd[5];

  //--- event selection loop ---//
  if (q0nucleus <= 0){
    //nothing transfered to nucleus thus no 2p2h
    xsec = 0.;
  } else {
    // energy was transfered to the nucleus
    modkprime = std::sqrt(TMath::Power(v4lep.E(),2)-TMath::Power(Mlep,2));
    v4lep.SetX(modkprime*Sintheta); 
    v4lep.SetY(0);
    v4lep.SetZ(modkprime*Costheta);

    //q: v4q = v4Nu - v4lep;
    v4q.SetE(q0nucleus);
    v4q.SetX(v4Nu.X() - v4lep.X());
    v4q.SetY(v4Nu.Y() - v4lep.Y());
    v4q.SetZ(v4Nu.Z() - v4lep.Z());
    //remember: qm = v4q.Vect().Mag()
  
    // keep this for step-by-step debugging of the calculation.
    // this overrides with specific kinematics I can compare to FORTRAN routine
    //std::cout << "RIK special " << q0nucleus << std::endl;
    //for(int i=0 ; i < 5; i++){
    //  wtotd[i]=HadTensor[i]->Evaluate(0.0190,0.000173);
    //  std::cout << "RIK special " << i << wtotd[i] << std::endl;
    //}

    for (int i=0 ; i < 5; i++){
      wtotd[i]=HadTensor[i]->Evaluate(v4q.Vect().Mag(),v4q.E());
    }

    // calculate hadron tensor components
    // these are footnote 2 of Nieves PRC 70 055503
    double W00 = wtotd[0]; double W0Z = wtotd[1]; double WXX = wtotd[2];
    double WXY = wtotd[3]; double WZZ = wtotd[4];

    w1=WXX/2.;
    w2=(W00+WXX+(q0*q0/(v4q.Vect().Mag()*v4q.Vect().Mag())
			   *(WZZ-WXX))
	  -(2.*q0/v4q.Vect().Mag()*W0Z))/2.;
    w3=WXY/v4q.Vect().Mag();
    w4=(WZZ-WXX)/(2.*v4q.Vect().Mag()*v4q.Vect().Mag());
    w5=(W0Z-(q0/v4q.Vect().Mag()*(WZZ-WXX)))/v4q.Vect().Mag();
    //w6 we have no need for w6, noted at the end of IIA.
    
    // adjust for anti neutrinos
    
    if (nupdg < 0) w3 = -1. * w3;

    // calculate cross section, in parts
    double xw1 = w1*Costheta;
    double xw2 = w2/2.*Costheta;
    double xw3 = w3/2.*(v4lep.E()+modkprime-(v4lep.E()+v4Nu.E())*Costheta);
    double xw4 = w4/2.*(Mlep*Mlep*Costheta+2.*v4lep.E()*(v4lep.E()+modkprime)*Sinh*Sinh);
    double xw5 = w5*(v4lep.E()+modkprime)/2.;
    part1 = xw1 - xw2 + xw3 + xw4 - xw5;

    double yw1 = 2.*w1*Sinh*Sinh;
    double yw2 = w2*Cosh*Cosh;
    double yw3 = w3*(v4lep.E()+v4Nu.E())*Sinh*Sinh;
    double yw4 = Mlep*Mlep*part1/(v4lep.E()*(v4lep.E()+modkprime));
    part2 = yw1 + yw2 - yw3 + yw4;

    xsec = modkprime*v4lep.E()*kGF2*2./kPi*part2*facconv;

    if( ! (xsec >= 0.0 ) ){
      //  !(xsec >= 0.0) traps negative numbers and nan's.
      //  Sometimes Costheta is just over 1.0 due to fluke or numerical precision, so Sintheta would be undefined.
      //  I guess, take the std::cout version below and turn it into an info or warning.
      std::cout << "RIK is bogus xsec b ? " << xsec << " " << part1 << " " << part2
		<< " w[0] " << w1 << " " << w2 << " " << w3 << " " << w4 << " " << w5
		<< " wtotd[] " << wtotd[0] << " " << wtotd[1] << " " << wtotd[2] << " " << wtotd[3] << " " << wtotd[4] << " "
		<< " vec " << v4q.Vect().Mag() << " " << q0  << " " << v4q.Px() << " " << v4q.Py() << " " << v4q.Pz()
		<< " v4qX " << v4Nu.X() << " " << v4lep.X() << " " << Costheta << " " << Sintheta << " " << modkprime
		<< std::endl;
    } 
    /* else {
      // But keep like this with these print lines for now.  They are
      // for Rik to check bit by bit with the Nieves-Vicente FORTRAN calculation.
      if(fabs(v4q.P() - 0.042) < 0.0001){
	std::cout << "RIK calc check1 " << nupdg << " " << v4Nu.E() << " " << v4lep.E() << " " << Costheta << " " << q0 << " " << xsec << std::endl;
	std::cout << "RIK calc check2 " << part1 << " " << part2 << " " << w1 << " " << w2 << " " << w3 << " " << w4 << " " << w5 << std::endl;
	std::cout << "RIK calc check3 " << wtotd[0] << " " << wtotd[1] << " " << wtotd[2] << " " << wtotd[3] << " " << wtotd[4] << std::endl;
	std::cout << "RIK calc check4 " << modkprime << " " << Cosh << " " << Sinh << " " << v4lep.E() << std::endl;
	std::cout << "RIK calc check5 " << part1 << " " << xw1 << " " << xw2 << " " << xw3 << " " << xw4 << " " << xw5 << std::endl;
	std::cout << "RIK calc check6 " << part2 << " " << yw1 << " " << yw2 << " " << yw3 << " " << yw4 << std::endl;
	} 
	
     } */

  } 

  //LOG("MECLoadHadTensor", pINFO) << "xsec=  " << xsec * (1.0E-41 * units::cm2) << " (1.0E-41 cm2)";
  
  return xsec * (1.0E-41 * units::cm2);
}

//___________________________________________________________________________
double MECLoadHadTensor::TotalXsecAtE(int targetpdg, int nupdg, double Enu){

  // JS: this needs to move to the GSL Integrator

  // this does a numerical integral in q0 q3 space using the known
  // 1.2 GeV limits of the model, so is the natural, optimal place to integrate.
  // If a future models are different, then a more asbtract version is needed
  // and/or push this to the GSL Integrator

  // at 3.0 GeV neutrino Carbon, this delivers 1467.16e-41 cm2 
  // This is good to 0.5% with the 1200, 0.001 step size integration.

  // if I'm calling this with the wrong nucleus, return 0.0 instead of crash.
  // not sure who would be calling this with the wrong nucleus.
  if(!KnownTarget(targetpdg))return 0.0;
   
  double maxstep = 1200;
  double stepsize = 0.001;
  long double totalXsec = 0.0;  // need more precision for the integration.
  double maxXsec = -0.0;

  int lpdg = 0;
  if(TMath::Abs(nupdg)==16)lpdg=15;
  else if(TMath::Abs(nupdg)==14)lpdg = 13;
  else if(TMath::Abs(nupdg)==12)lpdg = 11;
  else {
    // oh for dumb, this is not a neutrino?
    LOG("MEC", pERROR) << "Trying to get XSec for not a neutrino " << lpdg;
    assert(false);
  }
  
  double lmass = PDGLibrary::Instance()->Find( lpdg )->Mass();

  for(int iq3=0; iq3<=maxstep; iq3++){
    for(int iq0=0; iq0<=iq3; iq0++){
      double dq3 = stepsize * (double)iq3;
      double dq0 = stepsize * (double)iq0;

      double tmu = 0.0;
      double cost = 0.0;
      double area = 0.0;

      GetTmuCostFromq0q3(dq0,dq3,Enu,lmass,tmu,cost,area);

      if(tmu < 0.0)continue;
      if(area <= 0.0)continue;
      if(cost < -1.0)continue;

      double XsecArray[4];
      this->XSecFour(targetpdg, nupdg, Enu, tmu, cost, XsecArray, false);
      // called this without the Delta.

      double myXsec = XsecArray[0];
      
      totalXsec += (long double) ( myXsec * area  * stepsize * stepsize);
      if(myXsec > maxXsec)maxXsec = myXsec;
      
    }
  }

  LOG("MECLoadHadTensor", pINFO) << "maxXsec " << maxXsec << " " << Enu ;

  
  return totalXsec; 
  
}

//_______________________________________________________________________
double MECLoadHadTensor::GetTmuCostFromq0q3(double dq0, double dq3, double Enu, double lmass, double &tmu, double &cost, double &area){

  tmu = Enu - dq0 - lmass;
  if(tmu < 0.0){
    cost = -999;
    tmu = -999;
    area = 0.0;
    return -999;
  }

  double thisE = tmu + lmass;
  double thisp = sqrt( thisE * thisE - lmass * lmass);
  double numerator =  Enu * Enu + thisp * thisp - dq3 * dq3; 
  double denominator = 2.0 * thisp * Enu;
  if(denominator <= 0.0 ){
    cost = 0.0;
    if(denominator < 0.0){
      return -999;
    }
  }
  else cost = numerator / denominator;

  if(TMath::Abs(numerator) > TMath::Abs(denominator)){
    cost = -999;
    tmu = -999;
    area = 0.0;
    return -999;
  }
  
  // xCrossSect is not yet in the right units for this particular case.
  // need areaElement to go dsigma/dTmudcost to dsigma/dq0dq3
  // Recompute the area element jacobian
  // is this coded somewhere else in GENIE already ?
  
  // dT/dq0 x dC/dq3 - dT/dq3 x dC/dq0
  double areaElement = 0.0;
  //double veryCloseToZero = 0.000001;  // in GeV, this should be safe.
  double insqrt = 0.0;
  numerator = 0.0;
  insqrt = Enu * Enu - 2.0 * Enu * dq0 + dq0 * dq0 - lmass * lmass;
  numerator = dq3 / Enu;
  if(insqrt < 0.0)areaElement=0.0;
  else areaElement = numerator / TMath::Sqrt(insqrt);
  area = areaElement;
  
  return 0;
  
}



void MECLoadHadTensor::XSecFour(int targetpdg, int nupdg, double Enu, double Tmu, double CosTheta, double* FourXSec, bool delta){

  // All calls for a cross section must come through this method
  // To generate XSec for nuclei other than those with hadron tensors
  // We need to pull both the Full XS and the pn initial state fraction
  // Non-isoscalar nuclei are beyond the original published Valencia model
  // and scale with A according to the number of pp, pn, or nn pairs
  // the probe is expected to find.

  // There is some by-hand optimization here, skipping the delta part when
  // the spline calls are requesting the total cross section only.
  // Possible future models without a Delta had tensor would also use that
  // flag to call this without computing the Delta part.
    
  // passing a pointer double* FourXSec seems unsafe.
  // Have a real programmer consider this.
  
  int Arequest = pdg::IonPdgCodeToA(targetpdg);
  int Zrequest = pdg::IonPdgCodeToZ(targetpdg);

  int tensorpdg = targetpdg;

  if( !KnownTensor(tensorpdg) ){

    // Code which had tensor to use for different ranges of A.
    // If someone sends in a new model, this might need to be
    // abstracted to a model-specific configuration file.
    // Also needs to match the KnownTarget method.
    
    if (Arequest < 10){
      LOG("MEC", pWARN) << "Asked to scale to a nucleus " << targetpdg << " which we don't know yet.";
      FourXSec[0] = 0.0;  FourXSec[1] = 0.0; FourXSec[2] = 0.0;  FourXSec[3] = 0.0;
      return;
      // refuse to do He, Li, Be
    }else if( Arequest >= 10 && Arequest < 15){
      tensorpdg = kPdgTgtC12;
      //} else if( A >= 14 && A < 15){ AND CHANGE <=14 to <14.
      // tensorpdg = kPdgTgtN14;      
    } else if( Arequest >= 15 && Arequest < 22){
      tensorpdg = kPdgTgtO16;      
    } else if( Arequest >= 22 && Arequest < 33){
      // of special interest, this gets Al27 and Si28
      tensorpdg = 1000140280;
    } else if(Arequest >= 33 && Arequest < 50){
      // of special interest, this gets Ar40 and Ti48
      tensorpdg = kPdgTgtCa40;      
    } else if( Arequest >= 50 && Arequest < 90){
      // pseudoFe56, also covers many other ferrometals and Ge
      tensorpdg = 1000280560;
    } else if( Arequest >= 90 && Arequest < 160){
      // use Ba112 = PseudoCd.  Row5 of Periodic table useless. Ag, Xe?
      tensorpdg = 1000561120;
    } else if( Arequest >= 160 ){
      // use Rf208 = pseudoPb
      tensorpdg = 1001042080;
    } else {
      LOG("MEC", pWARN) << "Asked to scale to a nucleus " << targetpdg << " which we don't know yet.";
      FourXSec[0] = 0.0;  FourXSec[1] = 0.0; FourXSec[2] = 0.0;  FourXSec[3] = 0.0;
      return;
    }
    
  }

  //----------
  // This block runs regardless of whether we scale
  
  FourXSec[0] = this->XSecFullAll(targetpdg, tensorpdg, nupdg, Enu, Tmu, CosTheta);
  FourXSec[1] = this->XSecFullpn(targetpdg, tensorpdg, nupdg, Enu, Tmu, CosTheta);
    // possible time saver, the spline only needs the Full categories.  Suppress Delta?
  if(delta)FourXSec[2] = this->XSecDeltaAll(targetpdg, tensorpdg, nupdg, Enu, Tmu, CosTheta);
  if(delta)FourXSec[3] = this->XSecDeltapn(targetpdg, tensorpdg, nupdg, Enu, Tmu, CosTheta);

  //LOG("MECLoadHadTensor", pINFO) << "FourXsec " << FourXSec[0] << " " << FourXSec[1] << " " << FourXSec[2] << " " << FourXSec[3] << " " << tensorpdg << " " << targetpdg; 


  
  //------
  // now again, run this only if targetpdg != tensorpdg.
  
  if( !KnownTensor(targetpdg) ){
    // if we need to scale, figure it out here.
    
    double PP = Zrequest;
    double NN = Arequest - PP;
    double P = pdg::IonPdgCodeToZ(tensorpdg);
    double N = pdg::IonPdgCodeToA(tensorpdg) - P;
    
    double scalepn = TMath::Sqrt( (PP*NN)/(P*N) );
    double scalepp = TMath::Sqrt( (PP * (PP - 1.)) / (P * (P - 1.)) );
    double scalenn = TMath::Sqrt( (NN * (NN - 1.)) / (N * (N - 1.)) );

    //LOG("MECLoadHadTensor", pINFO) << "scale pn pp nn for " << targetpdg << " " << tensorpdg << " : " << scalepn << " " << scalepp << " " << scalenn;
    
    // this is an approximation in at least three senses.
    // we are scaling from an isoscalar nucleus using p and n counting
    // we are not using the right qvalue in the had tensor
    // we are not scaling the Delta faster than the non-Delta.
    // the guess is that these are good approximations.
    // a test we could document is to scale from O16 to N14 or C12 using this algorithm
    // and see how many percent deviation we see from the full calculation.      
    
    
    double tempAll = FourXSec[0];
    double temppn = FourXSec[1]*scalepn;
    double tempDeltaAll = 0.0;
    double tempDeltapn = 0.0;
    if(delta){
      tempDeltaAll = FourXSec[2];
      tempDeltapn = FourXSec[3]*scalepn;
    }
    
    if(nupdg > 0){
      tempAll = FourXSec[1] * scalepn + (FourXSec[0] - FourXSec[1]) * scalenn;
      if(delta)tempDeltaAll = FourXSec[3] * scalepn + (FourXSec[2] - FourXSec[3]) * scalenn;
      // here, could add or shave a little Delta.
      // maybe test that the Delta is not more than the total?
    } else {
      tempAll = FourXSec[1] * scalepn + (FourXSec[0] - FourXSec[1]) * scalepp;
      if(delta)tempDeltaAll = FourXSec[3] * scalepn + (FourXSec[2] - FourXSec[3]) * scalepp;
    }
    
    // now overwrite the cross section with the scaled version
    FourXSec[0] = tempAll;
    FourXSec[1] = temppn;
    FourXSec[2] = tempDeltaAll;
    FourXSec[3] = tempDeltapn;
  }
  
}

//__________________________________________________________
// specific XSec from the four specific Hadron Tensors

double MECLoadHadTensor::XSecFullAll(int targetpdg, int tensorpdg, int nupdg, double Enu, double Tmu, double CosTheta) {
  if(!KnownTensor(tensorpdg)) return 0.0;
  if(!KnownTarget(targetpdg)) return 0.0;
  // return 0.0 or better return obviously error -1 ?

  return XSec( targetpdg, tensorpdg, nupdg, Enu, Tmu, CosTheta, fTargetTensorTables[tensorpdg].fTables[HadTensorTables::kFullAll] );

}

double MECLoadHadTensor::XSecFullpn(int targetpdg, int tensorpdg, int nupdg, double Enu, double Tmu, double CosTheta) {

  if(!KnownTensor(tensorpdg)) return 0.0;
  if(!KnownTarget(targetpdg)) return 0.0;
    // return 0.0 or better return obviously error -1 ?
  
  return XSec( targetpdg, tensorpdg, nupdg, Enu, Tmu, CosTheta, fTargetTensorTables[tensorpdg].fTables[HadTensorTables::kFullpn] );

}

double MECLoadHadTensor::XSecDeltaAll(int targetpdg, int tensorpdg, int nupdg, double Enu, double Tmu, double CosTheta) {

  if(!KnownTensor(tensorpdg)) return -0.0;
  if(!KnownTarget(targetpdg)) return 0.0;
  // return 0.0 or better return obviously error -1 ?
  
  return XSec( targetpdg, tensorpdg, nupdg, Enu, Tmu, CosTheta, fTargetTensorTables[tensorpdg].fTables[HadTensorTables::kDeltaAll] );

}

double MECLoadHadTensor::XSecDeltapn(int targetpdg, int tensorpdg, int nupdg, double Enu, double Tmu, double CosTheta) {

  if(!KnownTensor(tensorpdg)) return 0.0;
  if(!KnownTarget(targetpdg)) return 0.0;
  // return 0.0 or better return obviously error -1 ?
  
  return XSec( targetpdg, tensorpdg, nupdg, Enu, Tmu, CosTheta, fTargetTensorTables[tensorpdg].fTables[HadTensorTables::kDeltapn] );


}

