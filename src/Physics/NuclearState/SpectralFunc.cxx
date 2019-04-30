//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

         Steven Gardiner <gardiner \at fnal.gov>
         Fermi National Accelerator Laboratory

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ May 01, 2012 - CA
   Pick spectral function data from $GENIE/data/evgen/nucl/spectral_functions

 @ April 29, 2019 - SG
   Full rewrite to use finely-binned histograms (tightly-peaked distributions led
   to very slow rejection sampling), XML-based configuration of spectral function
   data files with lazy initialization
*/
//____________________________________________________________________________

#include <sstream>
#include <string>

#include "TF2.h"
#include "TGraph2D.h"
#include "TH2D.h"
#include "TNtupleD.h"
#include "TRandom.h"
#include "TSystem.h"

#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Messenger/Messenger.h"
#include "Physics/NuclearState/SpectralFunc.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/Numerical/RandomGen.h"

// Create a functor class in this anonymous namespace for now. We will use it
// to facilitate numerical integration of TGraph2D objects in genie::SpectralFunc.
// TODO: when GENIE moves to C++11, replace this nonsense with a lambda expression.
namespace {

  // Store a pointer to the TGraph2D to be integrated, but don't take ownership
  class FunctorThatShouldBeALambda {
    public:
      FunctorThatShouldBeALambda(TGraph2D* graph) : fGraph( graph ) {}

      double operator()(double* x, double*) {
        return fGraph->Interpolate( x[0], x[1] );
      }

    private:
      TGraph2D* fGraph;
  };

  // Replace this with std::to_string when we switch to C++11
  std::string replace_with_std_to_string(int an_integer) {
    std::ostringstream oss;
    oss << an_integer;
    return oss.str();
  }

}

using namespace genie;
using namespace genie::constants;
using namespace genie::controls;

//____________________________________________________________________________
SpectralFunc::SpectralFunc() : NuclearModelI("genie::SpectralFunc")
{

}
//____________________________________________________________________________
SpectralFunc::SpectralFunc(string config) :
  NuclearModelI("genie::SpectralFunc", config)
{

}
//____________________________________________________________________________
SpectralFunc::~SpectralFunc()
{
  // Delete the TH2D objects from the spectral function map
  std::map<int, TH2D*>::iterator begin = fSpectralFunctionMap.begin();
  std::map<int, TH2D*>::iterator end = fSpectralFunctionMap.end();
  for (std::map<int, TH2D*>::iterator iter = begin; iter != end; ++iter) {
    TH2D* hist = iter->second;
    if ( hist ) delete hist;
  }
}
//____________________________________________________________________________
bool SpectralFunc::GenerateNucleon(const Target& target) const
{
  TH2D* sf = this->SelectSpectralFunction( target );

  if ( !sf ) {
    fCurrRemovalEnergy = 0.;
    fCurrMomentum.SetXYZ(0., 0., 0.);
    return false;
  }

  double kmin    = sf->GetXaxis()->GetBinLowEdge( 1 ); // momentum range
  double kmax    = sf->GetXaxis()->GetBinLowEdge( fNumXBins + 1 );
  double wmin    = sf->GetYaxis()->GetBinLowEdge( 1 ); // removal energy range
  double wmax    = sf->GetYaxis()->GetBinLowEdge( fNumYBins + 1 );

  LOG("SpectralFunc", pINFO) << "Momentum range = ["   << kmin << ", " << kmax << "]";
  LOG("SpectralFunc", pINFO) << "Rmv energy range = [" << wmin << ", " << wmax << "]";

  // Temporarily replace ROOT's gRandom RNG with GENIE's so that our call to
  // TH2::GetRandom2() will use GENIE's random numbers
  TRandom* old_gRandom = gRandom;
  RandomGen* rnd = RandomGen::Instance();
  gRandom = &rnd->RndGen();

  // Pick a random nucleon momentum and removal energy using the spectral
  // function histogram
  double kc = 0.;
  double wc = 0.;
  sf->GetRandom2( kc, wc );

  // Restore ROOT's old gRandom now that we've had our fun
  gRandom = old_gRandom;

  // Log our choices for the nucleon momentum and removal energy
  LOG("SpectralFunc", pINFO) << "|p,nucleon| = " << kc;
  LOG("SpectralFunc", pINFO) << "|w,nucleon| = " << wc;

  // Generate momentum components
  double costheta = -1. + 2. * rnd->RndGen().Rndm();
  double sintheta = TMath::Sqrt(1.-costheta*costheta);
  double fi       = 2 * kPi * rnd->RndGen().Rndm();
  double cosfi    = TMath::Cos(fi);
  double sinfi    = TMath::Sin(fi);

  double kx = kc * sintheta * cosfi;
  double ky = kc * sintheta * sinfi;
  double kz = kc * costheta;

  // set generated values
  fCurrRemovalEnergy = wc;
  fCurrMomentum.SetXYZ(kx,ky,kz);

  return true;
}
//____________________________________________________________________________
double SpectralFunc::Prob(double p, double w, const Target & target) const
{
  TH2D* sf = this->SelectSpectralFunction( target );
  if ( !sf ) return 0.;

  // TODO: check whether we should return the probability mass or density here.
  // The previous implementation used the density. For now, we do that too.
  int bin_num = sf->FindBin( p, w );
  double prob_mass = sf->GetBinContent( bin_num );

  int bx, by, bz;
  sf->GetBinXYZ( bin_num, bx, by, bz );

  double prob_density = prob_mass / sf->GetXaxis()->GetBinWidth( bx )
    / sf->GetYaxis()->GetBinWidth( by );

  return prob_density;
}
//____________________________________________________________________________
void SpectralFunc::Configure(const Registry & config)
{
  Algorithm::Configure( config );
  this->LoadConfig();
}
//____________________________________________________________________________
void SpectralFunc::Configure(string param_set)
{
  Algorithm::Configure( param_set );
  this->LoadConfig();
}
//____________________________________________________________________________
void SpectralFunc::LoadConfig(void)
{
  // Determine the folder that should be searched to find the spectral function
  // data files
  std::string data_path_type;
  this->GetParamDef( "DataPathType", data_path_type, std::string("relative") );

  std::string data_path;
  this->GetParam( "DataPath", data_path );
  // Ensure that the path ends with '/', just in case
  data_path += '/';

  if ( data_path_type == "relative" ) {
    std::string genie_path = gSystem->Getenv("GENIE");
    genie_path += '/';
    data_path = genie_path + data_path;
  }

  fDataPath = data_path;

  // Get the number of bins to use for each axis of the spectral function
  // histogram
  this->GetParamDef( "NumMomentumBins", fNumXBins, 1000 );
  this->GetParamDef( "NumRemovalEnergyBins", fNumYBins, 1000 );

}
//____________________________________________________________________________
TGraph2D* SpectralFunc::Convert2Graph(TNtupleD& sfdata) const
{
  int np = sfdata.GetEntries();
  TGraph2D * sfgraph = new TGraph2D( np );

  sfdata.Draw("k:e:prob","","GOFF");
  assert( np == sfdata.GetSelectedRows() );
  double* k = sfdata.GetV1();
  double* e = sfdata.GetV2();
  double* p = sfdata.GetV3();

  for (int i = 0; i < np; ++i ) {
    double ki = k[i] * ( units::MeV / units::GeV ); // momentum
    double ei = e[i] * ( units::MeV / units::GeV ); // removal energy
    double pi = p[i] * TMath::Power(ki,2); // probabillity
    sfgraph->SetPoint(i, ki, ei, pi);
  }
  return sfgraph;
}
//____________________________________________________________________________
TH2D* SpectralFunc::SelectSpectralFunction(const Target& t) const
{
  // First check whether we've already built the requested spectral function
  int target_pdg = t.Pdg();
  if ( fSpectralFunctionMap.count(target_pdg) ) {
    return fSpectralFunctionMap.at( target_pdg );
  }

  // If not, attempt to build it
  std::string target_pdg_string = replace_with_std_to_string( target_pdg );
  RgKey sf_key( "SpectFuncTable@Pdg=" + target_pdg_string );
  std::string data_filename;
  this->GetParamDef( sf_key, data_filename, std::string() );

  if ( data_filename.empty() ) {
    LOG("SpectralFunc", pERROR) << "** The spectral function for target "
      << target_pdg << " isn't available";
    std::exit( 1 );
  }

  // Prepend the data path to the file name
  data_filename = fDataPath + data_filename;

  TNtupleD temp_sf_data( ("sfdata_" + target_pdg_string).c_str(), "",
    "k:e:prob" );

  temp_sf_data.ReadFile( data_filename.c_str() );

  LOG("SpectralFunc", pDEBUG) << "Loaded " << temp_sf_data.GetEntries()
    << " spectral function data points for pdg code " << target_pdg;

  TGraph2D* temp_sf_graph = this->Convert2Graph( temp_sf_data );
  temp_sf_graph->SetName( ("sf_" + target_pdg_string).c_str() );

  // Create a TF2 that we can use to easily integrate the TGraph2D over
  // histogram bins
  double x_min = temp_sf_graph->GetXmin();
  double x_max = temp_sf_graph->GetXmax();

  double y_min = temp_sf_graph->GetYmin();
  double y_max = temp_sf_graph->GetYmax();

  FunctorThatShouldBeALambda functor( temp_sf_graph );
  TF2 temp_sf_func("temp_sf_func", functor, x_min, x_max, y_min, y_max, 0);

  // Create a 2D histogram to represent the spectral function
  TH2D* sf_hist = new TH2D( ("sf_hist" + target_pdg_string).c_str(),
    ("spectral function for pdg code " + target_pdg_string
    + "; nucleon momentum (GeV); removal energy (GeV); probability").c_str(),
    fNumXBins, x_min, x_max, fNumYBins, y_min, y_max );

  // Set the directory to NULL so that this histogram is never auto-deleted
  // (the genie::SpectralFunc object will take ownership)
  sf_hist->SetDirectory( NULL );

  for ( int bx = 1; bx <= fNumXBins; ++bx ) {
    for ( int by = 1; by <= fNumYBins; ++by ) {

      // Get the current bin's boundaries
      double bin_x_min = sf_hist->GetXaxis()->GetBinLowEdge( bx );
      double bin_x_max = sf_hist->GetXaxis()->GetBinLowEdge( bx + 1 );

      double bin_y_min = sf_hist->GetYaxis()->GetBinLowEdge( by );
      double bin_y_max = sf_hist->GetYaxis()->GetBinLowEdge( by + 1 );

      // Probability mass for the current bin
      double pmf = temp_sf_func.Integral( bin_x_min, bin_x_max, bin_y_min,
        bin_y_max );

      // Get the global bin number for the current bin
      int bin_num = sf_hist->GetBin( bx, by );

      // Set its contents to match the probability mass function
      sf_hist->SetBinContent( bin_num, pmf );
    }
  }

  delete temp_sf_graph;

  // Store the new spectral function histogram for easy retrieval later
  fSpectralFunctionMap[ target_pdg ] = sf_hist;

  return sf_hist;
}
