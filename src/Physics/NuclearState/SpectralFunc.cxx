//____________________________________________________________________________
/*
 Copyright (c) 2003-2020, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Steven Gardiner <gardiner \at fnal.gov>
 Fermi National Accelerator Laboratory

 Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
 University of Liverpool & STFC Rutherford Appleton Laboratory

 For the class documentation see the corresponding header file.
*/
//____________________________________________________________________________

#include <fstream>
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

namespace {

  // TODO: Replace this with std::to_string when we switch to C++11
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

  double kmin    = sf->GetXaxis()->GetXmin(); // momentum range
  double kmax    = sf->GetXaxis()->GetXmax();
  double wmin    = sf->GetYaxis()->GetXmin(); // removal energy range
  double wmax    = sf->GetYaxis()->GetXmax();

  LOG("SpectralFunc", pDEBUG) << "Momentum range = ["   << kmin << ", " << kmax << "]";
  LOG("SpectralFunc", pDEBUG) << "Rmv energy range = [" << wmin << ", " << wmax << "]";

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
  // The previous implementation used the density. For now, we return the mass.
  double prob_mass = 0.;
  if ( w >= 0. ) {
    // In the normal case, find the bin corresponding to (p, w) and return the
    // probability based on that
    int bin_num = sf->FindBin( p, w );
    prob_mass = sf->GetBinContent( bin_num );
  }
  else {
    // For w < 0 (e.g., w = -1 used by QELEventGenerator), ignore w and find
    // the integrated probability of sampling from the appropriate p bin
    TH1D* temp_p_hist = sf->ProjectionX("temp_p_proj_x", 1,
      sf->GetXaxis()->GetNbins());

    int p_bin_num = temp_p_hist->FindBin( p );
    prob_mass = temp_p_hist->GetBinContent( p_bin_num );
    delete temp_p_hist;
  }

  return prob_mass;

  //int bx, by, bz;
  //sf->GetBinXYZ( bin_num, bx, by, bz );

  //double prob_density = prob_mass / sf->GetXaxis()->GetBinWidth( bx )
  //  / sf->GetYaxis()->GetBinWidth( by )
  //  / std::pow(sf->GetXaxis()->GetBinLowEdge( bx ), 2)
  //  / (4. * kPi);

  //return prob_density;
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
  NuclearModelI::LoadConfig();

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
}
//____________________________________________________________________________
TH2D* SpectralFunc::SelectSpectralFunction(const Target& t) const
{
  // First check whether we've already built the requested spectral function
  int target_pdg = t.Pdg();
  std::map<int, TH2D*>::iterator my_iter = fSpectralFunctionMap.find( target_pdg );
  if ( my_iter != fSpectralFunctionMap.end() ) return my_iter->second;

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

  TH2D* sf_hist = this->LoadSFDataFile( data_filename, t.Z() );

  LOG("SpectralFunc", pNOTICE) << "Loaded spectral function data"
    " for target with PDG code " << target_pdg << " from the file "
    << data_filename;

  sf_hist->SetName( ("sf_" + target_pdg_string).c_str() );

  // Set the directory to NULL so that this histogram is never auto-deleted
  // (the genie::SpectralFunc object will take ownership)
  sf_hist->SetDirectory( NULL );

  // Store the new spectral function histogram for easy retrieval later
  fSpectralFunctionMap[ target_pdg ] = sf_hist;

  return sf_hist;
}
//____________________________________________________________________________
TH2D* SpectralFunc::LoadSFDataFile(const std::string& full_file_name,
  int targetZ ) const
{
  int num_E_bins, num_p_bins;
  double p_min, p_max, E_min, E_max;

  std::ifstream in_file( full_file_name );

  // Check that the file exists and is readable
  if ( !in_file.good() ) {
    LOG("SpectralFunc", pERROR) << "Could not read spectral function table"
      << " from the file " << full_file_name;
    std::exit( 1 );
  }

  // TODO: add I/O error handling

  // Read the histogram bin boundary information from the
  // spectral function data file header. Note that the current
  // format requires equally-spaced bins.
  in_file >> num_E_bins >> num_p_bins;
  in_file >> E_min >> p_min;
  in_file >> E_max >> p_max;

  // The input files use MeV while GENIE uses GeV, so make the
  // change to GENIE units now.
  E_min *= genie::units::MeV;
  E_max *= genie::units::MeV;

  p_min *= genie::units::MeV;
  p_max *= genie::units::MeV;

  // Build an empty histogram with the given bin boundaries
  TH2D* sf_hist = new TH2D("temp_sf_hist",
    "Spectral Function; nucleon momentum (GeV); removal energy (GeV)",
    num_p_bins, p_min, p_max, num_E_bins, E_min, E_max);

  // Now set the contents of each bin using the body of the input data file.
  // Note that the tabulated p and E values represent the bin centers.
  for ( int ip = 0; ip < num_p_bins; ++ip ) {

    double p;
    in_file >> p;

    // Convert from MeV to GeV
    p *= genie::units::MeV;

    for ( int ie = 0; ie < num_E_bins; ++ie ) {

      double E, prob_density;
      in_file >> E >> prob_density;

      // Convert from MeV to GeV
      E *= genie::units::MeV;

      // Convert from MeV^(-4) to GeV^(-4)
      prob_density *= std::pow( genie::units::MeV, -4 );

      // Remove the normalization factor of Z from the values tabulated in the
      // file
      prob_density /= targetZ;

      // Convert bin contents from probability density to probability
      // mass for easier sampling
      TAxis* p_axis = sf_hist->GetXaxis();
      int p_idx = p_axis->FindBin( p );
      double p_bin_width = p_axis->GetBinWidth( p_idx );

      TAxis* E_axis = sf_hist->GetYaxis();
      int E_idx = E_axis->FindBin( E );
      double E_bin_width = E_axis->GetBinWidth( E_idx );

      // This expression assumes that the p^2 * prob_density can
      // be treated as approximately constant over the bin of interest
      double prob_mass = prob_density * std::pow(p, 2) * 4. * kPi
        * p_bin_width * E_bin_width;

      sf_hist->SetBinContent( sf_hist->FindBin(p, E), prob_mass );
      LOG("SpectralFunc", pDEBUG) << "p = " << p << ", E = " << E << ", dens = " << prob_density
        << ", mass = " << prob_mass << ", p_bin_width = " << p_bin_width << ", E_bin_width = "
        << E_bin_width;
    }
  }

  return sf_hist;
}
