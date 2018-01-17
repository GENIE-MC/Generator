////////////////////////////////////////////////////////////////////////
/// \file  GFlavorMap.h
/// \class genie::flux::GFlavorMap
/// \brief GENIE interface for flavor modification
///
///        Concrete instance of GFlavorMixerI that maps from
///        one flavor to another independent of energy or distance.
///        Users specify the transition probability from one flavor 
///        to any of the PDG codes { 0, 12, 14, 16, -12, -14, -16 } 
///        where 0 represents the complete disappearance (decay, sterile, ...).
///
///        Probability is expected to be normalized (that is, the sum
///        of all possible outcomes, including 0, must be 1).
///
///        Supported config string formats:
///        1)  " swap pdg1:pdg2  pdg3:pdg4 "
///            Map all neutrinos of flavor "pdg1" to "pdg2",  "pdg3" to "pdg4"
///            Use PDG values { 0, 12, 14, 16, -12, -14, -16 }
///              for { sterile, nu_e, nu_mu, nu_tau, nu_e_bar, ...}
///            Unnamed initial flavors are left unmodified.
///            Use numeric values with spaces only between pairs.
///            Must start with the literal "swap"
///            (FMWK note:  param must be surrounded by explicit "'s)
///        2)  " fixedfrac {pdg1:f0,f12,f14,f16,f-12,f-14,f-16} ..."
///            For each group delineated by {}'s the map the "pdg1"
///            by each pdg by the fraction given [0-1, sum=1].
///            So {12:0.5,0.5,0,0,0,0,0} means nu_e => 50/50% nu_e/nu_mu.
///            Each list *must* have an int + 7 fractions.
///
/// \author  Robert Hatcher <rhatcher \at fnal.gov>
///          Fermi National Accelerator Laboratory
///
/// \created 2010-10-31
/// \version $Id: GFlavorMap.h,v 1.1.1.1 2010/12/22 16:18:52 p-nusoftart Exp $
////////////////////////////////////////////////////////////////////////

#ifndef GENIE_FLUX_GFLAVORSWAP_H
#define GENIE_FLUX_GFLAVORSWAP_H

#include <string>
#include "FluxDrivers/GFlavorMixerI.h"

namespace genie {
namespace flux {

  class GFlavorMap : public GFlavorMixerI {
    
  public:
  
    GFlavorMap();
    ~GFlavorMap();

    //
    // implement the GFlavorMixerI interface:
    //

    /// each schema must take a string that configures it
    /// it is up to the individual model to parse said string
    /// and extract parameters (e.g. sin2th23, deltam12, etc)
    void      Config(std::string config);

    /// for any pair of PDG codes the model must calculate
    /// the transition probability.  This can also depend on
    /// neutrino energy (in GeV) and distance (in meters) from 
    /// the neutrino origin.
    double    Probability(int pdg_initial, int pdg_final, 
                          double energy, double dist);

    /// provide a means of printing the configuration
    void     PrintConfig(bool verbose=true);

  private:

    void         ParseMapString(std::string config);
    void         ParseFixedfracString(std::string config);

    int          PDG2Indx(int pdg);
    int          Indx2PDG(int indx);
    const char*  IndxName(int indx);
    const char*  NuName(int pdg) { return IndxName(PDG2Indx(pdg)); }

    double   fProb[7][7];

  };

} // namespace flux
} // namespace genie

//
//    Name        PDG   Indx
//    sterile       0   0
//    nu_e         12   1
//    nu_mu        14   2
//    nu_tau       16   3
//    nu_e_bar    -12   4
//    nu_mu_bar   -14   5
//    nu_tau_bar  -16   6
//
inline int genie::flux::GFlavorMap::PDG2Indx(int pdg)
{
  switch ( pdg ) {
  case  12: return 1; break;
  case  14: return 2; break;
  case  16: return 3; break;
  case -12: return 4; break;
  case -14: return 5; break;
  case -16: return 6; break;
  default:  return 0; break;
  }
  return 0;
}
inline int genie::flux::GFlavorMap::Indx2PDG(int indx)
{
  switch ( indx ) {
  case  1: return  12; break;
  case  2: return  14; break;
  case  3: return  16; break;
  case  4: return -12; break;
  case  5: return -14; break;
  case  6: return -16; break;
  default: return   0; break;
  }
  return 0;
}

#endif //GENIE_FLUX_GFLAVORSWAP_H
