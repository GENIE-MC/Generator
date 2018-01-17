////////////////////////////////////////////////////////////////////////
/// \file  GFlavorMixerI.h
/// \class genie::flux::GFlavorMixerI
/// \brief GENIE interface for flavor modification
///
///        Specific implementations of this class when used in conjuction
///        with the genie::flux::GFluxBlender class allow it to act as
///        an intermediate between a concrete flux generator and the
///        genie::GMCJDriver class.  Using this adapter allows one to
///        apply neutrino flavor changes using different models without
///        modifying either the concrete flux generator or GMCJDriver.
///
///        Concrete instances of this interface must be configurable
///        from a string, and provide a means of calculating the
///        transition probability from one flavor to any of the PDG
///        codes { 0, 12, 14, 16, -12, -14, -16 } where 0 represents
///        the complete disappearance (decay, sterile, ...).
///
///        Probability is expected to be normalized (that is, the sum
///        of all possible outcomes, including 0, must be 1).
///
/// \author  Robert Hatcher <rhatcher \at fnal.gov>
///          Fermi National Accelerator Laboratory
///
/// \created 2010-10-31
/// \version $Id: GFlavorMixerI.h,v 1.1.1.1 2010/12/22 16:18:52 p-nusoftart Exp $
////////////////////////////////////////////////////////////////////////

#ifndef GENIE_FLUX_GFLAVORMIXERI_H
#define GENIE_FLUX_GFLAVORMIXERI_H

#include <string>

namespace genie {
namespace flux {

  class GFlavorMixerI {
    
  public:
  
    GFlavorMixerI();
    virtual ~GFlavorMixerI();

    //
    // define the GFlavorMixerI interface:
    //

    /// each schema must take a string that configures it
    /// it is up to the individual model to parse said string
    /// and extract parameters (e.g. sin2th23, deltam12, etc)
    virtual void      Config(std::string config) = 0;

    /// for any pair of PDG codes the model must calculate
    /// the transition probability.  This can also depend on
    /// neutrino energy (in GeV) and distance (in meters) from 
    /// the neutrino origin.
    virtual double    Probability(int pdg_initial, int pdg_final, 
                                  double energy, double dist) = 0;

    /// provide a means of printing the configuration
    virtual void     PrintConfig(bool verbose=true) = 0;

  };

} // namespace flux
} // namespace genie

#endif //GENIE_FLUX_GFLAVORMIXERI_H
