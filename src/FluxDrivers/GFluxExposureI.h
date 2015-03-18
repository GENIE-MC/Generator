////////////////////////////////////////////////////////////////////////
/// \file  GFluxExposureI.h
/// \class genie::flux::GFluxExposureI
/// \brief GENIE interface for uniform flux exposure iterface
///
///        Unified flux exposure interface to be used by flux drivers
///        that can support such.
///
/// \author  Robert Hatcher <rhatcher \at fnal.gov>
///          Fermi National Accelerator Laboratory
///
/// \created 2015-03-17  initial version
/// \version $Id: $
////////////////////////////////////////////////////////////////////////

#ifndef GENIE_FLUX_GFLUXEXPOSUREI_H
#define GENIE_FLUX_GFLUXEXPOSUREI_H

#include <string>

namespace genie {
namespace flux {

  typedef enum EExposure {
    kUnknown   = 0,
    kPOTs      = 1,  // particles (protons) on target
    kSeconds   = 2,  // exposure duration
    kNExposureTypes  //
  } Exposure_t;

  class GFluxExposureI {
    
  public:
  
    GFluxExposureI(genie::flux::Exposure_t etype);
    virtual ~GFluxExposureI();

    genie::flux::Exposure_t GetExposureType() const;

    //
    // define the GFluxExposureI interface:
    //
    virtual double GetTotalExposure() const = 0;

  private:
    genie::flux::Exposure_t fEType; 

  };

} // namespace flux
} // namespace genie

#endif //GENIE_FLUX_GFLUXEXPOSUREI_H
