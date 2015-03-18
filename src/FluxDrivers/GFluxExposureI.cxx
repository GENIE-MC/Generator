////////////////////////////////////////////////////////////////////////
/// \file  GFluxExposureI.cxx
/// \brief GENIE interface for uniform flux exposure iterface
///
/// \version $Id: $
/// \author  Robert Hatcher <rhatcher \at fnal.gov>
///          Fermi National Accelerator Laboratory
///
/// \update  2015-03-17 initial version
////////////////////////////////////////////////////////////////////////

#include "FluxDrivers/GFluxExposureI.h"

namespace genie {
namespace flux {

  GFluxExposureI::GFluxExposureI(genie::flux::Exposure_t etype) 
  { fEType = etype; }

  GFluxExposureI::~GFluxExposureI() { ; }

  genie::flux::Exposure_t GFluxExposureI::GetExposureType() const
  { return fEType; }

} // namespace flux
} // namespace genie
