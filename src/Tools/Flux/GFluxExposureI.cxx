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

#include "TMath.h"
#include <cstring>

namespace genie {
namespace flux {

  const char* GFluxExposureI::AsString(genie::flux::Exposure_t etype)
  {
    switch (etype) {
    case kUnknown:  return "UnknownExposureUnits";   break;
    case kPOTs:     return "POT";                    break;
    case kSeconds:  return "Seconds";                break;
    default:        return "?UnknownExposureUnits?"; break;
    }
  }

    genie::flux::Exposure_t 
      GFluxExposureI::StringToEnum(const char* chars, int maxChar)
  {
    int len = strlen(chars);
    if (maxChar == 0  ) maxChar = len;
    if (maxChar >  len) maxChar = len;
    if (maxChar <= 0  ) return genie::flux::kUnknown;
    
    char* lowchars = new char [maxChar];
    for (int i=0; i<maxChar; ++i) lowchars[i] = tolower(chars[i]);
    
    genie::flux::Exposure_t etype = genie::flux::kUnknown;
    if        ( 0 == strncmp(lowchars,"pot",TMath::Min(maxChar,3))) {
      etype = genie::flux::kPOTs;
    } else if ( 0 == strncmp(lowchars,"sec",TMath::Min(maxChar,3))) {
      etype = genie::flux::kSeconds;
    }
    delete [] lowchars;
    return etype;
  }

  GFluxExposureI::GFluxExposureI(genie::flux::Exposure_t etype) 
  { fEType = etype; }

  GFluxExposureI::~GFluxExposureI() { ; }

  genie::flux::Exposure_t GFluxExposureI::GetExposureType() const
  { return fEType; }

  const char* GFluxExposureI::GetExposureUnits() const
  { return genie::flux::GFluxExposureI::AsString(fEType); }

} // namespace flux
} // namespace genie
