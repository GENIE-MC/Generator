//____________________________________________________________________________
/*!
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Robert Hatcher <rhatcher@fnal.gov>
 Fermi National Accelerator Laboratory
*/
//____________________________________________________________________________

#include "Tools/Flux/GFluxExposureI.h"

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
