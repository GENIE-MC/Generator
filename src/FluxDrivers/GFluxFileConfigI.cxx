////////////////////////////////////////////////////////////////////////
/// \file  GFluxFileConfigI.cxx
/// \brief GENIE interface for uniform flux exposure iterface
///
/// \version $Id: $
/// \author  Robert Hatcher <rhatcher \at fnal.gov>
///          Fermi National Accelerator Laboratory
///
/// \update  2015-03-17 initial version
////////////////////////////////////////////////////////////////////////

#include "FluxDrivers/GFluxFileConfigI.h"

namespace genie {
namespace flux {

  GFluxFileConfigI::GFluxFileConfigI() { ; }

  GFluxFileConfigI::~GFluxFileConfigI() { ; }

  void GFluxFileConfigI::SetXMLFileBase(std::string xmlbasename)
  { fXMLbasename = xmlbasename; }

  //___________________________________________________________________________
  void GFluxFileConfigI::LoadBeamSimData(const std::set<std::string>& fileset, 
                                         const std::string&           config )
  {
    // Loads a beam simulation root file into the GFluxFileConfig driver.
    // have a set<> want a vector<>
    std::vector<std::string> filevec;
    std::copy(fileset.begin(),fileset.end(),std::back_inserter(filevec));
    LoadBeamSimData(filevec,config); // call the one that takes a vector
  }

  //___________________________________________________________________________
  void GFluxFileConfigI::LoadBeamSimData(const std::string& filename, 
                                         const std::string& config )
  {
    // Loads a beam simulation root file into the GFluxFileConfig driver.
    std::vector<std::string> filevec;
    filevec.push_back(filename);
    LoadBeamSimData(filevec,config); // call the one that takes a vector
  }

} // namespace flux
} // namespace genie
