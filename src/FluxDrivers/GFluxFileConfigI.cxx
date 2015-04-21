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
#include "Messenger/Messenger.h"
#include "TMath.h"

namespace genie {
namespace flux {

  GFluxFileConfigI::GFluxFileConfigI()
    : fXMLbasename(""), fNCycles(0), fZ0(-3.4e38)
  { ; }

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

  //___________________________________________________________________________
  void GFluxFileConfigI::GetBranchInfo(std::vector<std::string>& branchNames,
                                       std::vector<std::string>& branchClassNames,
                                       std::vector<void**>&      branchObjPointers)
  {
    // allow flux driver to report back current status and/or ntuple entry 
    // info for possible recording in the output file by supplying
    // the class name, and a pointer to the object that will be filled
    // as well as a suggested name for the branch.

    // default is not to supply anything
  }
  TTree* GFluxFileConfigI::GetMetaDataTree()
  {
    return 0;
  }

  //___________________________________________________________________________
  void GFluxFileConfigI::SetUpstreamZ(double z0)
  {
    // The flux neutrino position (x,y) is given on the user specified 
    // flux window.  This method sets the preferred user coord starting z 
    // position upstream of detector face. Each flux neutrino will be 
    // backtracked from the initial flux window to the input z0.  
    // If the value is unreasonable (> 10^30) then the ray is left on 
    // the flux window.

    fZ0 = z0;
  }
  //___________________________________________________________________________
  void GFluxFileConfigI::SetNumOfCycles(long int ncycle)
  {
    // The flux ntuples can be recycled for a number of times to boost 
    // generated event statistics without requiring enormous beam simulation
    // statistics.
    // That option determines how many times the driver is going to cycle 
    // through the input flux ntuple.
    // With ncycle=0 the flux ntuple will be recycled an infinite amount of 
    // times so that the event generation loop can exit only on a POT or 
    // event num check.

    fNCycles = TMath::Max(0L, ncycle);
  }
  //___________________________________________________________________________
  void GFluxFileConfigI::SetFluxParticles(const PDGCodeList & particles)
  {
    if (!fPdgCList) {
      fPdgCList = new PDGCodeList;
    }
    fPdgCList->Copy(particles);
    
    LOG("Flux", pINFO)
      << "Declared list of neutrino species: " << *fPdgCList;
  }

} // namespace flux
} // namespace genie
