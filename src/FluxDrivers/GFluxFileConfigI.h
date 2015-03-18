////////////////////////////////////////////////////////////////////////
/// \file  GFluxFileConfigI.h
/// \class genie::flux::GFluxFileConfigI
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

#ifndef GENIE_FLUX_GFLUXFILECONFIGI_H
#define GENIE_FLUX_GFLUXFILECONFIGI_H

#include <string>
#include <vector>
#include <set>

namespace genie {
namespace flux {

  class GFluxFileConfigI {
    
  public:
  
    GFluxFileConfigI();
    virtual ~GFluxFileConfigI();

    //
    // define the GFluxFileConfigI interface:
    //
    /// first is primary method for loading root flux ntuple files and config
    /// others are alternatives that can be overloaded but have
    /// sensible defaults to fall back to calling the vector version

    virtual void  LoadBeamSimData(const std::vector<std::string>& filenames,
                                  const std::string&              det_loc) = 0;

    virtual void  LoadBeamSimData(const std::set<std::string>&    filenames,
                                  const std::string&              det_loc);

    virtual void  LoadBeamSimData(const std::string&              filename, 
                                  const std::string&              det_loc);

    virtual void         SetXMLFileBase(std::string xmlbasename="");
    virtual std::string  GetXMLFileBase() const { return fXMLbasename; }

  private:
  
    std::string   fXMLbasename;  ///< XML file that might hold config param_sets

  };

} // namespace flux
} // namespace genie

#endif //GENIE_FLUX_GFLUXFILECONFIGI_H
