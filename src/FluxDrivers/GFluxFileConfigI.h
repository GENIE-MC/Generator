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

#include "PDG/PDGCodeList.h"
class TTree;

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

    /// allow caller to copy current status / ntuple entry info
    /// in the output file by providing copies of internal info
    ///
    /// Assumes that branch object pointers will not change
    /// which may require either a copy be made or, if using the class
    /// directly for reading the branch, one must force ROOT to
    /// not autodelete: 
    ///   myns::MyClassType* fCurrMyClass = new myns::MyClassType;
    ///   myTree->SetBranchAddress("bname",&fCurMyClass);
    ///   //? TBranch* b = myTree->GetBranch("bname");
    ///   //? b->SetAutoDelete(false);
    ///
    /// ensure vectors are sized sufficiently (or use .push_back())
    ///  branchNames[i]       = "bname"
    ///  branchClassNames[i]  = "myns::MyClassType"
    ///  branchObjPointers[i] = (void**)&fCurMyClass;

    virtual void  GetBranchInfo(std::vector<std::string>& branchNames,
                                std::vector<std::string>& branchClassNames,
                                std::vector<void**>&      branchObjPointers);

    virtual TTree* GetMetaDataTree(); //

    /// print the current configuration
    virtual void         PrintConfig() = 0;

    /// specify list of flux neutrino species
    virtual void         SetFluxParticles(const PDGCodeList & particles);

    /// set flux neutrino initial z position (upstream of the detector)
    /// pushed back from the normal flux window
    virtual void         SetUpstreamZ(double z0);

    /// limit cycling through input files
    virtual void         SetNumOfCycles(long int ncycle);

  protected:  // visible to derived classes

    PDGCodeList * fPdgCList;     ///< list of neutrino pdg-codes to generate  
    PDGCodeList * fPdgCListRej;  ///< list of nu pdg-codes seen but rejected
    std::string   fXMLbasename;  ///< XML file that might hold config param_sets
    long int      fNCycles;      ///< # times to cycle through the ntuple(s)
    long int      fICycle;       ///< current file cycle
                                 ///< default 0 = infinitely
    double        fZ0;           ///< configurable starting z position for 
                                 ///< each flux neutrino (in detector coord system)
  };

} // namespace flux
} // namespace genie

#endif //GENIE_FLUX_GFLUXFILECONFIGI_H
