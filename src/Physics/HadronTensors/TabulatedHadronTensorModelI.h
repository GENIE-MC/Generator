//____________________________________________________________________________
/*!

\class    genie::TabulatedHadronTensorModelI

\brief    Creates hadron tensor objects for cross section calculations
          using precomputed data tables

\author   Steven Gardiner <gardiner \at fnal.gov>
          Liang Liu <liangliu \at fnal.gov>
          Fermi National Accelerator Laboratory

\created  April 26, 2019

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE

*/
//____________________________________________________________________________

#ifndef _TABULATED_HADRON_TENSOR_MODEL_H_
#define _TABULATED_HADRON_TENSOR_MODEL_H_

// standard library includes
#include <map>
#include <string>
#include <vector>

// GENIE includes
#include "Framework/Algorithm/Algorithm.h"
#include "Physics/HadronTensors/HadronTensorI.h"
#include "Physics/HadronTensors/HadronTensorModelI.h"

namespace genie {

class TabulatedHadronTensorModelI : public HadronTensorModelI {

public:
  virtual ~TabulatedHadronTensorModelI();

  // Implementation of Algorithm interface
  virtual void Configure(const Registry& config);
  virtual void Configure(std::string config);

  // Implementation of HadronTensorModelI interface
  virtual const HadronTensorI* GetTensor(int tensor_pdg, HadronTensorType_t type) const;

protected:

  TabulatedHadronTensorModelI();
  TabulatedHadronTensorModelI(std::string name);
  TabulatedHadronTensorModelI(std::string name, std::string config);

  /// Saves some basic XML config parameters to data members
  void LoadConfig();

  /// Looks up the full path when constructing hadron tensor objects that are
  /// based on a data file
  std::string FindTensorTableFile(const std::string& basename,
    bool& ok) const;

  /// Struct used to provide a unique ID for each tensor object
  struct HadronTensorID {
    HadronTensorID(int pdg = 0, HadronTensorType_t typ = kHT_Undefined)
      : target_pdg( pdg ), type( typ ) {}
    int target_pdg;
    HadronTensorType_t type;

    // Less than operator needed for sorting a map of these IDs
    bool operator<(const HadronTensorID& other) const {
      return (target_pdg < other.target_pdg)
        || (target_pdg == other.target_pdg && type < other.type);
    }
  };

  /// If true, logging messages will be issued when a requested hadron tensor
  /// file cannot be found
  bool fWarnIfMissing;

  /// Cache of hadron tensor objects that have been fully loaded into memory
  ///
  /// Keys are tensor IDs, values are pointers to hadron tensor objects
  mutable std::map< HadronTensorID, HadronTensorI* > fTensors;

  /// Paths to check when searching for hadron tensor data files
  std::vector<std::string> fDataPaths;

  /// Create a HadronTensorI object given a particular HadronTensorID
  const HadronTensorI* BuildTensor( const HadronTensorID& ht_id ) const;

  /// Loads the basename for a particular hadron tensor file from
  /// the configuration Registry
  std::string GetTensorFileBasename( const HadronTensorID& ht_id ) const;

  /// Parses the hadron tensor file (specified by its full file name,
  /// including the path) and returns a HadronTensorI* to it
  virtual HadronTensorI* ParseTensorFile( const std::string& full_file_name ) const = 0;
};

} // namespace genie

#endif // _TABULATED_HADRON_TENSOR_MODEL_H_
