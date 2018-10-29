//____________________________________________________________________________
/*!

\class    genie::HadronTensorPool

\brief    Singleton class to load and serve hadron tensor objects for cross
          section calculations

\author   Steven Gardiner <gardiner \at fnal.gov>
          Fermi National Accelerator Laboratory

\created  August 28, 2018

\cpright  Copyright (c) 2003-2018, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _HADRON_TENSOR_POOL_H_
#define _HADRON_TENSOR_POOL_H_

// standard library includes
#include <map>
#include <string>
#include <vector>

// GENIE includes
#include "Framework/Conventions/XmlParserStatus.h"
#include "Physics/HadronTensors/HadronTensorI.h"

namespace genie {

class HadronTensorPool {
  public:

  static HadronTensorPool& Instance(void);

  /// Retrieves a pointer to a hadron tensor object from the pool
  /// \param[in] tensor_pdg The PDG code for the nuclide described by the
  /// tensor
  /// \param[in] type The desired kind of hadron tensor
  /// \returns A pointer to the requested hadron tensor, or NULL if a match
  /// could not be found in the pool
  const HadronTensorI* GetTensor(int tensor_pdg, HadronTensorType_t type);

  private:

  HadronTensorPool();
  HadronTensorPool(const HadronTensorPool& htp);
  virtual ~HadronTensorPool();

  /// Looks up the full path when constructing hadron tensor objects that are
  /// based on a data file
  std::string FindTensorTableFile(const std::string& basename,
    bool& ok) const;

  bool LoadConfig(void);
  XmlParserStatus_t ParseXMLConfig(const std::string& filename,
    const std::string& table_to_use = "Default");

  /// Cache of hadron tensor objects that have been fully loaded into memory
  ///
  /// Keys are (PDG code, hadron tensor type) pairs, values are pointers
  /// to hadron tensor objects
  std::map< std::pair<int, HadronTensorType_t>, HadronTensorI* > fTensors;

  /// Paths to check when searching for hadron tensor data files
  std::vector<std::string> fDataPaths;

  /// Struct that stores the XML attributes describing each hadron tensor.
  ///
  /// It is used to enable lazy initialization of hadron tensor objects by the
  /// HadronTensorPool. There should be one std::string in this struct for each
  /// recognized attribute of the <tensor> XML tag in HadronTensors.xml except
  /// for "type", which will be automatically converted to a
  /// genie::HadronTensorType_t while parsing the configuration file.
  struct HadronTensorXMLAttributes {
    HadronTensorXMLAttributes() {}
    HadronTensorXMLAttributes(const std::string& the_calc,
      const std::string& the_file) : calc(the_calc), file(the_file) {}
    std::string calc;
    std::string file;
  };

  /// XML attributes for all known hadron tensors. Used to construct
  /// them on demand (lazy initialization)
  ///
  /// Keys are (PDG code, hadron tensor type) pairs, values are structs
  /// containing the XML attributes describing the tensor in the
  /// HadronTensors.xml configuration file
  std::map< std::pair<int, HadronTensorType_t>, HadronTensorXMLAttributes >
    fTensorAttributes;

  /// Build a hadron tensor object on demand using a set of XML attributes
  const genie::HadronTensorI* BuildTensor(
    const std::pair<int, genie::HadronTensorType_t>& tensor_id,
    const HadronTensorXMLAttributes& attributes);

};

} // genie namespace

#endif // _HADRON_TENSOR_POOL_H_
