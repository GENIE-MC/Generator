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
#include "Physics/Multinucleon/XSection/ValenciaHadronTensorI.h"

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
  const ValenciaHadronTensorI* GetTensor(int tensor_pdg,
    HadronTensorType_t type);

  private:

  HadronTensorPool();
  HadronTensorPool(const HadronTensorPool& htp);
  virtual ~HadronTensorPool();

  std::string FindTensorTableFile(const std::string& basename,
    bool& ok) const;

  bool LoadConfig(void);
  XmlParserStatus_t ParseXMLConfig(const std::string& filename,
    const std::string& table_to_use = "Default");

  // Keys are (PDG code, hadron tensor type) pairs, values are pointers
  // to hadron tensor objects
  std::map< std::pair<int, HadronTensorType_t>,
    ValenciaHadronTensorI* > fTensors;

  // Paths to check when searching for hadron tensor data files
  std::vector<std::string> fDataPaths;
};

} // genie namespace

#endif // _HADRON_TENSOR_POOL_H_
