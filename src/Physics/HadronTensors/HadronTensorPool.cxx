//____________________________________________________________________________
/*
 Copyright (c) 2003-2018, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Steven Gardiner <gardiner \at fnal.gov>
         Fermi National Accelerator Laboratory

 For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________

// standard library includes
#include <cstdlib>
#include <fstream>
#include <string>

// libxml2 includes
#include "libxml/parser.h"
#include "libxml/xmlmemory.h"

// ROOT includes
#include <TSystem.h>

// GENIE includes
#include "Framework/Messenger/Messenger.h"
#include "Physics/HadronTensors/HadronTensorPool.h"
#include "Physics/HadronTensors/TabulatedValenciaHadronTensor.h"
#include "Physics/HadronTensors/HadronTensorI.h"
#include "Framework/Utils/StringUtils.h"
#include "Framework/Utils/XmlParserUtils.h"

namespace {

  /// Helper function that retrieves an attribute of an xmlNodePtr and trims
  /// whitespace before returning it
  std::string get_trimmed_attribute(const xmlNodePtr& node,
    const std::string& name)
  {
    return genie::utils::str::TrimSpaces(
      genie::utils::xml::GetAttribute(node, name.c_str()));
  }

  /// Converts a string to a genie::HadronTensorType_t value. If the string
  /// does not correspond to a valid tensor type, kHT_Undefined
  /// is returned, and the ok flag is set to false
  genie::HadronTensorType_t string_to_tensor_type(const std::string& str,
    bool& ok)
  {
    if (str == "MEC_FullAll") return genie::HadronTensorType::kHT_MEC_FullAll;
    else if (str == "MEC_Fullpn")
      return genie::HadronTensorType::kHT_MEC_Fullpn;
    else if (str == "MEC_DeltaAll")
      return genie::HadronTensorType::kHT_MEC_DeltaAll;
    else if (str == "MEC_Deltapn")
      return genie::HadronTensorType::kHT_MEC_Deltapn;
    else if (str == "MEC_EM")
      return genie::HadronTensorType::kHT_MEC_EM;
    else if (str == "QE_EM")
      return genie::HadronTensorType::kHT_QE_EM;
    else {
      ok = false;
      return genie::HadronTensorType::kHT_Undefined;
    }
  }

  /// Returns true if a given file exists and is accessible, or false otherwise
  bool file_exists(const std::string& file_name) {
    return std::ifstream(file_name).good();
  }

  /// Helper function for comparing xmlNode names and string literals. It helps
  /// to reduce the number of ugly casts needed in the XML parsing code.
  bool name_equal(const xmlNodePtr& np, const std::string& str) {
    return xmlStrEqual( np->name,
      reinterpret_cast<const xmlChar*>(str.c_str()) );
  }

  /// Helper function that converts the content of an xmlNode to a std::string
  std::string get_node_content(const xmlNodePtr& np) {
    return reinterpret_cast<const char*>( xmlNodeGetContent(np) );
  }
}

//____________________________________________________________________________
genie::HadronTensorPool::HadronTensorPool()
{
  bool init_ok = LoadConfig();
  if ( !init_ok ) LOG("HadronTensorPool", pERROR) << "Failed to initialize"
    " the HadronTensorPool";
}

//____________________________________________________________________________
genie::HadronTensorPool::~HadronTensorPool()
{
  std::map< std::pair<int, genie::HadronTensorType_t>, HadronTensorI* >
    ::iterator it;
  for(it = fTensors.begin(); it != fTensors.end(); ++it) {
    HadronTensorI* t = it->second;
    if (t) delete t;
  }
  fTensors.clear();
}
//____________________________________________________________________________
genie::HadronTensorPool& genie::HadronTensorPool::Instance()
{
  static HadronTensorPool singleton_instance;
  return singleton_instance;
}

//____________________________________________________________________________
const genie::HadronTensorI* genie::HadronTensorPool::GetTensor(
  int tensor_pdg, genie::HadronTensorType_t type)
{
  std::pair<int, genie::HadronTensorType_t> temp_pair(tensor_pdg, type);

  // First check to see if the hadron tensor object already exists in memory.
  // If it does, return the existing object.
  if ( fTensors.count(temp_pair) ) return fTensors.at(temp_pair);
  // If not, check to see if we have XML attributes describing it.
  else if ( fTensorAttributes.count(temp_pair) ) {
    return BuildTensor( temp_pair, fTensorAttributes.at(temp_pair) );
  }

  // Otherwise, give up and return a null pointer
  return NULL;
}

//____________________________________________________________________________
bool genie::HadronTensorPool::LoadConfig(void)
{
  bool ok = true;

  // Find the XML configuration file
  std::string filename = genie::utils::xml::GetXMLFilePath("HadronTensors.xml");

  LOG("HadronTensorPool", pINFO)  << "Loading hadron tensor settings"
    << " from the file " << filename;

  if ( file_exists(filename) ) {
    XmlParserStatus_t status = this->ParseXMLConfig(filename);
    if (status != kXmlOK) {
      LOG("HadronTensorPool", pWARN) << "Error encountered while"
        << " attempting to parse the XMl file \"" << filename << "\"."
        << " XML parser status: " << XmlParserStatus::AsString(status);
      ok = false;
    }
  }
  else {
    LOG("HadronTensorPool", pWARN) << "Could not read from the file: "
      << filename;
    ok = false;
  }
  return ok;
};

//____________________________________________________________________________
std::string genie::HadronTensorPool::FindTensorTableFile(
  const std::string& basename, bool& ok) const
{
  for (size_t p = 0; p < fDataPaths.size(); ++p) {
    const std::string& path = fDataPaths.at(p);
    std::string full_name = path + '/' + basename;
    if ( file_exists(full_name) ) return full_name;
  }

  // A matching file could not be found
  ok = false;
  return std::string();
}

//____________________________________________________________________________
genie::XmlParserStatus_t genie::HadronTensorPool::ParseXMLConfig(
  const std::string& filename, const std::string& table_to_use)
{
  LOG("HadronTensorPool", pDEBUG) << "Reading XML file: " << filename;

  xmlDocPtr xml_doc = xmlParseFile( filename.c_str() );
  if ( !xml_doc ) return kXmlNotParsed;

  xmlNodePtr xml_root = xmlDocGetRootElement( xml_doc );
  if ( !xml_root )
  {
    xmlFreeDoc(xml_doc);
    return kXmlEmpty;
  }

  if ( !name_equal(xml_root, "hadron_tensor_config") )
  {
    LOG("HadronTensorPool", pERROR) << "Missing <hadron_tensor_config>"
      << " tag in the configuration file " << filename;
    xmlFreeDoc(xml_doc);
    return kXmlInvalidRoot;
  }

  xmlNodePtr xml_tensor_table = xml_root->xmlChildrenNode;

  // Flag that indicates whether the requested table of hadron tensors
  // could be found
  bool found_table = false;

  // loop over <tensor_table> and <data_paths> nodes
  while (xml_tensor_table) {

    if ( name_equal(xml_tensor_table, "tensor_table") ) {

      // Load only the requested tensor table (this allows for multiple tables
      // with different names to be placed in the same XML configuration file)
      std::string table_name = get_trimmed_attribute(xml_tensor_table, "name");
      if ( table_name != table_to_use ) {
        xml_tensor_table = xml_tensor_table->next;
        continue;
      }
      else found_table = true;

      // loop over the <data_paths> entries in the tensor table
      xmlNodePtr xml_data_paths = xml_tensor_table->xmlChildrenNode;
      while ( xml_data_paths ) {
        if ( name_equal(xml_data_paths, "data_paths") ) {
          xmlNodePtr xml_path = xml_data_paths->xmlChildrenNode;

          while ( xml_path ) {
            if ( name_equal(xml_path, "path") ) {
              std::string path = genie::utils::str::TrimSpaces(
                get_node_content(xml_path));

              std::string path_type = get_trimmed_attribute(xml_path, "type");
              if (path_type == "relative") {
                // paths may be specified relative to the $GENIE folder
                path = std::string( gSystem->Getenv("GENIE") ) + '/' + path;
              }

              LOG("HadronTensorPool", pINFO) << "The HadronTensorPool will"
                " search for data files in " << path;
              fDataPaths.push_back( path );
            } // <x> == <path>
            xml_path = xml_path->next;
          } // <path> loop
        } // <x> == <data_paths>
        xml_data_paths = xml_data_paths->next;
      } // <data_paths> loop

      xmlNodePtr xml_nuclide = xml_tensor_table->xmlChildrenNode;

      // loop over the <nuclide> entries in the tensor table
      while ( xml_nuclide ) {
        if ( name_equal(xml_nuclide, "nuclide") ) {

          std::string pdg_str = get_trimmed_attribute(xml_nuclide, "pdg");

          LOG("HadronTensorPool", pDEBUG) << "Reading hadron tensor"
            << " configuration for nuclide " << pdg_str;

          int pdg = std::atoi( pdg_str.c_str() );

          xmlNodePtr xml_tensor = xml_nuclide->xmlChildrenNode;

          std::string type_str("unknown");

          while (xml_tensor) {
            if ( name_equal(xml_tensor, "tensor") ) {

              // Flag to indicate if there was a problem processing
              // the record for the current tensor
              bool tensor_ok = true;

              type_str = get_trimmed_attribute(xml_tensor, "type");

              // Check that the declared hadron tensor type is recognized
              genie::HadronTensorType_t type
                = string_to_tensor_type(type_str, tensor_ok);

              std::string calc_str
                = get_trimmed_attribute(xml_tensor, "calc");

              std::string file_str
                = get_trimmed_attribute(xml_tensor, "file");

              if ( !tensor_ok ) {
                LOG("HadronTensorPool", pWARN) << "Problem parsing"
                  << " configuration for the hadron"
                  " tensor for nuclide " << pdg_str << " of type " << type_str
                  << ". It will be ignored";
              }
              else {
                // The tensor is ok, so add its XML attributes to the cache
                std::pair<int, genie::HadronTensorType_t> tensor_id(pdg, type);
                fTensorAttributes[ tensor_id ]
                  = HadronTensorXMLAttributes(calc_str, file_str);
              }
            } // <x> == <tensor>
            xml_tensor = xml_tensor->next;
          } // <tensor> loop
        } // <x> == <nuclide>
        xml_nuclide = xml_nuclide->next;
      } // <nuclide> loop
    } // <x> == <tensor_table>
    xml_tensor_table = xml_tensor_table->next;
  } // <tensor_table> loop

  if ( !found_table ) LOG("HadronTensorPool", pERROR) << "Could not find"
    " a hadron tensor table named \"" << table_to_use << "\" in the XML"
    " configuration file " << filename;

  xmlFreeDoc(xml_doc);

  return kXmlOK;
}
//____________________________________________________________________________
const genie::HadronTensorI* genie::HadronTensorPool::BuildTensor(
  const std::pair<int, genie::HadronTensorType_t>& tensor_id,
  const HadronTensorXMLAttributes& attributes)
{
  // Determine the right way to build the tensor based on its "calc" XML
  // attribute
  if ( attributes.calc == "table" ) {

    bool tensor_ok = true;

    // Tensor values are represented using a 2D grid that is stored in a data
    // file. Get the full path to the file, or an empty string if it could not
    // be found. Also set the tensor_ok flag to false if the file could not be
    // found.
    std::string full_file_name = FindTensorTableFile(attributes.file,
      tensor_ok);

    if ( tensor_ok ) {
      // Create the new hadron tensor object
      LOG("HadronTensorPool", pINFO) << "Loading the hadron"
        << " tensor data file " << full_file_name;
      genie::HadronTensorI* temp_ptr
        = new TabulatedValenciaHadronTensor(full_file_name);

      // Place a pointer to it in the map of loaded tensor objects for easy
      // retrieval.
      /// \todo Switch to using std::unique_ptr here for easy cleanup
      /// once C++11 features are allowed in GENIE
      fTensors[tensor_id] = temp_ptr;

      // Return a pointer to the newly-created hadron tensor object
      return temp_ptr;
    }
    else {
      LOG("HadronTensorPool", pERROR) << "The hadron tensor data file \""
        << attributes.file << "\" requested for target pdg = "
        << tensor_id.first << " and hadron tensor type " << tensor_id.second
        << " could not be found.";
    }

  } // attributes.calc == "table"
  else {
    // Unrecognized "calc" attribute for the currently-requested tensor
    LOG("HadronTensorPool", pERROR) << "Unable to create a hadron tensor"
      << " with calculation method \"" << attributes.calc
      << "\" for target pdg = "
      << tensor_id.first << " and hadron tensor type " << tensor_id.second;
  }

  // If there was a problem, return a null pointer
  return NULL;
}
