//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Steven Gardiner <gardiner \at fnal.gov>
         Liang Liu <liangliu \at fnal.gov>
         Fermi National Accelerator Laboratory

 For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________

// standard library includes
#include <cstdlib>
#include <fstream>
#include <string>

// GENIE includes
#include "Framework/Messenger/Messenger.h"
#include "Physics/HadronTensors/TabulatedHadronTensorModelI.h"
#include "Physics/HadronTensors/TabulatedLabFrameHadronTensor.h"
#include "Physics/HadronTensors/HadronTensorI.h"
#include "Framework/Utils/StringUtils.h"
#include "Framework/Utils/XmlParserUtils.h"

namespace {

  /// Converts a string to a genie::HadronTensorType_t value. If the string
  /// does not correspond to a valid tensor type, kHT_Undefined
  /// is returned, and the ok flag is set to false
  genie::HadronTensorType_t string_to_tensor_type(const std::string& str,
    bool& ok)
  {
    if (str == "MEC_FullAll") return genie::kHT_MEC_FullAll;
    else if (str == "MEC_Fullpn")
      return genie::kHT_MEC_Fullpn;
    else if (str == "MEC_DeltaAll")
      return genie::kHT_MEC_DeltaAll;
    else if (str == "MEC_Deltapn")
      return genie::kHT_MEC_Deltapn;
    else if (str == "MEC_EM")
      return genie::kHT_MEC_EM;
    else if (str == "MEC_EM_pn")
        return genie::kHT_MEC_EM_pn;
    else if (str == "MEC_EM_pp")
        return genie::kHT_MEC_EM_pp;
    else if (str == "MEC_EM_wImag")
      return genie::kHT_MEC_EM_wImag;
    else if (str == "QE_EM")
      return genie::kHT_QE_EM;
    else if (str == "QE_EM_proton")
      return genie::kHT_QE_EM_proton;
    else if (str == "QE_EM_neutron")
      return genie::kHT_QE_EM_neutron;
    else if (str == "MEC_FullAll_Param")
      return genie::kHT_MEC_FullAll_Param;
    else if (str == "MEC_FullAll_wImag")
      return genie::kHT_MEC_FullAll_wImag;
    else if (str == "QE_Full")
      return genie::kHT_QE_Full;

    else if (str == "QE_CRPA_Low")
      return genie::kHT_QE_CRPA_Low;
    else if (str == "QE_CRPA_Medium")
      return genie::kHT_QE_CRPA_Medium;
    else if (str == "QE_CRPA_High")
      return genie::kHT_QE_CRPA_High;

    else if (str == "QE_CRPA_anu_Low")
      return genie::kHT_QE_CRPA_anu_Low;
    else if (str == "QE_CRPA_anu_Medium")
      return genie::kHT_QE_CRPA_anu_Medium;
    else if (str == "QE_CRPA_anu_High")
      return genie::kHT_QE_CRPA_anu_High;

    else if (str == "QE_HF_Low")
      return genie::kHT_QE_HF_Low;
    else if (str == "QE_HF_Medium")
      return genie::kHT_QE_HF_Medium;
    else if (str == "QE_HF_High")
      return genie::kHT_QE_HF_High;

    else if (str == "QE_HF_anu_Low")
      return genie::kHT_QE_HF_anu_Low;
    else if (str == "QE_HF_anu_Medium")
      return genie::kHT_QE_HF_anu_Medium;
    else if (str == "QE_HF_anu_High")
      return genie::kHT_QE_HF_anu_High;


    else if (str == "QE_CRPAPW_Low")
      return genie::kHT_QE_CRPAPW_Low;
    else if (str == "QE_CRPAPW_Medium")
      return genie::kHT_QE_CRPAPW_Medium;
    else if (str == "QE_CRPAPW_High")
      return genie::kHT_QE_CRPAPW_High;

    else if (str == "QE_CRPAPW_anu_Low")
      return genie::kHT_QE_CRPAPW_anu_Low;
    else if (str == "QE_CRPAPW_anu_Medium")
      return genie::kHT_QE_CRPAPW_anu_Medium;
    else if (str == "QE_CRPAPW_anu_High")
      return genie::kHT_QE_CRPAPW_anu_High;

    else if (str == "QE_HFPW_Low")
      return genie::kHT_QE_HFPW_Low;
    else if (str == "QE_HFPW_Medium")
      return genie::kHT_QE_HFPW_Medium;
    else if (str == "QE_HFPW_High")
      return genie::kHT_QE_HFPW_High;

    else if (str == "QE_HFPW_anu_Low")
      return genie::kHT_QE_HFPW_anu_Low;
    else if (str == "QE_HFPW_anu_Medium")
      return genie::kHT_QE_HFPW_anu_Medium;
    else if (str == "QE_HFPW_anu_High")
      return genie::kHT_QE_HFPW_anu_High;

    else if (str == "QE_SuSABlend")
      return genie::kHT_QE_SuSABlend;
    else if (str == "QE_SuSABlend_anu")
      return genie::kHT_QE_SuSABlend_anu;

    else {
      ok = false;
      return genie::kHT_Undefined;
    }
  }

  /// Converts a genie::HadronTensorType_t value to a string
  std::string tensor_type_to_string(genie::HadronTensorType_t htt)
  {
    if ( htt == genie::kHT_MEC_FullAll ) return "MEC_FullAll";
    else if ( htt == genie::kHT_MEC_Fullpn ) return "MEC_Fullpn";
    else if ( htt == genie::kHT_MEC_DeltaAll ) return "MEC_DeltaAll";
    else if ( htt == genie::kHT_MEC_Deltapn ) return "MEC_Deltapn";
    else if ( htt == genie::kHT_MEC_EM ) return "MEC_EM";
    else if ( htt == genie::kHT_MEC_EM_pn ) return "MEC_EM_pn";
    else if ( htt == genie::kHT_MEC_EM_pp ) return "MEC_EM_pp";
    else if ( htt == genie::kHT_MEC_EM_wImag ) return "MEC_EM_wImag";
    else if ( htt == genie::kHT_QE_EM ) return "QE_EM";
    else if ( htt == genie::kHT_QE_EM_proton ) return "QE_EM_proton";
    else if ( htt == genie::kHT_QE_EM_neutron ) return "QE_EM_neutron";
    else if ( htt == genie::kHT_MEC_FullAll_Param ) return "MEC_FullAll_Param";
    else if ( htt == genie::kHT_MEC_FullAll_wImag ) return "MEC_FullAll_wImag";
    else if ( htt == genie::kHT_QE_Full ) return "QE_Full";

    else if ( htt == genie::kHT_QE_CRPA_Low ) return "QE_CRPA_Low";
    else if ( htt == genie::kHT_QE_CRPA_Medium ) return "QE_CRPA_Medium";
    else if ( htt == genie::kHT_QE_CRPA_High ) return "QE_CRPA_High";

    else if ( htt == genie::kHT_QE_CRPA_anu_Low ) return "QE_CRPA_anu_Low";
    else if ( htt == genie::kHT_QE_CRPA_anu_Medium ) return "QE_CRPA_anu_Medium";
    else if ( htt == genie::kHT_QE_CRPA_anu_High ) return "QE_CRPA_anu_High";

    else if ( htt == genie::kHT_QE_HF_Low ) return "QE_HF_Low";
    else if ( htt == genie::kHT_QE_HF_Medium ) return "QE_HF_Medium";
    else if ( htt == genie::kHT_QE_HF_High ) return "QE_HF_High";

    else if ( htt == genie::kHT_QE_HF_anu_Low ) return "QE_HF_anu_Low";
    else if ( htt == genie::kHT_QE_HF_anu_Medium ) return "QE_HF_anu_Medium";
    else if ( htt == genie::kHT_QE_HF_anu_High ) return "QE_HF_anu_High";

    else if ( htt == genie::kHT_QE_CRPAPW_Low ) return "QE_CRPAPW_Low";
    else if ( htt == genie::kHT_QE_CRPAPW_Medium ) return "QE_CRPAPW_Medium";
    else if ( htt == genie::kHT_QE_CRPAPW_High ) return "QE_CRPAPW_High";

    else if ( htt == genie::kHT_QE_CRPAPW_anu_Low ) return "QE_CRPAPW_anu_Low";
    else if ( htt == genie::kHT_QE_CRPAPW_anu_Medium ) return "QE_CRPAPW_anu_Medium";
    else if ( htt == genie::kHT_QE_CRPAPW_anu_High ) return "QE_CRPAPW_anu_High";

    else if ( htt == genie::kHT_QE_HFPW_Low ) return "QE_HFPW_Low";
    else if ( htt == genie::kHT_QE_HFPW_Medium ) return "QE_HFPW_Medium";
    else if ( htt == genie::kHT_QE_HFPW_High ) return "QE_HFPW_High";

    else if ( htt == genie::kHT_QE_HFPW_anu_Low ) return "QE_HFPW_anu_Low";
    else if ( htt == genie::kHT_QE_HFPW_anu_Medium ) return "QE_HFPW_anu_Medium";
    else if ( htt == genie::kHT_QE_HFPW_anu_High ) return "QE_HFPW_anu_High";

    else if ( htt == genie::kHT_QE_SuSABlend ) return "QE_SuSABlend";
    else if ( htt == genie::kHT_QE_SuSABlend_anu ) return "QE_SuSABlend_anu";

    else return "Undefined";
  }

  /// Returns true if a given file exists and is accessible, or false otherwise
  bool file_exists(const std::string& file_name) {
    return std::ifstream(file_name.c_str()).good();
  }

}

//____________________________________________________________________________
genie::TabulatedHadronTensorModelI::TabulatedHadronTensorModelI()
  : genie::HadronTensorModelI()
{

}

//____________________________________________________________________________
genie::TabulatedHadronTensorModelI::TabulatedHadronTensorModelI(std::string name)
  : genie::HadronTensorModelI( name )
{

}

//____________________________________________________________________________
genie::TabulatedHadronTensorModelI::TabulatedHadronTensorModelI(std::string name,
  std::string config) : genie::HadronTensorModelI(name, config)
{

}

//____________________________________________________________________________
void genie::TabulatedHadronTensorModelI::Configure(const Registry& config)
{
  HadronTensorModelI::Configure( config );
  this->LoadConfig();
}
//____________________________________________________________________________
void genie::TabulatedHadronTensorModelI::Configure(std::string config)
{
  HadronTensorModelI::Configure( config );
  this->LoadConfig();
}
//____________________________________________________________________________
void genie::TabulatedHadronTensorModelI::LoadConfig(void)
{
  GetParamDef( "WarnIfMissing", fWarnIfMissing, true );

  // Either a data path relative to the root GENIE folder
  // or an absolute path can be used. Find out which
  // option was chosen.
  std::string path_type;
  GetParamDef( "DataPathType", path_type, std::string("relative") );

  // Right now, there can only be a single data path
  // specified. We use a vector of paths to allow for
  // easy expansion later.
  std::string data_path;
  GetParam( "DataPath", data_path );

  // Convert the relative path to an absolute one if needed
  if ( path_type == "relative" ) {
    data_path = std::string( gSystem->Getenv("GENIE") ) + '/' + data_path;
  }

  fDataPaths.push_back( data_path );
}

//____________________________________________________________________________
genie::TabulatedHadronTensorModelI::~TabulatedHadronTensorModelI()
{
  std::map< HadronTensorID, HadronTensorI* >::iterator it;
  for (it = fTensors.begin(); it != fTensors.end(); ++it) {
    HadronTensorI* t = it->second;
    if ( t ) delete t;
  }
  fTensors.clear();
}
//____________________________________________________________________________
const genie::HadronTensorI* genie::TabulatedHadronTensorModelI::GetTensor(
  int tensor_pdg, genie::HadronTensorType_t type) const
{
  HadronTensorID temp_id(tensor_pdg, type);

  // First check to see if the hadron tensor object already exists in memory.
  // If it does, return the existing object.
  if ( fTensors.count(temp_id) ) return fTensors.find(temp_id)->second;

  // If not, try to create it
  const HadronTensorI* ht = this->BuildTensor( temp_id );

  if ( !ht && fWarnIfMissing ) {
    LOG("TabulatedHadronTensorModelI", pWARN) << "Unable to create a hadron tensor"
      << " for target pdg = " << temp_id.target_pdg
      << " and hadron tensor type " << temp_id.type;
  }

  return ht;
}

//____________________________________________________________________________
std::string genie::TabulatedHadronTensorModelI::FindTensorTableFile(
  const std::string& basename, bool& ok) const
{
  for (size_t p = 0; p < fDataPaths.size(); ++p) {
    const std::string& path = fDataPaths.at( p );
    std::string full_name = path + '/' + basename;
    if ( file_exists(full_name) ) return full_name;
  }

  // A matching file could not be found
  ok = false;
  return std::string();
}

//____________________________________________________________________________
const genie::HadronTensorI* genie::TabulatedHadronTensorModelI::BuildTensor(
  const HadronTensorID& tensor_id) const
{
  bool tensor_ok = true;

  std::string tensor_file_basename = this->GetTensorFileBasename( tensor_id );

  // Tensor values are represented using a 2D grid that is stored in a data
  // file. Get the full path to the file, or an empty string if it could not
  // be found. Also set the tensor_ok flag to false if the file could not be
  // found.
  std::string full_file_name = FindTensorTableFile(tensor_file_basename,
    tensor_ok);

  if ( tensor_ok ) {

    // Create the new hadron tensor object
    LOG("TabulatedHadronTensorModelI", pINFO) << "Loading the hadron"
      << " tensor data file " << full_file_name;

    genie::HadronTensorI* temp_ptr = this->ParseTensorFile( full_file_name );

    // Place a pointer to it in the map of loaded tensor objects for easy
    // retrieval.
    /// \todo Switch to using std::unique_ptr here for easy cleanup
    /// once C++11 features are allowed in GENIE
    fTensors[tensor_id] = temp_ptr;

    // Return a pointer to the newly-created hadron tensor object
    return temp_ptr;
  }

  else {
    // If we couldn't make the hadron tensor, store a nullptr to avoid
    // unsuccessful repeat attempts. These can otherwise slow things down
    // for no good reason.
    fTensors[tensor_id] = NULL;

    if ( fWarnIfMissing ) {
      LOG("TabulatedHadronTensorModelI", pERROR) << "The hadron tensor data file \""
        << full_file_name << "\" requested for target pdg = "
        << tensor_id.target_pdg << " and hadron tensor type "
        << tensor_id.type << " could not be found.";
    }
  }

  // If there was a problem, return a null pointer
  return NULL;
}

//____________________________________________________________________________
std::string genie::TabulatedHadronTensorModelI::GetTensorFileBasename(
  const HadronTensorID& ht_id ) const
{

  std::string tgt_string;
  std::stringstream ss;
  ss <<  ht_id.target_pdg;
  tgt_string = ss.str();

  RgKey key = tensor_type_to_string( ht_id.type ) + "@Pdg="
    + tgt_string;

  std::string basename;
  GetParamDef( key, basename, std::string("TENSOR_FILE_NOT_FOUND") );

  return basename;
}
