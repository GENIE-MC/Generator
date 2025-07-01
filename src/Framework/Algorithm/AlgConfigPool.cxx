//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Costas Andreopoulos <c.andreopoulos \at cern.ch>
 University of Liverpool
*/
//____________________________________________________________________________

#include <cstring>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <fstream>

#include "libxml/xmlmemory.h"
#include "libxml/parser.h"

#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Registry/RegistryItemTypeDef.h"
#include "Framework/Utils/XmlParserUtils.h"

#include "Framework/Utils/StringUtils.h"

using std::setw;
using std::setfill;
using std::endl;
using std::ostringstream;

using namespace genie;

//____________________________________________________________________________
namespace genie {
  ostream & operator<<(ostream & stream, const AlgConfigPool & config_pool)
  {
    config_pool.Print(stream);
    return stream;
  }
}
//____________________________________________________________________________
AlgConfigPool * AlgConfigPool::fInstance = 0;
//____________________________________________________________________________
AlgConfigPool::AlgConfigPool()
{
  if( ! this->LoadAlgConfig() )
  LOG("AlgConfigPool", pERROR) << "Could not load XML config file";
  fInstance =  0;
}
//____________________________________________________________________________
AlgConfigPool::~AlgConfigPool()
{
// Clean up and report the most important physics params used in this instance.
// Don't clutter output if exiting in err.

  if(!gAbortingInErr) {
/*
    cout << "AlgConfigPool singleton dtor: "
         << "Deleting all owned algorithm configurations" << endl;
*/
  }
  map<string, Registry *>::iterator citer;
  for(citer = fRegistryPool.begin(); citer != fRegistryPool.end(); ++citer) {
    string key = citer->first;
    Registry * config = citer->second;
    if(config) {
      delete config;
      config = 0;
    }
  }
  fRegistryPool.clear();
  fConfigFiles.clear();
  fConfigKeyList.clear();
  fInstance = 0;
}
//____________________________________________________________________________
AlgConfigPool * AlgConfigPool::Instance()
{
  if(fInstance == 0) {
    static AlgConfigPool::Cleaner cleaner;
    cleaner.DummyMethodAndSilentCompiler();
    fInstance = new AlgConfigPool;
  }
  return fInstance;
}
//____________________________________________________________________________
bool AlgConfigPool::LoadAlgConfig(void)
{
// Loads all algorithm XML configurations and creates a map with all loaded
// configuration registries

  SLOG("AlgConfigPool", pINFO)
        << "AlgConfigPool late initialization: Loading all XML config. files";

  //-- read the global parameter lists
  if(!this->LoadGlobalParamLists()) 
  {
    SLOG("AlgConfigPool", pERROR)
      << "Global parameter lists not available "
      "This can be normal if you are running "
      "Comparisons/Reweight in some cases";
  }

  //-- read the MASTER_CONFIG XML file
  if(!this->LoadMasterConfigs()) 
  {
    SLOG("AlgConfigPool", pERROR)
      << "Master config file not available";
  }

  //-- read Tune Generator List for the tune, if available
  if( ! LoadTuneGeneratorList() ) {
    SLOG( "AlgConfigPool", pWARN ) << "Tune generator List not available" ;
  }

  //-- loop over all XML config files and read all named configuration
  //   sets for each algorithm
  map<string, string>::const_iterator conf_file_iter;

  for(conf_file_iter = fConfigFiles.begin();
                  conf_file_iter != fConfigFiles.end(); ++conf_file_iter) {

    string alg_name    = conf_file_iter->first;
    string file_name   = conf_file_iter->second;

    SLOG("AlgConfigPool", pINFO)
         << setfill('.') << setw(40) << alg_name << " -> " << file_name;

    string full_path = utils::xml::GetXMLFilePath(file_name);
    SLOG("AlgConfigPool", pNOTICE)
      << "*** GENIE XML config file " << full_path;
    bool ok = this->LoadSingleAlgConfig(alg_name, full_path);
    if(!ok) {
      SLOG("AlgConfigPool", pERROR)
           << "Error in loading config sets for algorithm = " << alg_name;
    }
  }
  return true;
};
//____________________________________________________________________________
bool AlgConfigPool::LoadMasterConfig(std::string configname)
{
// Loads the master config XML file: the file that specifies which XML config
// file to load for each algorithm

  //-- get the master config XML file using GXMLPATH + default locations
  // fMasterConfig = utils::xml::GetXMLFilePath("master_config.xml");
  fMasterConfig = utils::xml::GetXMLFilePath(configname);

  bool is_accessible = ! (gSystem->AccessPathName( fMasterConfig.c_str() ));
  if (!is_accessible) {
     SLOG("AlgConfigPool", pERROR)
       << "The XML doc doesn't exist! (filename : " << fMasterConfig << ")";
     return false;
  }

  xmlDocPtr xml_doc = xmlParseFile(fMasterConfig.c_str());
  if(xml_doc==NULL) {
     SLOG("AlgConfigPool", pERROR)
      << "The XML doc can't be parsed! (filename : " << fMasterConfig << ")";
     return false;
  }

  xmlNodePtr xml_root = xmlDocGetRootElement(xml_doc);
  if(xml_root==NULL) {
     SLOG("AlgConfigPool", pERROR)
             << "The XML doc is empty! (filename : " << fMasterConfig << ")";
     xmlFreeDoc(xml_doc);
     return false;
  }

  if( xmlStrcmp(xml_root->name, (const xmlChar *) "genie_config") ) {
     SLOG("AlgConfigPool", pERROR)
              << "The XML doc has invalid root element! "
                                   << "(filename : " << fMasterConfig << ")";
     xmlFreeDoc(xml_doc);
     return false;
  }

  // loop over all xml tree nodes (<alg_config>) that are children of the
  // root node and read the config file name for each registered algorithm
  xmlNodePtr xml_ac = xml_root->xmlChildrenNode;
  while (xml_ac != NULL) {
    if( (!xmlStrcmp(xml_ac->name, (const xmlChar *) "config")) ) {

       string alg_name  = utils::str::TrimSpaces(
                     utils::xml::GetAttribute(xml_ac, "alg"));
       string config_file = utils::xml::TrimSpacesClean(
                xmlNodeListGetString(xml_doc, xml_ac->xmlChildrenNode, 1));

       pair<string, string> alg_conf(alg_name, config_file);
       fConfigFiles.insert(alg_conf);
    }
    xml_ac = xml_ac->next;
  }
  xmlFreeNode(xml_ac);
  xmlFreeDoc(xml_doc);
  return true;
}

bool AlgConfigPool::LoadMasterConfigs(void) {
  auto main = LoadMasterConfig("master_config.xml");
  if (std::getenv("GENIE_REWEIGHT")) {
    auto rew_main = LoadMasterConfig("reweight_master_config.xml");
    return main && rew_main;
  } else
    return main;
}
//____________________________________________________________________________
bool AlgConfigPool::LoadGlobalParamLists(void)
{
// Load the global parameter list (a list of physics constants at a given MC
// job, that is allowed to be modified to fine tune the generator output)
//
  SLOG("AlgConfigPool", pINFO) << "Loading global parameter lists";

  // -- get the user config XML file using GXMLPATH + default locations
  string glob_params = utils::xml::GetXMLFilePath("ModelConfiguration.xml");

  // fixed key prefix
  string key_prefix = "GlobalParameterList";

  // load and report status
  return this->LoadRegistries(key_prefix, glob_params, "global_param_list");
}
//____________________________________________________________________________
bool AlgConfigPool::LoadCommonLists( const string & file_id )
{
// Load the common parameter list
//
  SLOG("AlgConfigPool", pINFO) << "Loading Common " << file_id << " lists";

  // -- get the user config XML file using GXMLPATH + default locations
  std::string xml_name = "Common" + file_id + ".xml" ;
  string full_path = utils::xml::GetXMLFilePath( xml_name );

  // fixed key prefix
  string key_prefix = "Common" + file_id + "List";

  // load and report status
  if ( ! this->LoadRegistries(key_prefix, full_path, "common_"+file_id+"_list") ) {

	  SLOG("AlgConfigPool", pERROR) << "Failed to load Common " << file_id ;
	  return false ;
  }

  return true ;
}
//____________________________________________________________________________
bool AlgConfigPool::LoadTuneGeneratorList(void)
{
// Load the common parameter list
//
  SLOG("AlgConfigPool", pINFO) << "Loading Tune Gerator List";

  // -- get the user config XML file using GXMLPATH + default locations
  string generator_list_file = utils::xml::GetXMLFilePath("TuneGeneratorList.xml");

  // fixed key prefix
  string key_prefix = "TuneGeneratorList";

  // load and report status
  return this->LoadRegistries(key_prefix, generator_list_file, "tune_generator_list");
}
//____________________________________________________________________________

bool AlgConfigPool::LoadSingleAlgConfig(string alg_name, string file_name)
{
// Loads all configuration sets for the input algorithm that can be found in
// the input XML file

  // use the algorithm name as the key prefix
  string key_prefix = alg_name;

  // load and report status
  return this->LoadRegistries(key_prefix, file_name, "alg_conf");
}
//____________________________________________________________________________
bool AlgConfigPool::LoadRegistries(
                             string key_prefix, string file_name, string root)
{
// Loads all the configuration registries from the input XML file

  SLOG("AlgConfigPool", pDEBUG) << "[-] Loading registries:";

  bool is_accessible = ! (gSystem->AccessPathName(file_name.c_str()));
  if (!is_accessible) {
     SLOG("AlgConfigPool", pERROR)
       << "The XML doc doesn't exist! (filename : " << file_name << ")";
     return false;
  }

  xmlDocPtr xml_doc = xmlParseFile( file_name.c_str() );
  if(xml_doc==NULL) {
     SLOG("AlgConfigPool", pERROR)
      << "The XML document can't be parsed! (filename : " << file_name << ")";
     return false;
  }

  xmlNodePtr xml_cur = xmlDocGetRootElement( xml_doc );
  if(xml_cur==NULL) {
     SLOG("AlgConfigPool", pERROR)
             << "The XML document is empty! (filename : " << file_name << ")";
     xmlFreeDoc(xml_doc);
     return false;
  }
  if( xmlStrcmp(xml_cur->name, (const xmlChar *) root.c_str()) ) {
     SLOG("AlgConfigPool", pERROR)
              << "The XML document has invalid root element! "
                                        << "(filename : " << file_name << ")";
     xmlFreeNode(xml_cur);
     xmlFreeDoc(xml_doc);
     return false;
  }

  // loop over all xml tree nodes that are children of the root node
  xml_cur = xml_cur->xmlChildrenNode;
  while (xml_cur != NULL) {
    // enter everytime you find an 'param_set' tag
    if( (!xmlStrcmp(xml_cur->name, (const xmlChar *) "param_set")) ) {

      string param_set  = utils::str::TrimSpaces(
                           utils::xml::GetAttribute(xml_cur, "name"));

      // build the registry key
      ostringstream key;
      key << key_prefix << "/" << param_set;

      // store the key in the key list
      fConfigKeyList.push_back(key.str());

      // create a new Registry and fill it with the configuration params
      Registry * config = new Registry(param_set,false);

      xmlNodePtr xml_param = xml_cur->xmlChildrenNode;
      while (xml_param != NULL) {
        if( (!xmlStrcmp(xml_param->name, (const xmlChar *) "param")) ) {

            string param_type =
                   utils::str::TrimSpaces(
                       utils::xml::GetAttribute(xml_param, "type"));
            string param_name =
                   utils::str::TrimSpaces(
                       utils::xml::GetAttribute(xml_param, "name"));
            string param_value =
                    utils::xml::TrimSpacesClean(
                               xmlNodeListGetString(
                                 xml_doc, xml_param->xmlChildrenNode, 1));


	    if ( param_type.find( "vec-" ) == 0 ) {

	      param_type = param_type.substr( 4 ) ;

	      string delim = utils::str::TrimSpaces( utils::xml::GetAttribute(xml_param, "delim"));

	      this -> AddParameterVector( config, param_type, param_name, param_value, delim ) ;
	    }
            else if (param_type.find( "mat-" ) == 0) {
              param_type = param_type.substr( 4 ) ;
              LOG("AlgConfigPool", pNOTICE) << "Liang Liu" << param_type ;
                string importfile = utils::str::TrimSpaces( utils::xml::GetAttribute(xml_param, "importfile"));
                if(!importfile.compare("false")){
                  string rowdelim = utils::str::TrimSpaces( utils::xml::GetAttribute(xml_param, "rowdelim"));
                  string coldelim = utils::str::TrimSpaces( utils::xml::GetAttribute(xml_param, "coldelim"));
                  this -> AddParameterMatrix(config, param_type, param_name, param_value, rowdelim, coldelim);
                }
                else if(!importfile.compare("true")){
                  this -> AddParameterMatrix(config, param_type, param_name, param_value );
                }
            }
	    else this->AddConfigParameter( config,
					   param_type, param_name,
					   param_value);
        }
        xml_param = xml_param->next;
      }
      //xmlFree(xml_param);
      xmlFreeNode(xml_param);
      config->SetName(param_set);
      config->Lock();

      pair<string, Registry *> single_reg(key.str(), config);
      fRegistryPool.insert(single_reg);

      SLOG("AlgConfigPool", pDEBUG) << " |---o " << key.str();
    }
    xml_cur = xml_cur->next;
  }
  //xmlFree(xml_cur);
  xmlFreeNode(xml_cur);
  //xmlFree(xml_doc);
  xmlFreeDoc(xml_doc);

  return true;
}
//____________________________________________________________________________
int  AlgConfigPool::AddParameterVector  (Registry * r, string pt, string pn, string pv,
					 const string & delim ) {

  // Adds a configuration parameter vector
  // It is simply add a number of entries in the Registy
  // The name scheme starts from the name and it goes like
  // 'N'+pn+'s' that will be an integer with the number of entries.
  // Each entry will be named pn+"-i" where i is replaced by the number

  SLOG("AlgConfigPool", pDEBUG)
    << "Adding Parameter Vector [" << pt << "]: Key = "
    << pn << " -> Value = " << pv;

  vector<string> bits = utils::str::Split( pv, delim ) ;

  string n_name = Algorithm::BuildParamVectSizeKey( pn ) ;

  std::stringstream n_value ;
  n_value << bits.size() ;

  this->AddConfigParameter(r, "int", n_name, n_value.str() );

  for ( unsigned int i = 0 ; i < bits.size() ; ++i ) {

    std::string name = Algorithm::BuildParamVectKey( pn, i ) ;

    this -> AddConfigParameter( r, pt, name, utils::str::TrimSpaces( bits[i] ) );

  }

  return bits.size() ;

}

//____________________________________________________________________________
int  AlgConfigPool::AddParameterMatrix  (Registry * r, string pt, string pn, string pv,
					 const string & rowdelim, const string & coldelim ) {

  // Adds a configuration parameter matrix
  // It is simply add a number of entries in the Registy
  // The name scheme starts from the name and it goes like
  // 'Nrow'+pn+'s' for the size of rows and 
  // 'Ncol'+pn+'s' for the size of columns
  // that will be an integer with the number of entries.
  // Each entry will be named pn+"-i"+"-j" where i and j are 
  // index of row and column and replaced by the number

  if(!rowdelim.compare(coldelim) || rowdelim.empty() || coldelim.empty()) {
    LOG("AlgConfigPool", pFATAL) << "row and column have wrong delims: " << rowdelim << "  " << coldelim ;
    exit(1);
  }
  SLOG("AlgConfigPool", pDEBUG)
    << "Adding Parameter Matrix [" << pt << "]: Key = "
    << pn << " -> Value = " << pv;

  vector<string> mat_row = utils::str::Split( pv, rowdelim ) ;

  string r_name = Algorithm::BuildParamMatRowSizeKey( pn ) ;


  unsigned int n_row = 0, n_col = 0;
  std::stringstream r_value ;
  r_value << mat_row.size() ;
  n_row = mat_row.size() ;

  this->AddConfigParameter(r, "int", r_name, r_value.str() );

  for ( unsigned int i = 0 ; i < mat_row.size() ; ++i ) {
    vector<string> bits = utils::str::Split( mat_row[i], coldelim ) ;
    if(i == 0){
      string c_name = Algorithm::BuildParamMatColSizeKey( pn ) ;
      std::stringstream c_value ;
      c_value << bits.size();
      n_col = bits.size();
      this->AddConfigParameter(r, "int", c_name, c_value.str() );
    }
    else{
      if(n_col != bits.size()){
        LOG("AlgConfigPool", pFATAL) << "wrong size of matrix in row: " << i;
        exit(1);
      }
    }

    for ( unsigned int j = 0 ; j < bits.size() ; ++j ) {

      std::string name = Algorithm::BuildParamMatKey( pn, i, j ) ;

      this -> AddConfigParameter( r, pt, name, utils::str::TrimSpaces( bits[j] ) );

    }
  }
  return n_row * n_col;

}
//____________________________________________________________________________
int  AlgConfigPool::AddParameterMatrix  (Registry * r, string pt, string pn, string pv ) {

  // Adds a configuration parameter matrix
  // It is simply add a number of entries in the Registy
  // The name scheme starts from the name and it goes like
  // 'Nrow'+pn+'s' for the size of rows and 
  // 'Ncol'+pn+'s' for the size of columns
  // that will be an integer with the number of entries.
  // Each entry will be named pn+"-i"+"-j" where i and j are 
  // index of row and column and replaced by the number

  char * GENIE_PATH = std::getenv("GENIE");
  string filepath = std::string(GENIE_PATH) + "/" + pv;
  std::ifstream file(filepath);
  if (!file.is_open()) {
    throw std::runtime_error("Could not open file");
  }

  std::string line;
  int n_row = 0, n_col = 0;
  int i_row = 0;
  while (getline(file, line)) {
    size_t start = line.find_first_not_of(" \t");
    std::string trimmedLine = (start == std::string::npos) ? "" : line.substr(start);
    if(trimmedLine.empty()) continue;
    if(!trimmedLine.empty() && trimmedLine[0] == '%') continue;

    std::istringstream iss(line);
    double value;
    int i_col = 0;
    while (iss >> value) {
      std::string name = Algorithm::BuildParamMatKey( pn, i_row, i_col ) ;
      this -> AddConfigParameter( r, pt, name, utils::str::TrimSpaces( std::to_string(value) ) );
      i_col++;
    }
    if(i_row == 0)
      n_col = i_col;
    else{
      if(n_col != i_col){
         LOG("AlgConfigPool", pFATAL) << "wrong size of matrix in row: " << i_row;
         exit(1);
      }
    }
    i_row++;
  }
  n_row = i_row;
  std::stringstream r_value ;
  r_value << n_row ;
  string r_name = Algorithm::BuildParamMatRowSizeKey( pn ) ;
  this->AddConfigParameter(r, "int", r_name, r_value.str() );
  std::stringstream c_value ;
  c_value << n_col ;
  string c_name = Algorithm::BuildParamMatColSizeKey( pn ) ;
  this->AddConfigParameter(r, "int", c_name, c_value.str() );
  return n_row*n_col;
}



//____________________________________________________________________________
void AlgConfigPool::AddConfigParameter( Registry * r,
					string ptype, string pname, string pvalue)
{
// Adds a configuration parameter with type = ptype, key = pname and value =
// pvalue at the input configuration registry r

  SLOG("AlgConfigPool", pDEBUG)
    << "Adding Parameter [" << ptype << "]: Key = "
    << pname << " -> Value = " << pvalue;

  bool isRootObjParam = (strcmp(ptype.c_str(), "h1f")    == 0) ||
                        (strcmp(ptype.c_str(), "Th2f")   == 0) ||
                        (strcmp(ptype.c_str(), "tree")   == 0);
  bool isBasicParam   = (strcmp(ptype.c_str(), "int")    == 0) ||
                        (strcmp(ptype.c_str(), "bool")   == 0) ||
                        (strcmp(ptype.c_str(), "double") == 0) ||
                        (strcmp(ptype.c_str(), "string") == 0) ||
                        (strcmp(ptype.c_str(), "alg")    == 0);


  if     (isBasicParam)   this->AddBasicParameter  (r, ptype, pname, pvalue);
  else if(isRootObjParam) this->AddRootObjParameter(r, ptype, pname, pvalue);
  else {
    SLOG("AlgConfigPool", pERROR)
      << "Parameter [" << ptype << "]: Key = " << pname
      << " -> Value = " << pvalue << " could not be added";
  }

}
//____________________________________________________________________________
void AlgConfigPool::AddBasicParameter(
                     Registry * r, string ptype, string pname, string pvalue)
{
  RgKey key = pname;

  if (ptype=="double") {
    RgDbl item = (double) atof(pvalue.c_str());
    r->Set(key, item);
  }
  else if (ptype=="int") {
    RgInt item = (int) atoi(pvalue.c_str());
    r->Set(key, item);
  }
  else if (ptype=="bool") {
    if      (pvalue=="true" ) r->Set(key, true );
    else if (pvalue=="TRUE" ) r->Set(key, true );
    else if (pvalue=="1"    ) r->Set(key, true );
    else if (pvalue=="false") r->Set(key, false);
    else if (pvalue=="FALSE") r->Set(key, false);
    else if (pvalue=="0"    ) r->Set(key, false);
    else {
      LOG("AlgConfigPool", pERROR)
        << "Could not set bool param: " << key;
    }
  }
  else if (ptype=="string") {
    RgStr item = pvalue;
    r->Set(key, item);
  }
  else if (ptype=="alg") {
    string name, config;
    vector<string> algv = utils::str::Split(pvalue, "/");
    if (algv.size()==2) {
      name   = algv[0];
      config = algv[1];
    }
    else if (algv.size()==1) {
      name   = algv[0];
      config = "Default";
    } else {
       LOG("AlgConfigPool", pFATAL)
             << "Unrecognized algorithm id: " << pvalue;
       exit(1);
    }
    RgAlg item(name,config);
    r->Set(key, item);
  }
  else {
   LOG("AlgConfigPool", pERROR)
        << "Config. parameter: " << key
                 << "has unrecognized type: " << ptype;
  }
}
//____________________________________________________________________________
void AlgConfigPool::AddRootObjParameter(
                     Registry * r, string ptype, string pname, string pvalue)
{
  // the ROOT object is given in the XML config file as
  // <param> object_name@root_file_name </param>
  vector<string> rootobjv = utils::str::Split(pvalue, "@");

  if(rootobjv.size() != 2) {
    SLOG("AlgConfigPool", pWARN)
        << "ROOT objects are added in XML config files as: "
             << "object-name@file-name. Wrong syntax in: [" << pvalue << "]";
    SLOG("AlgConfigPool", pERROR)
            << "Parameter [" << ptype << "]: Key = " << pname
                        << " -> Value = " << pvalue << " could not be added";
  }

  string rootobj  = rootobjv[0];
  string rootfile = rootobjv[1];

  TFile f(rootfile.c_str(), "read");

  if (ptype=="h1f") {
    TH1F * h  = (TH1F*) f.Get(rootobj.c_str());
    if(h) {
      TH1F * ch = new TH1F(*h); // clone
		ch->SetDirectory(0);
      r->Set(pname,ch);
    } else {
      SLOG("AlgConfigPool", pERROR)
         << "No TH1F named = " << rootobj << " in ROOT file = " << rootfile;
    }
  } else if (ptype=="h2f") {
    TH2F * h2  = (TH2F*) f.Get(rootobj.c_str());
    if(h2) {
      TH2F * ch2 = new TH2F(*h2); // clone
		ch2->SetDirectory(0);
      r->Set(pname,ch2);
    } else {
      SLOG("AlgConfigPool", pERROR)
         << "No TH2F named = " << rootobj << " in ROOT file = " << rootfile;
    }
  } else if (ptype=="tree") {
    TTree * t  = (TTree*) f.Get(rootobj.c_str());
    if(t) {
      //TTree * ct = new TTree(*t); // clone
      TTree * ct = t->CopyTree("1");
      r->Set(pname,ct);
    } else {
      SLOG("AlgConfigPool", pERROR)
        << "No TTree named = " << rootobj << " in ROOT file = " << rootfile;
    }
  }
  else {}
}
//____________________________________________________________________________
Registry * AlgConfigPool::FindRegistry(const Algorithm * algorithm) const
{
  string key = algorithm->Id().Key();
  return this->FindRegistry(key);
}
//____________________________________________________________________________
Registry * AlgConfigPool::FindRegistry(const AlgId & algid) const
{
  string key = algid.Key();
  return this->FindRegistry(key);
}
//____________________________________________________________________________
Registry* AlgConfigPool::FindRegistry(string alg_name, string param_set) const
{
  AlgId id(alg_name,param_set);
  string key = id.Key();
  return this->FindRegistry(key);
}
//____________________________________________________________________________
Registry * AlgConfigPool::FindRegistry(string key) const
{
  LOG("AlgConfigPool", pDEBUG) << "Searching for registry with key " << key;

  if( fRegistryPool.count(key) == 1 ) {
     map<string, Registry *>::const_iterator config_entry =
                                                   fRegistryPool.find(key);
     return config_entry->second;
  } else {
     LOG("AlgConfigPool", pDEBUG) << "No config registry for key " << key;
     return 0;
  }
  return 0;
}
//____________________________________________________________________________
Registry * AlgConfigPool::GlobalParameterList(void) const
{
  string glob_param_set = (gSystem->Getenv("GUSERPHYSOPT")) ?
                        string(gSystem->Getenv("GUSERPHYSOPT")) : "Default";
  ostringstream key;
  key << "GlobalParameterList/" << glob_param_set;

  return this->FindRegistry(key.str());
}
//____________________________________________________________________________
Registry * AlgConfigPool::CommonList( const string & file_id, const string & set_name ) const
{

  ostringstream key;
  key << "Common" << file_id << "List/" << set_name;

  if ( ! this->FindRegistry(key.str()) ) {
	const_cast<AlgConfigPool*>( this ) -> LoadCommonLists( file_id ) ;
  }

  return this->FindRegistry(key.str()) ;
}
//____________________________________________________________________________
Registry * AlgConfigPool::TuneGeneratorList( void ) const
{

  ostringstream key;
  key << "TuneGeneratorList/Default";

  return this->FindRegistry(key.str());
}
//____________________________________________________________________________
const vector<string> & AlgConfigPool::ConfigKeyList(void) const
{
  return fConfigKeyList;
}
//____________________________________________________________________________
void AlgConfigPool::Print(ostream & stream) const
{
  string frame(100,'~');

  typedef map<string, Registry *>::const_iterator  sregIter;
  typedef map<string, Registry *>::size_type       sregSize;

  sregSize size = fRegistryPool.size();

  stream << frame
         << endl << "Algorithm Configuration Pool ("
         << size << " configuration sets found)"
         << endl << frame << endl;

  sregIter iter = fRegistryPool.begin();
  for( ; iter != fRegistryPool.end(); iter++) {
     const Registry & reg = *(iter->second);
     stream << iter->first;
     stream << reg;
     stream << frame << endl;
  }
}
//____________________________________________________________________________
