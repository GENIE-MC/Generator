//____________________________________________________________________________
/*!

\class    genie::AlgConfigPool

\brief    A singleton class which holds all configuration registries assembled
          from XML configuration files for all (agorithm-name, parameter-set)
          pairs. Any algorithmic object can get an instance of the algorithm
          config. pool and query it to learn its configuration parameters.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 06, 2004

*/
//____________________________________________________________________________

#include <sstream>

#include "libxml/xmlmemory.h"
#include "libxml/parser.h"

#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>

#include "Algorithm/AlgConfigPool.h"
#include "Messenger/Messenger.h"
#include "Utils/XmlParserUtils.h"

using namespace genie;

using std::endl;
using std::ostringstream;

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

  //-- get base GENIE directory from $GENIE environmental variable & build
  //   the GENIE config dir name
  string base_dir   = string( gSystem->Getenv("GENIE") );
  string config_dir = base_dir + string("/config");

  //-- read the MASTER_CONFIG XML file
  fMasterConfig = config_dir + string("/master_conf.xml");
  if(!this->LoadMasterConfig()) return false;

  //-- loop over all XML config files and read all named configuration
  //   sets for each algorithm
  map<string, string>::const_iterator conf_file_iter;

  for(conf_file_iter = fConfigFiles.begin();
                  conf_file_iter != fConfigFiles.end(); ++conf_file_iter) {

    string alg_name    = conf_file_iter->first;
    string config_file = config_dir + "/" + conf_file_iter->second;

    SLOG("AlgConfigPool", pINFO) << alg_name << " ---> " << config_file;
    this->LoadSingleAlgConfig(alg_name, config_file);
  }
  return true;
};
//____________________________________________________________________________
bool AlgConfigPool::LoadMasterConfig(void)
{
// Loads the master config XML file: the file that specifies which XML config
// file to load for each algorithm

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
     return false;
  }

  if( xmlStrcmp(xml_root->name, (const xmlChar *) "genie_config") ) {
     SLOG("AlgConfigPool", pERROR)
              << "The XML doc has invalid root element! "
                                   << "(filename : " << fMasterConfig << ")";
     return false;
  }

  // loop over all xml tree nodes (<alg_config>) that are children of the
  // root node and read the config file name for each registered algorithm

  xmlNodePtr xml_ac = xml_root->xmlChildrenNode;
  while (xml_ac != NULL) {
    if( (!xmlStrcmp(xml_ac->name, (const xmlChar *) "alg_config")) ) {

       string alg_name  = string_utils::TrimSpaces(
                    XmlParserUtils::GetAttribute(xml_ac, "alg_name"));

       // loop over all xml tree nodes that are children of the <alg_config>
       // node, find the <config_file> node and read the filename.
       string alg_filename="";
       xmlNodePtr xml_cf = xml_ac->xmlChildrenNode;
       while (xml_cf != NULL) {
          if( (!xmlStrcmp(xml_cf->name, (const xmlChar *) "config_file")) ) {
             alg_filename = XmlParserUtils::TrimSpaces(
                xmlNodeListGetString(xml_doc, xml_cf->xmlChildrenNode, 1));
          }
          xml_cf = xml_cf->next;
       }
       xmlFree(xml_cf);
       pair<string, string> alg_conf(alg_name, alg_filename);
       fConfigFiles.insert(alg_conf);
    }
    xml_ac = xml_ac->next;
  }
  xmlFree(xml_ac);
  xmlFree(xml_root);
  return true;
}
//____________________________________________________________________________
bool AlgConfigPool::LoadSingleAlgConfig(string alg_name, string file_name)
{
// Loads all configuration sets for the input algorithm using the input XML
// file
  ostringstream allconf;
  allconf << "loading configs:";

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
     return false;
  }
  if( xmlStrcmp(xml_cur->name, (const xmlChar *) "alg_conf") ) {
     SLOG("AlgConfigPool", pERROR)
              << "The XML document has invalid root element! "
                                        << "(filename : " << file_name << ")";
     return false;
  }

  // loop over all xml tree nodes that are children of the root node

  xml_cur = xml_cur->xmlChildrenNode;
  while (xml_cur != NULL) {
    // enter everytime you find an 'param_set' tag
    if( (!xmlStrcmp(xml_cur->name, (const xmlChar *) "param_set")) ) {

      string param_set  = string_utils::TrimSpaces(
                           XmlParserUtils::GetAttribute(xml_cur, "name"));
      string config_key = alg_name + "/" + param_set;

      // create a new Registry and fill it with the configuration params
      Registry * config = new Registry();

      xmlNodePtr xml_param = xml_cur->xmlChildrenNode;
      while (xml_param != NULL) {
        if( (!xmlStrcmp(xml_param->name, (const xmlChar *) "param")) ) {

            string param_type =
                   string_utils::TrimSpaces(
                         XmlParserUtils::GetAttribute(xml_param, "type"));
            string param_name =
                   string_utils::TrimSpaces(
                         XmlParserUtils::GetAttribute(xml_param, "name"));
            string param_value =
                   XmlParserUtils::TrimSpaces(
                               xmlNodeListGetString(
                                 xml_doc, xml_param->xmlChildrenNode, 1));
            this->AddConfigParameter(
                             config, param_type, param_name, param_value);
        }
        xml_param = xml_param->next;
      }
      xmlFree(xml_param);

      config->SetName(param_set);
      config->Lock();

      pair<string, Registry *> alg_single_conf(config_key, config);
      fRegistryPool.insert(alg_single_conf);

      allconf << " [" << param_set << "]";
    }
    xml_cur = xml_cur->next;
  }
  xmlFree(xml_cur);
  SLOG("AlgConfigPool", pINFO) << allconf.str();
  return true;
}
//____________________________________________________________________________
void AlgConfigPool::AddConfigParameter(
                      Registry * r, string ptype, string pname, string pvalue)
{
// Adds a configuration parameter with type = ptype, key = pname and value =
// pvalue at the input configuration registry r

  SLOG("AlgConfigPool", pDEBUG)
          << "Adding Parameter [" << ptype << "]: Key = "
                                         << pname << " -> Value = " << pvalue;

  bool isRootObjParam = (strcmp(ptype.c_str(), "TH1F")   == 0) ||
                        (strcmp(ptype.c_str(), "TH2F")   == 0) ||
                        (strcmp(ptype.c_str(), "TTree")  == 0);
  bool isBasicParam   = (strcmp(ptype.c_str(), "int")    == 0) ||
                        (strcmp(ptype.c_str(), "bool")   == 0) ||
                        (strcmp(ptype.c_str(), "double") == 0) ||
                        (strcmp(ptype.c_str(), "string") == 0);

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
  if (strcmp(ptype.c_str(), "double") == 0)
  {
    r->Set(pname, (double) atof(pvalue.c_str()));
  }
  else if (strcmp(ptype.c_str(), "int") == 0)
  {
    r->Set(pname, (int) atoi(pvalue.c_str()));
  }
  else if (strcmp(ptype.c_str(), "bool") == 0)
  {
    if      (strcmp(pvalue.c_str(), "true" ) == 0) r->Set(pname, true );
    else if (strcmp(pvalue.c_str(), "TRUE" ) == 0) r->Set(pname, true );
    else if (strcmp(pvalue.c_str(), "1"    ) == 0) r->Set(pname, true );
    else if (strcmp(pvalue.c_str(), "false") == 0) r->Set(pname, false);
    else if (strcmp(pvalue.c_str(), "FALSE") == 0) r->Set(pname, false);
    else if (strcmp(pvalue.c_str(), "0"    ) == 0) r->Set(pname, false);
    else { }
  }
  else if (strcmp(ptype.c_str(), "string") == 0)
  {
    r->Set(pname, pvalue);
  }
  else {}
}
//____________________________________________________________________________
void AlgConfigPool::AddRootObjParameter(
                     Registry * r, string ptype, string pname, string pvalue)
{
  // the ROOT object is given in the XML config file as
  // <param> object_name@root_file_name </param>
  vector<string> rootobjv = string_utils::Split(pvalue, "@");

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

  if (strcmp(ptype.c_str(), "TH1F") == 0)
  {
    TH1F * h  = (TH1F*) f.Get(rootobj.c_str());
    if(h) {
      TH1F * ch = new TH1F(*h); // clone
      r->Set(pname,ch);
    } else {
      SLOG("AlgConfigPool", pERROR)
         << "No TH1F named = " << rootobj << " in ROOT file = " << rootfile;
    }
  }
  else if (strcmp(ptype.c_str(), "TH2F") == 0)
  {
    TH2F * h2  = (TH2F*) f.Get(rootobj.c_str());
    if(h2) {
      TH2F * ch2 = new TH2F(*h2); // clone
      r->Set(pname,ch2);
    } else {
      SLOG("AlgConfigPool", pERROR)
         << "No TH2F named = " << rootobj << " in ROOT file = " << rootfile;
    }
  }
  else if (strcmp(ptype.c_str(), "TTree") == 0)
  {
    TTree * t  = (TTree*) f.Get(rootobj.c_str());
    if(t) {
      TTree * ct = new TTree(*t); // clone
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
  return FindRegistry( algorithm->Name(), algorithm->ParamSet() );
}
//____________________________________________________________________________
Registry* AlgConfigPool::FindRegistry(string alg_name, string param_set) const
{
  string key = alg_name + "/" + param_set;

  LOG("AlgConfigPool", pDEBUG) << "Searching for registry with key " << key;

  if( fRegistryPool.count(key) == 1 ) {

     map<string, Registry *>::const_iterator config_entry =
                                                   fRegistryPool.find(key);
     return config_entry->second;

  } else {
     LOG("AlgConfigPool", pWARN) << "No config registry for key " << key;
     return 0;
  }
}
//____________________________________________________________________________
void AlgConfigPool::Print(ostream & stream) const
{
  stream << "\n[-] ALGORITHM CONFIGURATION POOL: ";

  typedef map<string, Registry *>::const_iterator  sregIter;
  typedef map<string, Registry *>::size_type       sregSize;

  sregSize size = fRegistryPool.size();

  stream << size << " configuration sets found" << endl;
  for(sregIter iter = fRegistryPool.begin();
                                       iter != fRegistryPool.end(); iter++) {
     stream << iter->first     << endl;
     stream << *(iter->second) << endl;
  }
}
//____________________________________________________________________________

