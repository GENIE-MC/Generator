//____________________________________________________________________________
/*!

\class    genie::FermiMomentumTablePool

\brief    Singleton class to load & serve tables of Fermi momentum constants

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  August 18, 2005

*/
//____________________________________________________________________________

#include <iostream>
#include <string>

#include "libxml/parser.h"
#include "libxml/xmlmemory.h"

#include <TSystem.h>

#include "Messenger/Messenger.h"
#include "Nuclear/FermiMomentumTablePool.h"
#include "Nuclear/FermiMomentumTable.h"
#include "Utils/StringUtils.h"
#include "Utils/XmlParserUtils.h"

using std::string;
using std::cout;
using std::endl;

using namespace genie;

//____________________________________________________________________________
FermiMomentumTablePool * FermiMomentumTablePool::fInstance = 0;
//____________________________________________________________________________
FermiMomentumTablePool::FermiMomentumTablePool()
{
  if( ! this->LoadTables() ) {
   LOG("FermiP", pERROR) << "FermiMomentumTablePool initialization failed!";
  }
  fInstance =  0;
}
//____________________________________________________________________________
FermiMomentumTablePool::~FermiMomentumTablePool()
{
  cout << "FermiMomentumTablePool singleton dtor: "
                               << "Deleting all Fermi momenta tables" << endl;
  map<string, FermiMomentumTable *>::iterator titer;
  for(titer = fKFSets.begin(); titer != fKFSets.end(); ++titer) {
    FermiMomentumTable * t = titer->second;
    if(t) delete t;
    t=0;
  }
  fKFSets.clear();
  fInstance = 0;
}
//____________________________________________________________________________
FermiMomentumTablePool * FermiMomentumTablePool::Instance()
{
  if(fInstance == 0) {

    LOG("FermiP", pINFO) << "FermiMomentumTablePool late initialization";

    static FermiMomentumTablePool::Cleaner cleaner;
    cleaner.DummyMethodAndSilentCompiler();

    fInstance = new FermiMomentumTablePool;
  }
  return fInstance;
}
//____________________________________________________________________________
const FermiMomentumTable * FermiMomentumTablePool::GetTable(string name)
{
  string defopt = "Default"; // default option

  map<string, FermiMomentumTable *>::const_iterator table_iter;
  FermiMomentumTable * table = 0;

  if(fKFSets.count(name) == 1) {
     table_iter = fKFSets.find(name);
     table = table_iter->second;
     if(table) return table;
  }
  if(fKFSets.count(defopt) == 1) {
     LOG("FermiP", pWARN)
         << "Fermi momentum table: [" << name << "] was not found! "
                               << "Switching to table: [" << defopt << "]";
     table_iter = fKFSets.find(defopt);
     table = table_iter->second;
     if(table) return table;
  }
  LOG("FermiP", pERROR)
            << "Not even the default Fermi momentum table was not found!";
  return 0;
}
//____________________________________________________________________________
bool FermiMomentumTablePool::LoadTables(void)
{
  bool loaded = true;

  //-- get base GENIE config directory from the environment
  //   (search for $GALGCONF or use the default: $GENIE/config)
  string config_dir = (gSystem->Getenv("GALGCONF")) ?
            string(gSystem->Getenv("GALGCONF")) :
            string(gSystem->Getenv("GENIE")) + string("/config");

  //-- Build the path Fermi momenta sets XML file
  string filename = config_dir + string("/FermiMomentumTables.xml");

  LOG("FermiP", pINFO)  << "Loading Fermi momenta from file: " << filename;

  bool is_accessible = ! (gSystem->AccessPathName( filename.c_str() ));

  if ( is_accessible ) {
      XmlParserStatus_t status = this->ParseXMLTables(filename);
      if(status != kXmlOK) {
         string sst = XmlParserStatus::AsString(status);
         LOG("FermiP", pWARN)
                  << "XML parser status: " << sst
                                 << " - Couldn't read file: " << filename;
         loaded = false;
      }
  } else {
      LOG("FermiP", pWARN) << "Not accessible file: " << filename;
      loaded = false;
  }
  return loaded;
};
//____________________________________________________________________________
XmlParserStatus_t FermiMomentumTablePool::ParseXMLTables(string filename)
{
  LOG("FermiP", pDEBUG) << "Reading XML file: " << filename;

  xmlDocPtr xml_doc = xmlParseFile(filename.c_str());
  if(xml_doc == NULL) return kXmlNotParsed;

  xmlNodePtr xml_root = xmlDocGetRootElement(xml_doc);
  if(xml_root==NULL) return kXmlEmpty;

  const xmlChar * xml_root_name = (const xmlChar *)"fermi_momentum_const";
  if( xmlStrcmp(xml_root->name, xml_root_name) ) return kXmlInvalidRoot;

  xmlNodePtr xml_kft = xml_root->xmlChildrenNode; // <kf_table>'s

  // loop over <kf_table> nodes
  while (xml_kft != NULL) {
    if( (!xmlStrcmp(xml_kft->name, (const xmlChar *) "kf_table")) ) {

       string name = utils::str::TrimSpaces(
                      XmlParserUtils::GetAttribute(xml_kft, "name"));

       LOG("FermiP", pDEBUG) << "Reading Fermi momenta table: " << name;

       FermiMomentumTable * kftable = new FermiMomentumTable; // new table

       // loop over <kf> nodes
       xmlNodePtr xml_kf = xml_kft->xmlChildrenNode; // <kf_table> children
       while (xml_kf != NULL) {
         if( (!xmlStrcmp(xml_kf->name, (const xmlChar *) "kf")) ) {

           string spdgc = utils::str::TrimSpaces(
                       XmlParserUtils::GetAttribute(xml_kf, "nucleus_pdgc"));
           int    pdgc = atoi (spdgc.c_str());

           xmlNodePtr xml_cur = xml_kf->xmlChildrenNode; // <kf> children
           const xmlChar * ntag = (const xmlChar *)"n";
           const xmlChar * ptag = (const xmlChar *)"p";
           double kfp=0, kfn=0, kf=0;
           // loop over <kf> children
           while (xml_cur != NULL) {
             bool isp = !xmlStrcmp(xml_cur->name, ptag);
             bool isn = !xmlStrcmp(xml_cur->name, ntag);
             if(isn || isp) {
               string skf = XmlParserUtils::TrimSpaces(
                  xmlNodeListGetString(xml_doc, xml_cur->xmlChildrenNode, 1));
               kf = atof(skf.c_str());
             }
             if(isp) kfp = kf;
             if(isn) kfn = kf;
             xml_cur = xml_cur->next;
           }
           xmlFree(xml_cur);

           KF_t kft;
           kft.p = kfp;
           kft.n = kfn;

           LOG("FermiP", pDEBUG)
              << "Add KF table entry: PDGC = " << pdgc
                  << " --> " << "kf(p) = " << kft.p << ", kf(n) = " << kft.n;
           kftable->AddTableEntry(pdgc,kft);
         } //<x> == <kf>
         xml_kf = xml_kf->next;
      } //<kf> loop
      xmlFree(xml_kf);

      fKFSets.insert(
                map<string, FermiMomentumTable *>::value_type(name,kftable));

    } //<x> == <kf_table>
    xml_kft = xml_kft->next;
  } //<kf_table> loop
  xmlFree(xml_kft);
  xmlFree(xml_root);

  return kXmlOK;
}
//____________________________________________________________________________

