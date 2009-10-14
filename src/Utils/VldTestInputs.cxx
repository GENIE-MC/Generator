//____________________________________________________________________________
/*
 Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Oct 13, 2009 - CA
   This class was first added in version 2.5.1
*/
//____________________________________________________________________________

#include <cstdlib>

#include "libxml/parser.h"
#include "libxml/xmlmemory.h"
#include "libxml/xmlreader.h"

#include "Messenger/Messenger.h"
#include "Utils/VldTestInputs.h"
#include "Utils/StringUtils.h"

using namespace genie;

//____________________________________________________________________________
VldTestInputs::VldTestInputs(const int nmaxmodels)
{
  Init(nmaxmodels);
}
//____________________________________________________________________________
VldTestInputs::~VldTestInputs(void)
{

}
//____________________________________________________________________________
int VldTestInputs::NModels(void) const
{
  return fNModels;
}
//____________________________________________________________________________
string VldTestInputs::ModelTag(int imodel) const
{
  return (*fModelTag)[imodel];
}
//____________________________________________________________________________
TFile*  VldTestInputs::XSecFile(int imodel) const
{
  return (*fXSecFile)[imodel];
}
//____________________________________________________________________________
TChain* VldTestInputs::EvtChain(int imodel) const
{
  return (*fEvtChain)[imodel];
}
//____________________________________________________________________________
bool VldTestInputs::LoadFromFile(string xmlfile)
{
  LOG("VldTestInputs", pNOTICE) << "Loading: " << xmlfile;

  vector<string>  & model_tag  = *fModelTag;
  vector<TFile*>  & xsec_file  = *fXSecFile;
  vector<TChain*> & evt_chain  = *fEvtChain;

  xmlTextReaderPtr reader = xmlNewTextReaderFilename(xmlfile.c_str());

  const int kNodeTypeStartElement = 1;
//const int kNodeTypeEndElement   = 15;

  int  imodel       = -1;
  bool is_xsec_file = false;
  bool is_evt_file  = false;

  if (reader != NULL) {
     int ret = xmlTextReaderRead(reader);
     while (ret == 1) {
       xmlChar * name  = xmlTextReaderName     (reader);
       xmlChar * value = xmlTextReaderValue    (reader);
       int       type  = xmlTextReaderNodeType (reader);
       int       depth = xmlTextReaderDepth    (reader);

       bool start_element = (type==kNodeTypeStartElement);

       if(depth==0 && start_element) {
         LOG("VldTestInputs", pDEBUG) << "Root element = " << name;
         if(xmlStrcmp(name, (const xmlChar *) "vld_inputs")) {
           LOG("VldTestInputs", pERROR)
             << "\nXML doc. has invalid root element! [filename: " 
             << xmlfile << "]";
           return false;
         }
       }
       if( (!xmlStrcmp(name, (const xmlChar *) "model")) && start_element) {
         imodel++;
         is_xsec_file      = false;
         is_evt_file       = false;
         xmlChar * xname = xmlTextReaderGetAttribute(reader,(const xmlChar*)"name");
         string sname    = utils::str::TrimSpaces((const char *)xname);
         model_tag[imodel] = sname;
         LOG("VldTestInputs", pNOTICE) 
            << "Adding files for model ID: " 
            << imodel << " (" << model_tag[imodel] << ")";
         xmlFree(xname);
       }
       if( (!xmlStrcmp(name, (const xmlChar *) "xsec_file")) && start_element) {
         is_xsec_file = true;
         is_evt_file  = false;
       }
       if( (!xmlStrcmp(name, (const xmlChar *) "evt_file")) && start_element) {
         is_xsec_file = false;
         is_evt_file  = true;
       }

       if( (!xmlStrcmp(name, (const xmlChar *) "#text")) && depth==3) {
         string filename = utils::str::TrimSpaces((const char *)value);
         if(is_evt_file) {
           LOG("VldTestInputs", pNOTICE) 
                << " * Adding event file: " << filename;
           if(!evt_chain[imodel]) {
              evt_chain[imodel] = new TChain("gst");
           }
           evt_chain[imodel]->Add(filename.c_str());
         }
         if(is_xsec_file) {
           LOG("VldTestInputs", pNOTICE) 
                << " * Adding cross section file: " << filename;
           xsec_file[imodel] = new TFile(filename.c_str(), "recreate");
           if(!xsec_file[imodel]) {
             exit(1);
           }
         }
       }

       xmlFree(name);
       xmlFree(value);
       ret = xmlTextReaderRead(reader);
     }//ret==1

     xmlFreeTextReader(reader);

  }//reader!=null

  fNModels = imodel+1;

  return true;
}
//____________________________________________________________________________
void VldTestInputs::Init(const int nmaxmodels)
{
  fModelTag = new vector<string>  (nmaxmodels);
  fXSecFile = new vector<TFile*>  (nmaxmodels);
  fEvtChain = new vector<TChain*> (nmaxmodels);

  for(int i=0; i<nmaxmodels; i++) {
    (*fModelTag) [i] = "";
    (*fXSecFile) [i] = 0;
    (*fEvtChain) [i] = 0;
  }
}
//____________________________________________________________________________

