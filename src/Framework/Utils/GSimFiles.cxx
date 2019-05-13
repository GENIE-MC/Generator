//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Laboratory

 For the class documentation see the corresponding header file.
*/
//____________________________________________________________________________

#include <cstdlib>

#include "libxml/parser.h"
#include "libxml/xmlmemory.h"
#include "libxml/xmlreader.h"

#include "Framework/Messenger/Messenger.h"
#include "Framework/Utils/GSimFiles.h"
#include "Framework/Utils/StringUtils.h"

using std::endl;

using namespace genie;

//____________________________________________________________________________
namespace genie {
  ostream & operator << (ostream & stream, const GSimFiles & f)
  {
    f.Print(stream);
    return stream;
  }
}
//____________________________________________________________________________
GSimFiles::GSimFiles(bool chain, const int nmaxmodels)
{
  fDoChain = chain;
  this->Init(nmaxmodels);
}
//____________________________________________________________________________
GSimFiles::~GSimFiles(void)
{
  this->CleanUp();
}
//____________________________________________________________________________
int GSimFiles::NModels(void) const
{
  return fNModels;
}
//____________________________________________________________________________
int GSimFiles::FindModelID(string tag) const
{
  int imodel = 0;
  vector<string>::const_iterator it = fModelTag->begin();
  for( ; it != fModelTag->end(); ++it) {
    if(*it == tag) return imodel;
    imodel++;
  }
  return -1;
}
//____________________________________________________________________________
string GSimFiles::ModelTag(int imodel) const
{
  return (*fModelTag)[imodel];
}
//____________________________________________________________________________
TFile* GSimFiles::XSecFile(int imodel) const
{
  return (*fXSecFile)[imodel];
}
//____________________________________________________________________________
string GSimFiles::XSecFileName(int imodel) const
{
  return (*fXSecFileName)[imodel];
}
//____________________________________________________________________________
TChain* GSimFiles::EvtChain(int imodel) const
{
  return (*fEvtChain)[imodel];
}
//____________________________________________________________________________
vector<string> & GSimFiles::EvtFileNames(int imodel) const
{
  return (*fEvtFileNames)[imodel];
}
//____________________________________________________________________________
const string & GSimFiles::PathToXMLFile(void) const
{
   return fPath2XMLFile;
}
//____________________________________________________________________________
bool GSimFiles::LoadFromFile(string xmlfile)
{
  LOG("GSimFiles", pNOTICE) << "Loading: " << xmlfile;

  vector<string>  &          model_tag      = *fModelTag;
  vector<TFile*>  &          xsec_file      = *fXSecFile;
  vector<string>  &          xsec_filename  = *fXSecFileName;
  vector<TChain*> &          evt_chain      = *fEvtChain;
  vector< vector<string> > & evt_filenames  = *fEvtFileNames;
 
  xmlTextReaderPtr reader = xmlNewTextReaderFilename(xmlfile.c_str());
  if(reader == NULL) {
    return false;
  }

  const int kNodeTypeStartElement = 1;
  const int kNodeTypeEndElement   = 15;

  int  imodel            = 0;
  bool is_xsec_file      = false;
  bool is_evt_file       = false;
  bool is_ghep_evt_file  = false;
  bool is_gst_evt_file   = false;
  bool have_ghep_files   = false;
  bool have_gst_files    = false;

  if (reader != NULL) {
     int ret = xmlTextReaderRead(reader);
     while (ret == 1) {
       xmlChar * name  = xmlTextReaderName     (reader);
       xmlChar * value = xmlTextReaderValue    (reader);
       int       type  = xmlTextReaderNodeType (reader);
       int       depth = xmlTextReaderDepth    (reader);

       bool start_element = (type==kNodeTypeStartElement);
       bool end_element   = (type==kNodeTypeEndElement);

       if(depth==0 && start_element) {
         LOG("GSimFiles", pDEBUG) << "Root element = " << name;
         if(xmlStrcmp(name, (const xmlChar *) "genie_simulation_outputs")) {
           LOG("GSimFiles", pERROR)
             << "\nXML doc. has invalid root element! [filename: " 
             << xmlfile << "]";
           return false;
         }
       }

       if( (!xmlStrcmp(name, (const xmlChar *) "model")) && start_element) {
         xmlChar * xname = xmlTextReaderGetAttribute(reader,(const xmlChar*)"name");
         string sname    = utils::str::TrimSpaces((const char *)xname);
         model_tag[imodel] = sname;
         LOG("GSimFiles", pNOTICE) 
            << "Adding files for model ID: " 
            << imodel << " (" << model_tag[imodel] << ")";
         xmlFree(xname);
       }
       if( (!xmlStrcmp(name, (const xmlChar *) "model")) && end_element) {
         LOG("GSimFiles", pNOTICE) 
            << "Done adding files for model ID: " << imodel;
         imodel++;
       }
       if( (!xmlStrcmp(name, (const xmlChar *) "xsec_file")) && start_element) {
         is_xsec_file = true;
       }
       if( (!xmlStrcmp(name, (const xmlChar *) "xsec_file")) && end_element) {
         is_xsec_file = false;
       }
       if( (!xmlStrcmp(name, (const xmlChar *) "evt_file")) && start_element) {
         is_evt_file       = true;
         is_ghep_evt_file  = false;
         is_gst_evt_file   = false;
         xmlChar * xfmt = xmlTextReaderGetAttribute(reader,(const xmlChar*)"format");
         string sfmt = utils::str::TrimSpaces((const char *)xfmt);
         if (sfmt.find("gst")  != string::npos) 
         { 
            is_gst_evt_file  = true; 
            if(!have_gst_files) { have_gst_files = true; }
         }
         else 
         if (sfmt.find("ghep") != string::npos) 
         { 
            is_ghep_evt_file = true; 
            if(!have_ghep_files) { have_ghep_files = true; }
         } 
         if(have_gst_files && have_ghep_files) {
            LOG("GSimFiles", pFATAL) 
               << "Oops! You shouldn't mix GHEP and GST event files in GSimFiles";
            LOG("GSimFiles", pFATAL) 
               << "Please correct XML file: " << xmlfile;
            gAbortingInErr = true;;
            exit(1);
         }
         xmlFree(xfmt);
       }
       if( (!xmlStrcmp(name, (const xmlChar *) "evt_file")) && end_element) {
         is_evt_file  = false;
       }
       if( (!xmlStrcmp(name, (const xmlChar *) "#text")) && depth==3) {
         string filename = utils::str::TrimSpaces((const char *)value);
         if(is_evt_file) {
           LOG("GSimFiles", pNOTICE) 
                << " * Adding event file: " << filename;
           // chain the event trees, if requested
           if(fDoChain) {
             if(!evt_chain[imodel] && is_gst_evt_file) {
                 evt_chain[imodel] = new TChain("gst");
             } else 
             if(!evt_chain[imodel] && is_ghep_evt_file) {
                 evt_chain[imodel] = new TChain("gtree");
             }
             if(evt_chain[imodel]) {
                evt_chain[imodel]->Add(filename.c_str());
             }
           }//chain?
           evt_filenames[imodel].push_back(filename);
         }
         if(is_xsec_file) {
           LOG("GSimFiles", pNOTICE) 
                << " * Adding cross section file: " << filename;
           xsec_file    [imodel] = new TFile(filename.c_str(), "read");
           xsec_filename[imodel] = filename;
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

  fNModels = imodel;
  
  fPath2XMLFile = xmlfile;

  return true;
}
//____________________________________________________________________________
void GSimFiles::Print(ostream & stream) const
{
  stream << endl;
  stream << "loaded from path: " << fPath2XMLFile << endl;
  for(int imodel=0; imodel < this->NModels(); imodel++) {
     stream << "model tag: [" << this->ModelTag(imodel) << "]" << endl;
     if(this->XSecFile(imodel)) {
        stream << "   xsec file  : " << this->XSecFileName(imodel) << endl;
     }
     const vector<string> & filenames = this->EvtFileNames(imodel);  
     vector<string>::const_iterator iter = filenames.begin();
     for( ; iter != filenames.end(); ++iter) {
        string filename = *iter;
        stream << "   event file : " << filename << endl;
     }
  }
}
//____________________________________________________________________________
void GSimFiles::Init(const int nmaxmodels)
{
  fNModels      = 0;
  fModelTag     = new vector<string>          (nmaxmodels);
  fXSecFile     = new vector<TFile*>          (nmaxmodels);
  fXSecFileName = new vector<string>          (nmaxmodels);
  fEvtChain     = new vector<TChain*>         (nmaxmodels);
  fEvtFileNames = new vector<vector<string> > (nmaxmodels);

  for(int i=0; i<nmaxmodels; i++) {
    (*fModelTag) [i] = "";
    (*fXSecFile) [i] = 0;
    (*fEvtChain) [i] = 0;
  }

  fPath2XMLFile = "";

}
//____________________________________________________________________________
void GSimFiles::CleanUp(void)
{

   fPath2XMLFile = "";

}
//____________________________________________________________________________

