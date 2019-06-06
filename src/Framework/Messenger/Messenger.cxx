//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab - June 16, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Jun 01, 2008 - CA
   At Configure(), if the GPRODMODE environmental variable is set then use
   mesg thresholds from Messenger_production.xml rather than Messenger.xml.
   That minimizes verbosity during production jobs.
 @ Dec 06, 2008 - CA
   Adding gAbortingInErr to be set if GENE is exiting in err so as to prevent 
   output clutter (caused by reporting singletons) and help spotting the fatal
   mesg more easily. In PriorityFromString() re-ordered the if-statements, 
   putting the most commonly used priorities on top, so to improve performance.
 @ Aug 25, 2009 - RH
   Use the GetXMLFilePath() to search the potential XML config file locations
   and return the first actual file that can be found. Adapt code to use the
   utils::xml namespace.
 @ Jan 31, 2013 - CA
   The $GMSGCONF var is no longer used. Instead, call 
   Messenger::SetPrioritiesFromXmlFile(string filename) explicitly.

*/
//____________________________________________________________________________

#include <iostream>
#include <vector>
#include <iomanip>

#include "libxml/parser.h"
#include "libxml/xmlmemory.h"
#include "log4cpp/SimpleLayout.hh"

#include <TSystem.h>

#include "Framework/Messenger/Messenger.h"
#include "Framework/Utils/StringUtils.h"
#include "Framework/Utils/PrintUtils.h"
#include "Framework/Utils/XmlParserUtils.h"

using std::setw;
using std::setfill;
using std::cout;
using std::endl;
using std::vector;

using namespace genie;

bool genie::gAbortingInErr = false;

//____________________________________________________________________________
Messenger * Messenger::fInstance = 0;
//____________________________________________________________________________
Messenger::Messenger()
{
  fInstance =  0;
}
//____________________________________________________________________________
Messenger::~Messenger()
{
  fInstance = 0;
}
//____________________________________________________________________________
Messenger * Messenger::Instance()
{
  if(fInstance == 0) {

    // the first thing that get's printed in a GENIE session is the banner
    utils::print::PrintBanner();
	
    static Messenger::Cleaner cleaner;
    cleaner.DummyMethodAndSilentCompiler();

    fInstance = new Messenger;

    log4cpp::Appender * appender;
    appender = new log4cpp::OstreamAppender("default", &cout);
    const char* layoutenv = gSystem->Getenv("GMSGLAYOUT");
    std::string layoutstr = (layoutenv) ? string(layoutenv) : "BASIC";
    if ( layoutstr == "SIMPLE" ) 
      appender->setLayout(new log4cpp::SimpleLayout());
    else
      appender->setLayout(new log4cpp::BasicLayout());


    log4cpp::Category & MSG = log4cpp::Category::getRoot();

    MSG.setAdditivity(false);
    MSG.addAppender(appender);

    fInstance->Configure(); // set user-defined priority levels
  }
  return fInstance;
}
//____________________________________________________________________________
log4cpp::Category & Messenger::operator () (const char * stream)
{
  log4cpp::Category & MSG = log4cpp::Category::getInstance(stream);

  return MSG;
}
//____________________________________________________________________________
void Messenger::SetPriorityLevel(
   const char * stream, log4cpp::Priority::Value priority)
{
  log4cpp::Category & MSG = log4cpp::Category::getInstance(stream);

  MSG.setPriority(priority);
}
//____________________________________________________________________________
void Messenger::Configure(void)
{
// The Configure() method will look for priority level xml config files, read
// the priority levels and set them.
// By default, first it will look for the $GENIE/config/Messenger.xml file.
// Look at this file for the  XML schema.

  string msg_config_file = utils::xml::GetXMLFilePath("Messenger.xml");

  // parse & set the default priority levels
  bool ok = this->SetPrioritiesFromXmlFile(msg_config_file);
  if(!ok) {
    SLOG("Messenger", pERROR)
      << "Priority levels from: " << msg_config_file << " were not set!";
  }

}
//____________________________________________________________________________
bool Messenger::SetPrioritiesFromXmlFile(string filenames)
{
// Reads input XML config file and sets the priority levels.
// The input can be a colection of XML files, delimited with a ":".
// The full path for each file should be given.
// In case of multiple files, all files are read.
// The later each file is, the higher priority it has (if the same mesg 
// stream is listed twice with conflicting priority, the last one is used.).
//

  if(filenames.size()==0) return false;

  vector<string> filename_vec = utils::str::Split(filenames, ":");
  if(filename_vec.size() == 0) return false;

  vector<string>::const_iterator filename_iter;
  for(filename_iter  = filename_vec.begin();
      filename_iter != filename_vec.end(); ++filename_iter) 
  {
     // look in all known XML locations, just like primary file
     string filename = utils::xml::GetXMLFilePath(*filename_iter);

    SLOG("Messenger", pNOTICE)
      << "Reading msg stream priorities from XML file: " << filename;

    xmlDocPtr xml_doc = xmlParseFile(filename.c_str());
    if(xml_doc==NULL) {
       SLOG("Messenger", pERROR)
           << "XML file could not be parsed! [file: " << filename << "]";
       return false;
    }

    xmlNodePtr xml_root = xmlDocGetRootElement(xml_doc);
    if(xml_root==NULL) {
       SLOG("Messenger", pERROR)
           << "XML doc. has null root element! [file: " << filename << "]";
       xmlFreeDoc(xml_doc);
       return false;
    }

    if( xmlStrcmp(xml_root->name, (const xmlChar *) "messenger_config") ) {
       SLOG("Messenger", pERROR)
         << "XML doc. has invalid root element! [file: " << filename << "]";
       xmlFreeNode(xml_root);
       xmlFreeDoc(xml_doc);
       return false;
    }

    xmlNodePtr xml_msgp = xml_root->xmlChildrenNode; // <priority>'s

    // loop over all xml tree nodes that are children of the <spline> node
    while (xml_msgp != NULL) {

      // enter everytime you find a <priority> tag
      if( (!xmlStrcmp(xml_msgp->name, (const xmlChar *) "priority")) ) {

         string msgstream = utils::str::TrimSpaces(
                  utils::xml::GetAttribute(xml_msgp, "msgstream"));
         string priority =
                utils::xml::TrimSpaces( xmlNodeListGetString(
                               xml_doc, xml_msgp->xmlChildrenNode, 1));
         log4cpp::Priority::Value pv = this->PriorityFromString(priority);
         this->SetPriorityLevel(msgstream.c_str(), pv);
         SLOG("Messenger", pINFO)
                  << "Set priority level: " << setfill('.')
                          << setw(24) << msgstream << " --> " << priority;
      }
      xml_msgp = xml_msgp->next;
    }//xml_msgp != NULL

    //xmlFree(xml_msgp);
    xmlFreeNode(xml_msgp);
    xmlFreeDoc(xml_doc);
  }//filename_iter
  
  return true;
}
//____________________________________________________________________________
log4cpp::Priority::Value Messenger::PriorityFromString(string p)
{
  if ( p.find("DEBUG")  != string::npos ) return log4cpp::Priority::DEBUG;
  if ( p.find("INFO")   != string::npos ) return log4cpp::Priority::INFO;
  if ( p.find("NOTICE") != string::npos ) return log4cpp::Priority::NOTICE;
  if ( p.find("WARN")   != string::npos ) return log4cpp::Priority::WARN;
  if ( p.find("ERROR")  != string::npos ) return log4cpp::Priority::ERROR;
  if ( p.find("CRIT")   != string::npos ) return log4cpp::Priority::CRIT;
  if ( p.find("ALERT")  != string::npos ) return log4cpp::Priority::ALERT;
  if ( p.find("FATAL")  != string::npos ) return log4cpp::Priority::FATAL;

  SLOG("Messenger", pWARN)
                    << "Unknown priority = " << p << " - Setting to INFO";
  return log4cpp::Priority::INFO;
}
//____________________________________________________________________________
