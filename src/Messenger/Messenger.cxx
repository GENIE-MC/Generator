//____________________________________________________________________________
/*
 Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory - June 16, 2004

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
*/
//____________________________________________________________________________

#include <iostream>
#include <vector>
#include <iomanip>

#include "libxml/parser.h"
#include "libxml/xmlmemory.h"

#include <TSystem.h>

#include "Messenger/Messenger.h"
#include "Utils/StringUtils.h"
#include "Utils/PrintUtils.h"
#include "Utils/XmlParserUtils.h"

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
// Clean up. Don't clutter output if exiting in err.

  if(!gAbortingInErr) {
     cout << "Messenger singleton dtor" << endl;
  }
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
// By default, first it will look for the $GENIE/config/messenger.xml file.
// Then it will look any messenger configuration xml file defined in the
// GMSGCONF env variable. This variable can contain more than one XML files
// (that should be delimited with a ';'). The full path for each file should
// be given. See the $GENIE/config/messenger.xml for the XML schema.
// The later each file is, the higher priority it has - eg. if the same stream
// is listed twice with conflicting priority then the one found last is used.

  bool ok = false;

  string filename = gSystem->Getenv("GPRODMODE") ? 
                    "Messenger_production.xml" : "Messenger.xml";

  string msg_config_file = utils::xml::GetXMLFilePath(filename);

  // parse & set the default priority levels
  ok = this->SetPrioritiesFromXmlFile(msg_config_file);
  if(!ok) {
    SLOG("Messenger", pERROR)
           << "Priority levels from: " << msg_config_file << " were not set!";
  }

  //-- checkout the GMSGCONF conf for additional messenger configuration files
  string gmsgconf = (gSystem->Getenv("GMSGCONF") ?
                                            gSystem->Getenv("GMSGCONF") : "");
  SLOG("Messenger", pINFO) << "$GMSGCONF env.var = " << gmsgconf;

  if(gmsgconf.size()>0) {
     //-- check for multiple files delimited with a ":"
     vector<string> conf_xmlv = utils::str::Split(gmsgconf, ":");

     //-- loop over messenger config files -- parse & set priorities
     vector<string>::const_iterator conf_iter;
     for(conf_iter = conf_xmlv.begin();
                                 conf_iter != conf_xmlv.end(); ++conf_iter) {
          string conf_xml = *conf_iter;
          ok = this->SetPrioritiesFromXmlFile(conf_xml);
          if(!ok) {
            SLOG("Messenger", pERROR)
                << "Priority levels from: " << conf_xml << " were not set!";
            }
     }
  } else {
    SLOG("Messenger", pINFO)
                  << "No additional messenger config XML file was specified";
  }
}
//____________________________________________________________________________
bool Messenger::SetPrioritiesFromXmlFile(string filename)
{
// Reads the XML config file and sets the priority levels
//
  SLOG("Messenger", pINFO)
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
    return false;
  }

  if( xmlStrcmp(xml_root->name, (const xmlChar *) "messenger_config") ) {
    SLOG("Messenger", pERROR)
      << "XML doc. has invalid root element! [file: " << filename << "]";
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
         SLOG("Messenger", pINFO)
                  << "Setting priority level: " << setfill('.')
                          << setw(24) << msgstream << " --> " << priority;

         log4cpp::Priority::Value pv = this->PriorityFromString(priority);
         this->SetPriorityLevel(msgstream.c_str(), pv);
     }
     xml_msgp = xml_msgp->next;
  }
  xmlFree(xml_msgp);

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
