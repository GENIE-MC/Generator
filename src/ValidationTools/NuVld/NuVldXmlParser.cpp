//____________________________________________________________________________
/*
 Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory - May 06, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Aug 25, 2009 - CA
   Removed redundant versions of ParserUtils.h and ParserStatus.h in favor of
   the ones in $GENIE/Conventions and $GENIE/Utils. Updated code accordingly.

*/
//____________________________________________________________________________

#include "Messenger/Messenger.h"
#include "Utils/XmlParserUtils.h"
#include "Utils/StringUtils.h"
#include "ValidationTools/NuVld/NuVldXmlParser.h" 
#include "ValidationTools/NuVld/XmlNuXSecRecord.h"
#include "ValidationTools/NuVld/XmlElDiffXSecRecord.h"
#include "ValidationTools/NuVld/XmlSFRecord.h"

namespace genie {
namespace nuvld {
  
//__________________________________________________________________________
NuVldXmlParser::NuVldXmlParser()
{
  fDataSet = new XmlDataSet;
} 
//__________________________________________________________________________
NuVldXmlParser::~NuVldXmlParser()
{

}
//__________________________________________________________________________
void NuVldXmlParser::ParseXmlDocument(const char * filename)
{
  fXmlFilename = string(filename);   
  fXmlDoc      = xmlParseFile(fXmlFilename.c_str());

  if((fXmlPStatus = VerifyParsing()) == kXmlOK) FillDataSet();
}
//__________________________________________________________________________
const XmlDataSet & NuVldXmlParser::GetDataSet(void) const
{
  return *fDataSet;
}
//__________________________________________________________________________
XmlParserStatus_t NuVldXmlParser::VerifyParsing(void)
{
 SLOG("NuVld", pINFO)  << "NuVldXmlParser::VerifyParsing() :......: ";

 XmlParserStatus_t status;
 
 if(fXmlDoc==NULL) { 
   status = kXmlNotParsed;
   SLOG("NuVld", pERROR)  
     << "***" << XmlParserStatus::AsString(status);
   return status;
 }
 
 xmlNodePtr xml_cur = xmlDocGetRootElement(fXmlDoc);

 if(xml_cur==NULL) {
   status = kXmlEmpty;
   SLOG("NuVld", pERROR)  
     << "***" << XmlParserStatus::AsString(status);
   return status;
 }

 if( xmlStrcmp(xml_cur->name, (const xmlChar *) "nuscat_data") ) {

   status = kXmlInvalidRoot;
   xmlFreeDoc(fXmlDoc);
   SLOG("NuVld", pERROR)  
     << "***" << XmlParserStatus::AsString(status);
   return status;
 }

 status = kXmlOK;
 
 SLOG("NuVld", pINFO) 
   << XmlParserStatus::AsString(status);
 
 return status;
}
//__________________________________________________________________________
void NuVldXmlParser::FillDataSet(void)
{
 SLOG("NuVld", pDEBUG) << "Starts parsing the XML file";
 
 xmlChar *  xml_exp_name = 0;
 xmlNodePtr xml_cur      = xmlDocGetRootElement(fXmlDoc);
 
 xml_cur = xml_cur->xmlChildrenNode;

 // loop over all xml tree nodes that are children of the root node
 while (xml_cur != NULL) {

   // enter everytime you find an 'experiment' tag
   if( (!xmlStrcmp(xml_cur->name, (const xmlChar *) "experiment")) ) {

     // get experiment attributes (name,tag) - the last can be absent
     // if there is only one experiment entry with the given name

     xml_exp_name = xmlGetProp(xml_cur, (const xmlChar *)"name");
     string name  = utils::xml::TrimSpaces( xml_exp_name );

     xmlFree(xml_exp_name);     

     // parse all measurements associated with the given experiment
     SLOG("NuVld", pDEBUG) 
        << "Parsing entries for experiment: " << name;     

     XmlExperimentMeasurements * meas_list = ParseExperiment(xml_cur, name);

     // add the entry to the top-level STL map
     fDataSet->Add( name, meas_list );
   }
   
   xml_cur = xml_cur->next;
 }

 SLOG("NuVld", pDEBUG) << "Finished parsing XML data file";
}
//__________________________________________________________________________
XmlExperimentMeasurements * NuVldXmlParser::ParseExperiment(
                                            xmlNodePtr xml_cur, string name)
{
 XmlExperimentMeasurements * meas_list = new XmlExperimentMeasurements;
 
 fMeasurementTag = 0; 
 
 xml_cur = xml_cur->xmlChildrenNode;

 while (xml_cur != NULL) {

    //
    //  ***  parse 'experiment-info' section 
    //

    if( (!xmlStrcmp(xml_cur->name, (const xmlChar *) "exp_info")) ) {

       SLOG("NuVld", pDEBUG) 
          << "Parsing experiment info";
       
       XmlExperimentInfo * e_info = ParseXmlExperimentInfo(xml_cur);
       e_info->Add("name", name);
       cout << *e_info << endl; // send the experiment info to stdout       
       meas_list->Add( e_info );
    }

    // ***  parse 'beam-spectrum' section 

    if( (!xmlStrcmp(xml_cur->name, (const xmlChar *) "flux_spectrum")) ) {
       
       string beam    = utils::xml::GetAttribute(xml_cur, "beam");       
       string F_units = utils::xml::GetAttribute(xml_cur, "flux_units");
       string E_units = utils::xml::GetAttribute(xml_cur, "E_units");
       string E_frame = utils::xml::GetAttribute(xml_cur, "E_frame");

       SLOG("NuVld", pDEBUG) << "Parsing flux spectrum for beam: " << beam;

       XmlBeamFluxSpectrum * flux = ParseXmlBeamFluxSpectrum(xml_cur);
       
       flux->SetFluxUnits(F_units);
       flux->SetEnergyUnits(E_units);
       flux->SetEnergyFrame(E_frame);

       cout << *flux << endl;
       
       meas_list->Add( beam, flux );
    }

    //
    // *** parse 'measurement' section 
    //
    if( (!xmlStrcmp(xml_cur->name, (const xmlChar *) "measurement")) ) {

       // understand what kind of measurement this is...
       string obs = utils::xml::GetAttribute(xml_cur, "observable");
       XmlObservable_t obst = XmlObservable::GetXmlObservable(obs); 

       // read & add this measurement
       XmlMeasurement * meas = ParseXmlMeasurement(xml_cur, obst);       
       meas_list->Add( meas );
    }

    xml_cur = xml_cur->next; // move to next xml node
    
 } // loop over xml nodes
 
 return meas_list;
}
//__________________________________________________________________________
XmlExperimentInfo * NuVldXmlParser::ParseXmlExperimentInfo(xmlNodePtr xml_cur)
{
 xmlChar *  xml_string;
 string     std_string;

 XmlExperimentInfo * e_info = new XmlExperimentInfo;
 
 string comment = utils::xml::TrimSpaces( 
             xmlNodeListGetString(fXmlDoc, xml_cur->xmlChildrenNode, 1) );

 e_info->Add("comment", utils::str::FilterString(",",comment));

 xml_cur = xml_cur->xmlChildrenNode;
  
 while (xml_cur != NULL) {
    for(int itag = 0; itag < c_exp_info_ntags; itag++) {

       // the element names are defined in XmlExperimentInfo.h       
       if( (!xmlStrcmp(xml_cur->name, 
                       (const xmlChar *) c_exp_info_tag[itag].c_str())) ) {

          xml_string = xmlNodeListGetString(fXmlDoc, xml_cur->xmlChildrenNode, 1); 
          std_string = utils::xml::TrimSpaces(xml_string);

          SLOG("NuVld", pINFO)
                << "Adding " << c_exp_info_tag[itag] << ": " << std_string;
                
          e_info->Add(c_exp_info_tag[itag], 
                      utils::str::FilterString(",",std_string)); 

          xmlFree(xml_string);     
        }
    }
    xml_cur = xml_cur->next;
 } // over xml nodes
 
 return e_info;
}
//__________________________________________________________________________
XmlBeamFluxSpectrum * 
      NuVldXmlParser::ParseXmlBeamFluxSpectrum(xmlNodePtr xml_cur)
{
  XmlBeamFluxSpectrum * spectrum = new XmlBeamFluxSpectrum;

  xml_cur = xml_cur->xmlChildrenNode;

  while (xml_cur != NULL) {
    if( (!xmlStrcmp(xml_cur->name, (const xmlChar *) "bin")) ) {
    
        XmlBeamFluxBin * spectrum_bin = ParseXmlBeamFluxBin(xml_cur);
        spectrum->Add( spectrum_bin );
    }
    xml_cur = xml_cur->next;
  }  
  return spectrum;
}
//__________________________________________________________________________
XmlBeamFluxBin * NuVldXmlParser::ParseXmlBeamFluxBin(xmlNodePtr xml_cur)
{
  xml_cur = xml_cur->xmlChildrenNode;

  xmlChar *  xml_string;
  string     std_string;

  vector<string> energies;
  vector<string> fluxes;
  
  while (xml_cur != NULL) {
    if( (!xmlStrcmp(xml_cur->name, (const xmlChar *) "E")) ) {
        xml_string = xmlNodeListGetString(fXmlDoc, xml_cur->xmlChildrenNode, 1); 
        std_string = utils::xml::TrimSpaces(xml_string);

        energies = utils::str::Split( std_string, string(",") );
    }
    if( (!xmlStrcmp(xml_cur->name, (const xmlChar *) "Flux")) ) {
        xml_string = xmlNodeListGetString(fXmlDoc, xml_cur->xmlChildrenNode, 1); 
        std_string = utils::xml::TrimSpaces(xml_string);

        fluxes = utils::str::Split( std_string, string(",") );
    }
    xml_cur = xml_cur->next;
  }  
  XmlBeamFluxBin * spectrum_bin = new XmlBeamFluxBin(
       energies[0], energies[1], energies[2], 
       fluxes[0],   fluxes[1],   fluxes[2]   );

  return spectrum_bin;
}
//__________________________________________________________________________
XmlMeasurement * NuVldXmlParser::ParseXmlMeasurement(
     xmlNodePtr xml_cur, XmlObservable_t obs)
{
  SLOG("NuVld", pINFO) 
    << "Start parsing a <<" 
    << XmlObservable::AsString(obs) << ">> measurement";

  xml_cur = xml_cur->xmlChildrenNode;

  XmlMeasurement *       meas   = new XmlMeasurement;
  XmlMeasurementHeader * header = 0;

  while (xml_cur != NULL) {

    // parse the measurement header
    if( (!xmlStrcmp(xml_cur->name, (const xmlChar *) "header")) ) {
        header = ParseXmlMeasurementHeader(xml_cur);

        header->Add("observable", XmlObservable::Code(obs));
        header->Add("tag",        utils::str::IntAsString(fMeasurementTag));

        meas->Add( header ); 
        cout << *header << endl;
    }

    // parse every occurence of 'point' within the current measurement
    if( (!xmlStrcmp(xml_cur->name, (const xmlChar *) "point")) ) {
        XmlRecordBase * rec = ParsePoint(xml_cur, obs);  
        cout << *rec;
        meas->Add(rec);
    }
    xml_cur = xml_cur->next;
  }
  fMeasurementTag++;
  
  return meas;
}
//__________________________________________________________________________
XmlMeasurementHeader * NuVldXmlParser::ParseXmlMeasurementHeader(xmlNodePtr xml_cur)
{
  XmlMeasurementHeader * header = new XmlMeasurementHeader;

  string comment = utils::xml::TrimSpaces( 
             xmlNodeListGetString(fXmlDoc, xml_cur->xmlChildrenNode, 1) );

  header->Add("comment", utils::str::FilterString(",",comment));

  xml_cur = xml_cur->xmlChildrenNode;

  while (xml_cur != NULL) {
   // Add a reference that is listed in relation with this measurement.
   // There is no limit in the number of references that can be added for a
   // single measurement
   
   if( (!xmlStrcmp(xml_cur->name, (const xmlChar *)"ref")) ) {

        XmlCitation * ref = ParseReference(xml_cur);
        header->Add(ref);

   // or parse normal values
   } else {
    for(int itag = 0; itag < c_meas_header_ntags; itag++) {
       if( (!xmlStrcmp(xml_cur->name, 
                  (const xmlChar *) c_meas_header_tag[itag].c_str())) ) {
          header->Add(c_meas_header_tag[itag], 
          utils::xml::TrimSpaces(xmlNodeListGetString(
              fXmlDoc, xml_cur->xmlChildrenNode, 1)) );
        }
    }
   }  
  
   xml_cur = xml_cur->next;  
  } //over nodes

  return header;
}
//__________________________________________________________________________
XmlCitation * NuVldXmlParser::ParseReference(xmlNodePtr xml_cur)
{
 XmlCitation * ref = new XmlCitation;

 xml_cur = xml_cur->xmlChildrenNode;
 
 while (xml_cur != NULL) {
    for(int itag = 0; itag < c_ref_ntags; itag++) {
       if( (!xmlStrcmp(xml_cur->name, 
                            (const xmlChar *) c_ref_tag[itag].c_str())) ) {

          string value = utils::xml::TrimSpaces( 
                xmlNodeListGetString(fXmlDoc, xml_cur->xmlChildrenNode, 1));

          ref->Add( c_ref_tag[itag], utils::str::FilterString(",",value) );
        }	
    } 
    xml_cur = xml_cur->next;
 }
 return ref;
}
//__________________________________________________________________________
XmlRecordBase * 
   NuVldXmlParser::ParsePoint(xmlNodePtr xml_cur, XmlObservable_t obs)
{
 //cout << "NuVldXmlParser::parse_point() : begin" << endl;

 XmlRecordBase * rec = CreateNewXmlRecord(obs);
 if(!rec) {
     SLOG("NuVld", pERROR)
        << "I can not create a record for observable: " 
        << XmlObservable::AsString(obs);
 }
 
 const vector<string> elements = rec->GetElements(); 
 vector<string>::const_iterator element_iter;
 
 string filtered("+/-");
 
 xmlChar *  xml_element = 0;
 
 xml_cur = xml_cur->xmlChildrenNode;

 for(element_iter = elements.begin(); 
    element_iter != elements.end(); ++element_iter) {

    xml_element = (xmlChar *) (*element_iter).c_str(); 

    const vector<string> attributes = rec->GetAttributes( *element_iter ); 
    vector<string>::const_iterator attribute_iter;

    xmlNodePtr xml_node = xml_cur;

    while (xml_node != NULL) {
        if( (!xmlStrcmp(xml_node->name, xml_element)) ) {
  
          // get the value associated with this element (tag)
          string value =  utils::xml::TrimSpaces(
               xmlNodeListGetString(fXmlDoc, xml_node->xmlChildrenNode, 1) ); 

          // expand / filter elements that represent errors
          if( (*element_iter).find("err") != string::npos ) {
            if( value.find("+") != string::npos ) 
                       rec->Add( *element_iter + "+", 
                       utils::str::FilterString("+/-", value));
            if( value.find("-") != string::npos ) 
                       rec->Add( *element_iter + "-",
                       utils::str::FilterString("+/-", value));	    
          } else rec->Add( *element_iter, value);
  
          // now get the values for all attributes associated with this element  
          for(attribute_iter = attributes.begin(); 
                        attribute_iter != attributes.end(); ++attribute_iter) {	  
            string attrib = utils::xml::GetAttribute(xml_node, *attribute_iter);
            rec->Add( *element_iter + "_" + *attribute_iter, attrib);
          } // element attributes
        } // if element is matched
    
        xml_node = xml_node->next;
	
    } // xml node pointers    
 } // elements

 //xmlFree(xml_element);  
 
 return rec;  
}
//__________________________________________________________________________
XmlRecordBase * NuVldXmlParser::CreateNewXmlRecord(XmlObservable_t obs)
{
 XmlRecordBase * rec = 0;
  
 switch(obs) {
   case (e_tot_xsec):
   case (e_qes_xsec):
   case (e_spp_xsec):
   case (e_mpp_xsec):
   case (e_coh_xsec):
            rec = new XmlRecord<XmlNuXSecRecord>();
            break;
   case (e_electron_diff_xsec):
            rec = new XmlRecord<XmlElDiffXSecRecord>();
            break;
   case (e_f2):
   case (e_xf3):
            rec = new XmlRecord<XmlSFRecord>();
            break;
   default:
            rec = 0;         
 }
 return rec;
}
//__________________________________________________________________________

} // nuvld namespace
} // genie namespace
