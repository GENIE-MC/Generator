//_____________________________________________________________________________
/*!

\class    genie::nuvld::NuVldXmlParser

\brief    A libxml2 based parser for the NuValidator XML data files

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>

\created  August 2003          
*/
//_____________________________________________________________________________

#include "Messenger/Messenger.h"
#include "XmlParser/NuVldXmlParser.h" 
#include "XmlParser/ParserUtils.h" 
#include "XmlParser/vXSecRecord.h"
#include "XmlParser/eDiffXSecRecord.h"

namespace genie {
namespace nuvld {
  
//__________________________________________________________________________
NuVldXmlParser::NuVldXmlParser()
{
  _data = new XmlDataSet;
} 
//__________________________________________________________________________
NuVldXmlParser::~NuVldXmlParser()
{

}
//__________________________________________________________________________
void NuVldXmlParser::ParseXmlDocument(const char * filename)
{
  _xml_filename = string(filename);   
  _xml_doc      = xmlParseFile( _xml_filename.c_str() );

  if( (_parser_status = VerifyParsing()) == eXml_OK ) FillDataSet();
}
//__________________________________________________________________________
const XmlDataSet & NuVldXmlParser::GetDataSet(void) const
{
  return *_data;
}
//__________________________________________________________________________
XmlParserStatus_t NuVldXmlParser::VerifyParsing(void)
{
 SLOG("NuVld", pINFO)  << "NuVldXmlParser::VerifyParsing() :......: ";

 XmlParserStatus_t status;
 
 if(_xml_doc==NULL) { 
   status = eXml_NOT_PARSED;

   SLOG("NuVld", pERROR)  << "***" << ParserStatus::AsString(status);
   return status;
 }
 
 xmlNodePtr xml_cur = xmlDocGetRootElement(_xml_doc);

 if(xml_cur==NULL) {
   status = eXml_EMPTY;

   SLOG("NuVld", pERROR)  << "***" << ParserStatus::AsString(status);
   return status;
 }

 if( xmlStrcmp(xml_cur->name, (const xmlChar *) "nuscat_data") ) {

   status = eXml_INVALID_ROOT;

   xmlFreeDoc(_xml_doc);
   SLOG("NuVld", pERROR)  << "***" << ParserStatus::AsString(status);
   return status;
 }

 status = eXml_OK;
 
 SLOG("NuVld", pINFO) << ParserStatus::AsString(status);
 
 return status;
}
//__________________________________________________________________________
void NuVldXmlParser::FillDataSet(void)
{
 SLOG("NuVld", pDEBUG) << "Starts parsing the XML file";
 
 xmlChar *  xml_exp_name = 0;

 xmlNodePtr xml_cur = xmlDocGetRootElement(_xml_doc);
 
 xml_cur = xml_cur->xmlChildrenNode;

 // loop over all xml tree nodes that are children of the root node
 while (xml_cur != NULL) {

   // enter everytime you find an 'experiment' tag
   if( (!xmlStrcmp(xml_cur->name, (const xmlChar *) "experiment")) ) {

     // get experiment attributes (name,tag) - the last can be absent
     // if there is only one experiment entry with the given name
     
     xml_exp_name     = xmlGetProp(xml_cur, (const xmlChar *)"name");

     string name  = ParserUtils::trim_spaces( xml_exp_name );

     xmlFree(xml_exp_name);     

     SLOG("NuVld", pDEBUG) << "Parsing entries for experiment: " << name;

     // parse all measurements associated with the given experiment
     
     ExperimentMeasurements * meas_list = ParseExperiment(xml_cur, name);

     // add the entry to the top-level STL map
     _data->Add( name, meas_list );
   }
   
   xml_cur = xml_cur->next;
 }

 SLOG("NuVld", pDEBUG) << "Finished parsing XML data file";
}
//__________________________________________________________________________
ExperimentMeasurements * NuVldXmlParser::ParseExperiment(
                                            xmlNodePtr xml_cur, string name)
{
 ExperimentMeasurements * meas_list = new ExperimentMeasurements;
 
 _measurement_tag = 0; 
 
 xml_cur = xml_cur->xmlChildrenNode;

 while (xml_cur != NULL) {

    //------------------ parse 'experiment-info' section ------------------

    if( (!xmlStrcmp(xml_cur->name, (const xmlChar *) "exp_info")) ) {

        SLOG("NuVld", pDEBUG) << "Parsing experiment info";
       
       ExperimentInfo * e_info = ParseExperimentInfo(xml_cur);
       
       e_info->Add("name", name);
       
       cout << *e_info << endl; // send the experiment info to stdout
       
       meas_list->Add( e_info );
    }

    //-------------------- parse 'beam-spectrum' section ------------------

    if( (!xmlStrcmp(xml_cur->name, (const xmlChar *) "flux_spectrum")) ) {
       
       string beam    = ParserUtils::get_attribute(xml_cur, "beam");       
       string F_units = ParserUtils::get_attribute(xml_cur, "flux_units");
       string E_units = ParserUtils::get_attribute(xml_cur, "E_units");
       string E_frame = ParserUtils::get_attribute(xml_cur, "E_frame");

       SLOG("NuVld", pDEBUG) << "Parsing flux spectrum for beam: " << beam;

       BeamFluxSpectrum * flux = ParseBeamFluxSpectrum(xml_cur);
       
       flux->SetFluxUnits(F_units);
       flux->SetEnergyUnits(E_units);
       flux->SetEnergyFrame(E_frame);

       cout << *flux << endl;
       
       meas_list->Add( beam, flux );
    }

    //--------------------- parse 'measurement' section -------------------

    if( (!xmlStrcmp(xml_cur->name, (const xmlChar *) "measurement")) ) {

       //-------- understand what kind of measurement this is...
       string obs = ParserUtils::get_attribute(xml_cur, "observable");
              
       Observable_t obst = Observable::GetObservable(obs);
 
       Measurement * meas = ParseMeasurement(xml_cur, obst);
       
       meas_list->Add( meas );
    }

    xml_cur = xml_cur->next; // move to next xml node
    
 } // loop over xml nodes
 
 return meas_list;
}
//__________________________________________________________________________
ExperimentInfo * NuVldXmlParser::ParseExperimentInfo(xmlNodePtr xml_cur)
{
 xmlChar *  xml_string;
 string     std_string;

 ExperimentInfo * e_info = new ExperimentInfo;
 
 string comment = ParserUtils::trim_spaces( 
             xmlNodeListGetString(_xml_doc, xml_cur->xmlChildrenNode, 1) );

 e_info->Add("comment", ParserUtils::filter_string(",",comment));

 xml_cur = xml_cur->xmlChildrenNode;
  
 while (xml_cur != NULL) {

    for(int itag = 0; itag < c_exp_info_ntags; itag++) {

       // the element names are defined in ExperimentInfo.h       
       if( (!xmlStrcmp(xml_cur->name, 
                       (const xmlChar *) c_exp_info_tag[itag].c_str())) ) {

          xml_string = xmlNodeListGetString(_xml_doc, xml_cur->xmlChildrenNode, 1); 
          std_string = ParserUtils::trim_spaces(xml_string);

          SLOG("NuVld", pINFO)
                << "Adding " << c_exp_info_tag[itag] << ": " << std_string;
                
          e_info->Add(c_exp_info_tag[itag], ParserUtils::filter_string(",",std_string)); 

          xmlFree(xml_string);     
        }
    }
    xml_cur = xml_cur->next;
 } // over xml nodes
 
 return e_info;
}
//__________________________________________________________________________
BeamFluxSpectrum * NuVldXmlParser::ParseBeamFluxSpectrum(xmlNodePtr xml_cur)
{
  BeamFluxSpectrum * spectrum = new BeamFluxSpectrum;

  xml_cur = xml_cur->xmlChildrenNode;

  while (xml_cur != NULL) {

    if( (!xmlStrcmp(xml_cur->name, (const xmlChar *) "bin")) ) {
    
        BeamFluxBin * spectrum_bin = ParseBeamFluxBin(xml_cur);
        spectrum->Add( spectrum_bin );
    }

    xml_cur = xml_cur->next;
  }  
  return spectrum;
}
//__________________________________________________________________________
BeamFluxBin * NuVldXmlParser::ParseBeamFluxBin(xmlNodePtr xml_cur)
{
  xml_cur = xml_cur->xmlChildrenNode;

  xmlChar *  xml_string;
  string     std_string;

  vector<string> energies;
  vector<string> fluxes;
  
  while (xml_cur != NULL) {

    if( (!xmlStrcmp(xml_cur->name, (const xmlChar *) "E")) ) {

        xml_string = xmlNodeListGetString(_xml_doc, xml_cur->xmlChildrenNode, 1); 
        std_string = ParserUtils::trim_spaces(xml_string);

        energies = ParserUtils::split( std_string, string(",") );

    }
    
    if( (!xmlStrcmp(xml_cur->name, (const xmlChar *) "Flux")) ) {

        xml_string = xmlNodeListGetString(_xml_doc, xml_cur->xmlChildrenNode, 1); 
        std_string = ParserUtils::trim_spaces(xml_string);

        fluxes = ParserUtils::split( std_string, string(",") );
    }
    
    xml_cur = xml_cur->next;
  }  

  BeamFluxBin * spectrum_bin = new BeamFluxBin(
                                    energies[0], energies[1], energies[2], 
                                    fluxes[0],   fluxes[1],   fluxes[2]   );

  return spectrum_bin;
}
//__________________________________________________________________________
Measurement * NuVldXmlParser::ParseMeasurement(
                                      xmlNodePtr xml_cur, Observable_t obs)
{
  SLOG("NuVld", pINFO) << "Start parsing a <<" 
                           << Observable::AsString(obs) << ">> measurement";

  xml_cur = xml_cur->xmlChildrenNode;

  Measurement *       meas   = new Measurement;
  MeasurementHeader * header = 0;

  while (xml_cur != NULL) {

    // parse the measurement header

    if( (!xmlStrcmp(xml_cur->name, (const xmlChar *) "header")) ) {
               
        header = ParseMeasurementHeader(xml_cur);

        header->Add("observable", Observable::Code(obs));
        header->Add("tag",        ParserUtils::int_as_string(_measurement_tag));

        meas->Add( header ); 
        cout << *header << endl;
    }

    // parse every occurence of 'point' within the current measurement
    
    if( (!xmlStrcmp(xml_cur->name, (const xmlChar *) "point")) ) {
               
        RecordBase * rec = ParsePoint(xml_cur, obs);  

        cout << *rec;

        meas->Add(rec);
    }

    xml_cur = xml_cur->next;
  }

  _measurement_tag++;
  
  return meas;
}
//__________________________________________________________________________
MeasurementHeader * NuVldXmlParser::ParseMeasurementHeader(xmlNodePtr xml_cur)
{
  MeasurementHeader * header = new MeasurementHeader;

  string comment = ParserUtils::trim_spaces( 
             xmlNodeListGetString(_xml_doc, xml_cur->xmlChildrenNode, 1) );

  header->Add("comment", ParserUtils::filter_string(",",comment));

  xml_cur = xml_cur->xmlChildrenNode;

  while (xml_cur != NULL) {

   // Add a reference that is listed in relation with this measurement.
   // There is no limit in the number of references that can be added for a
   // single measurement
   
   if( (!xmlStrcmp(xml_cur->name, (const xmlChar *)"ref")) ) {

        Citation * ref = ParseReference(xml_cur);
        header->Add(ref);

   // or parse normal values
   } else {

    for(int itag = 0; itag < c_meas_header_ntags; itag++) {

       if( (!xmlStrcmp(xml_cur->name, 
                  (const xmlChar *) c_meas_header_tag[itag].c_str())) ) {

          header->Add(c_meas_header_tag[itag], 
          ParserUtils::trim_spaces(xmlNodeListGetString(
                 _xml_doc, xml_cur->xmlChildrenNode, 1)) );
        }
    }
   }  
  
   xml_cur = xml_cur->next;  
  } //over nodes

  return header;
}
//__________________________________________________________________________
Citation * NuVldXmlParser::ParseReference(xmlNodePtr xml_cur)
{
 Citation * ref = new Citation;

 xml_cur = xml_cur->xmlChildrenNode;
 
 while (xml_cur != NULL) {

    for(int itag = 0; itag < c_ref_ntags; itag++) {

       if( (!xmlStrcmp(xml_cur->name, 
                            (const xmlChar *) c_ref_tag[itag].c_str())) ) {


        string value = ParserUtils::trim_spaces( xmlNodeListGetString(
                                    _xml_doc, xml_cur->xmlChildrenNode, 1));

          ref->Add( c_ref_tag[itag], ParserUtils::filter_string(",",value) );
        }	
    } 
    xml_cur = xml_cur->next;
 }
 return ref;
}
//__________________________________________________________________________
RecordBase * NuVldXmlParser::ParsePoint(xmlNodePtr xml_cur, Observable_t obs)
{
 //cout << "NuVldXmlParser::parse_point() : begin" << endl;

 RecordBase * rec = CreateNewRecord(obs);

 if(!rec) {
     SLOG("NuVld", pERROR)
             << "I can not create a record for observable: " 
                                              << Observable::AsString(obs);
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

          string value =  ParserUtils::trim_spaces(
               xmlNodeListGetString(_xml_doc, xml_node->xmlChildrenNode, 1) ); 

          // expand / filter elements that represent errors

          if( (*element_iter).find("err") != string::npos ) {
            
            if( value.find("+") != string::npos ) 
                       rec->Add( *element_iter + "+", 
                       ParserUtils::filter_string("+/-", value));
            if( value.find("-") != string::npos ) 
                       rec->Add( *element_iter + "-",
                       ParserUtils::filter_string("+/-", value));	    
          } else rec->Add( *element_iter, value);
  
          // now get the values for all attributes associated with this element
  
          for(attribute_iter = attributes.begin(); 
                        attribute_iter != attributes.end(); ++attribute_iter) {	  

                  string attrib = ParserUtils::get_attribute(
                                              xml_node, *attribute_iter);

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
RecordBase * NuVldXmlParser::CreateNewRecord(Observable_t obs)
{
 RecordBase * rec = 0;
  
 switch(obs) {
   
   case (e_tot_xsec):
   case (e_qes_xsec):
   case (e_spp_xsec):
   case (e_mpp_xsec):
   case (e_coh_xsec):

            rec = new Record<vXSecRecord>();
            break;
                               
   case (e_electron_diff_xsec):

            rec = new Record<eDiffXSecRecord>();
            break;

   default:
            rec = 0;         
 }
 return rec;
}
//__________________________________________________________________________

} // nuvld namespace
} // genie namespace
