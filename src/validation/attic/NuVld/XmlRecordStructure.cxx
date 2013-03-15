//____________________________________________________________________________
/*
 Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Dec 08, 2008 - CA
   NuValidator package refurbishment. Removed the neugen3 dependency.
   Moved all sources to $GENIE/src/ValidationTools/NuVld.
   Some clases have been renamed.
*/
//____________________________________________________________________________ 

#include <cstdlib>
#include <fstream>
#include <iostream>

#include "ValidationTools/NuVld/XmlRecordStructure.h"

using std::ifstream;
using std::cout;
using std::cerr;
using std::endl;
using std::ios;

namespace genie {
namespace nuvld {
  
//_____________________________________________________________________________
ostream & operator << (ostream & stream, const XmlRecordStructure & recbase)
{
  map<string, vector<string> >::const_iterator element_iter;
  vector<string>::const_iterator               attribute_iter;

  stream << "[record structure]" << endl;
  for(element_iter = recbase._elements.begin(); 
                    element_iter != recbase._elements.end(); ++element_iter) {

      stream << "[-] element: " << (*element_iter).first << endl;;
      vector<string> attributes = (*element_iter).second;
      for(attribute_iter = attributes.begin(); 
                       attribute_iter != attributes.end(); ++attribute_iter) {

         stream << " |----> with attribute: " << (*attribute_iter) << endl;
      }
  }
  return stream;
}
//_____________________________________________________________________________
XmlRecordStructure::XmlRecordStructure()
{

}
//_____________________________________________________________________________
XmlRecordStructure::XmlRecordStructure(const char * filename)
{
  ReadElements(filename);
}
//_____________________________________________________________________________
XmlRecordStructure::~XmlRecordStructure()
{

}
//_____________________________________________________________________________
int XmlRecordStructure::ReadElements(const char * filename)
{
  ifstream inp(filename,ios::in);

  if( inp.fail() ) {
    cerr << "Failed to open filename: " << filename << endl;
    exit(1);
  }

  string         name;
  vector<string> names;
  while ( inp >> name ) names.push_back(name);

  vector<string>::reverse_iterator name_iter;
  vector<string> attributes;

  for(name_iter = names.rbegin(); name_iter != names.rend(); ++name_iter) {
     if( (*name_iter)[0] == '+') {
         string element = (*name_iter).erase(0,1);
         _elements.insert( 
               map<string, vector<string> >::value_type(element, attributes));
         attributes.clear();                    
     } else {
        attributes.push_back( (*name_iter).erase(0,1) );
     }
  }

  inp.close();

  return 0;
}
//_____________________________________________________________________________

} // nuvld namespace
} // genie namespace
