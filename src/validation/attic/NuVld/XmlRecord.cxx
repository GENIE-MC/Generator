//____________________________________________________________________________
/*
 Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory - Aug 01, 2003

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Dec 08, 2008 - CA
   NuValidator package refurbishment. Removed the neugen3 dependency.
   Moved all sources to $GENIE/src/ValidationTools/NuVld.
   Some clases have been renamed.
*/
//____________________________________________________________________________ 

#include <typeinfo>

#include "ValidationTools/NuVld/XmlRecord.h"
#include "ValidationTools/NuVld/XmlNuXSecRecord.h"
#include "ValidationTools/NuVld/XmlElDiffXSecRecord.h"
#include "ValidationTools/NuVld/XmlSFRecord.h"

using std::cout;
using std::endl;

namespace genie {
namespace nuvld {
  
template ostream & operator 
                   << (ostream & stream, const XmlRecord<XmlNuXSecRecord> & rec);
template ostream & operator 
               << (ostream & stream, const XmlRecord<XmlElDiffXSecRecord> & rec);
template ostream & operator 
                      << (ostream & stream, const XmlRecord<XmlSFRecord> & rec);

//____________________________________________________________________________
template <class T>
ostream & operator << (ostream & stream, const XmlRecord<T> & rec)
{
  try {
    const XmlRecordBase & rec_base = dynamic_cast<const XmlRecordBase &> (rec);
    stream << rec_base;
  }
  catch( std::bad_cast ) {
   stream
        << "warning: could not dynamic_cast to a reference of the base object"
        << endl;
  }
  return stream;
}
//____________________________________________________________________________
template <class T> XmlRecord<T>::XmlRecord() : XmlRecordBase()
{
  _structure = T::Instance();
}
//____________________________________________________________________________
template <class T> XmlRecord<T>::~XmlRecord()
{

}
//____________________________________________________________________________
template <class T> void XmlRecord<T>::Add(string key, string value)
{
  _rec.insert(map<string, string>::value_type(key, value));
}
//____________________________________________________________________________
template <class T> const string XmlRecord<T>::Get(string key) const
{
  map<string, string>::const_iterator iter = _rec.find(key);

  if(iter != _rec.end() ) return (*iter).second;
  else                    return "0";
}
//____________________________________________________________________________
template <class T> void XmlRecord<T>::PrintStructure(void) const 
{ 
  cout << *_structure << endl; 
}
//____________________________________________________________________________
template <class T> const vector<string> XmlRecord<T>::GetElements(void) const
{
  vector<string> elements;

  const map<string, vector<string> > struct_map = _structure->get_structure();
  map<string, vector<string> >::const_iterator struct_iter;

  for(struct_iter = struct_map.begin(); 
                             struct_iter != struct_map.end(); ++struct_iter) {
     elements.push_back( (*struct_iter).first );
  }

  return elements;
}
//____________________________________________________________________________
template <class T> 
           const vector<string> XmlRecord<T>::GetAttributes(string element) const
{
  vector<string> attrib;

  const map<string, vector<string> > struct_map = _structure->get_structure();
  map<string, vector<string> >::const_iterator struct_iter 
                                                   = struct_map.find(element);
  attrib = (*struct_iter).second;

  return attrib;
}
//____________________________________________________________________________

// template specializations:
//
template class XmlRecord<XmlNuXSecRecord>;
template class XmlRecord<XmlElDiffXSecRecord>;
template class XmlRecord<XmlSFRecord>;

} // nuvld namespace
} // genie namespace
