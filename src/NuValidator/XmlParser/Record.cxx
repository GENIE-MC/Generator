//_____________________________________________________________________________
/*!

\class    Record

\brief

\author   Costas Andreopoulos (Rutherford Lab.)  <costas.andreopoulos \at stfc.ac.uk>

\created  August 2003         
*/
//_____________________________________________________________________________

#include "Record.h"

using std::cout;
using std::endl;

#include "vXSecRecord.h"
#include "eDiffXSecRecord.h"
#include "SFRecord.h"

namespace genie {
namespace nuvld {
  
template ostream & operator 
                   << (ostream & stream, const Record<vXSecRecord> & rec);
template ostream & operator 
               << (ostream & stream, const Record<eDiffXSecRecord> & rec);
template ostream & operator 
                      << (ostream & stream, const Record<SFRecord> & rec);

//____________________________________________________________________________
template <class T>
ostream & operator << (ostream & stream, const Record<T> & rec)
{
  try {
    const RecordBase & rec_base = dynamic_cast<const RecordBase &> (rec);
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
template <class T> Record<T>::Record() : RecordBase()
{
  _structure = T::Instance();
}
//____________________________________________________________________________
template <class T> Record<T>::~Record()
{

}
//____________________________________________________________________________
template <class T> void Record<T>::Add(string key, string value)
{
  _rec.insert(map<string, string>::value_type(key, value));
}
//____________________________________________________________________________
template <class T> const string Record<T>::Get(string key) const
{
  map<string, string>::const_iterator iter = _rec.find(key);

  if(iter != _rec.end() ) return (*iter).second;
  else                    return "0";
}
//____________________________________________________________________________
template <class T> void Record<T>::PrintStructure(void) const 
{ 
  cout << *_structure << endl; 
}
//____________________________________________________________________________
template <class T> const vector<string> Record<T>::GetElements(void) const
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
           const vector<string> Record<T>::GetAttributes(string element) const
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
template class Record<vXSecRecord>;
template class Record<eDiffXSecRecord>;
template class Record<SFRecord>;

} // nuvld namespace
} // genie namespace
