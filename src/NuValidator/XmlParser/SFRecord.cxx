//_____________________________________________________________________________
/*!

\class    genie::nuvld::SFRecord

\brief    XML record structure for structure function data points

\author   Costas Andreopoulos (Rutherford Lab.)  <costas.andreopoulos \at stfc.ac.uk>

\created  August 2003          
*/
//_____________________________________________________________________________

#include <cstdlib>
#include <string>

#include "SFRecord.h"

using std::string;
using std::endl;
using std::cerr;
using std::cout;

namespace genie {
namespace nuvld {
  
//____________________________________________________________________________
ostream & operator<<(ostream & stream, const SFRecord & rec_xsec) 
{  
  try {
   const RecordStructure & rec_str = 
                            dynamic_cast<const RecordStructure &> (rec_xsec);
   stream << rec_str;
  } 
  catch( std::bad_cast ) {
   stream 
        << "warning: could not dynamic_cast to a reference of the base object"
        << endl;
  }
   return stream;
}
//____________________________________________________________________________
SFRecord * SFRecord::_instance = 0; 
//____________________________________________________________________________
SFRecord::SFRecord() : RecordStructure()
{
  string path = string( getenv("GENIE") ) + 
                      string("/src/NuValidator/XmlParser/record_sf.elements");

  ReadElements(path.c_str());
  
  _instance = 0;
}
//____________________________________________________________________________
SFRecord::~SFRecord()
{
  _instance = 0;
} 
//____________________________________________________________________________
SFRecord * SFRecord::Instance()
{
  if(_instance == 0) _instance = new SFRecord();
  return _instance;
}
//____________________________________________________________________________

} // nuvld namespace
} // genie namespace
