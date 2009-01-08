//_____________________________________________________________________________
/*!

\class    genie::nuvld::vXSecRecord

\brief    XML record structure for neutrino cross section data points

\author   Costas Andreopoulos (Rutherford Lab.)  <costas.andreopoulos \at stfc.ac.uk>

\created  August 2003          
*/
//_____________________________________________________________________________

#include <cstdlib>
#include <string>

#include "vXSecRecord.h"

using std::string;
using std::endl;
using std::cerr;
using std::cout;

namespace genie {
namespace nuvld {
  
//____________________________________________________________________________
ostream & operator<<(ostream & stream, const vXSecRecord & rec_xsec) 
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
vXSecRecord * vXSecRecord::_instance = 0; 
//____________________________________________________________________________
vXSecRecord::vXSecRecord() : RecordStructure()
{
  string path = string( getenv("GENIE") ) + 
                  string("/src/NuValidator/XmlParser/record_v_xsec.elements");

  ReadElements(path.c_str());
  
  _instance = 0;
}
//____________________________________________________________________________
vXSecRecord::~vXSecRecord()
{
  _instance = 0;
} 
//____________________________________________________________________________
vXSecRecord * vXSecRecord::Instance()
{
  if(_instance == 0) _instance = new vXSecRecord();
  return _instance;
}
//____________________________________________________________________________

} // nuvld namespace
} // genie namespace
