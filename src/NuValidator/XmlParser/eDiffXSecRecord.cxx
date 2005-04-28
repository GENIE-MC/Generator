//_____________________________________________________________________________
/*!

\class    genie::nuvld::eDiffXSecRecord

\brief    XML record structure for eN differential cross section data points

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>

\created  August 2003          
*/
//_____________________________________________________________________________

#include <cstdlib>
#include <string>

#include "eDiffXSecRecord.h"

using std::string;
using std::endl;
using std::cerr;
using std::cout;

namespace genie {
namespace nuvld {
  
//____________________________________________________________________________
ostream & operator<<(ostream & stream, const eDiffXSecRecord & rec_xsec) 
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
eDiffXSecRecord * eDiffXSecRecord::_instance = 0; 
//____________________________________________________________________________
eDiffXSecRecord::eDiffXSecRecord() : RecordStructure()
{
  string path = string( getenv("GENIE") ) + 
              string("/src/NuValidator/XmlParser/record_e_diff_xsec.elements");

  ReadElements(path.c_str());
  
  _instance = 0;
}
//____________________________________________________________________________
eDiffXSecRecord::~eDiffXSecRecord()
{
  _instance = 0;
} 
//____________________________________________________________________________
eDiffXSecRecord * eDiffXSecRecord::Instance()
{
  if(_instance == 0) _instance = new eDiffXSecRecord();
  return _instance;
}
//____________________________________________________________________________

} // nuvld namespace
} // genie namespace
