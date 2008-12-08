//_____________________________________________________________________________
/*!

\class    genie::nuvld::XmlNuXSecRecord

\brief    XML record structure for neutrino cross section data points

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>

\created  August 2003
*/
//_____________________________________________________________________________

#ifndef _V_XSEC_RECORD_H_
#define _V_XSEC_RECORD_H_

#include <iostream>

#include "ValidationTools/NuVld/XmlRecordStructure.h"

using std::ostream;

namespace genie {
namespace nuvld {
  
class XmlNuXSecRecord : public XmlRecordStructure {

public:
  static XmlNuXSecRecord * Instance(void);

  friend ostream & operator<<(ostream & stream, const XmlNuXSecRecord & rec_xsec); 

private:

  XmlNuXSecRecord(); 
  XmlNuXSecRecord(const XmlNuXSecRecord & rec_xsec);
  virtual ~XmlNuXSecRecord();
  
  static XmlNuXSecRecord * _instance;

};

} // nuvld namespace
} // genie namespace

#endif // _V_XSEC_RECORD_H_
