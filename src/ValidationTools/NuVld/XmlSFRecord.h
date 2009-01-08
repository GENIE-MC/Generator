//_____________________________________________________________________________
/*!

\class    genie::nuvld::XmlSFRecord

\brief    XML record structure for structure function data points

\author   Costas Andreopoulos (Rutherford Lab.)  <costas.andreopoulos \at stfc.ac.uk>

\created  August 2003
*/
//_____________________________________________________________________________

#ifndef _SF_RECORD_H_
#define _SF_RECORD_H_

#include <iostream>

#include "ValidationTools/NuVld/XmlRecordStructure.h"

using std::ostream;

namespace genie {
namespace nuvld {
  
class XmlSFRecord : public XmlRecordStructure {

public:
  static XmlSFRecord * Instance(void);

  friend ostream & operator<<(ostream & stream, const XmlSFRecord & rec_xsec); 

private:

  XmlSFRecord(); 
  XmlSFRecord(const XmlSFRecord & rec_xsec);
  virtual ~XmlSFRecord();
  
  static XmlSFRecord * _instance;

};

} // nuvld namespace
} // genie namespace

#endif // _SF_RECORD_H_
