//_____________________________________________________________________________
/*!

\class    genie::nuvld::XmlElDiffXSecRecord

\brief    XML record structure for eN differential cross section data points

\author   Costas Andreopoulos (Rutherford Lab.)  <costas.andreopoulos \at stfc.ac.uk>

\created  August 2003
*/
//_____________________________________________________________________________

#ifndef _E_DIFF_XSEC_RECORD_H_
#define _E_DIFF_XSEC_RECORD_H_

#include <iostream>

#include "ValidationTools/NuVld/XmlRecordStructure.h"

using std::ostream;

namespace genie {
namespace nuvld {
  
class XmlElDiffXSecRecord : public XmlRecordStructure {

public:
  static XmlElDiffXSecRecord * Instance(void);

  friend ostream & operator<<(ostream & stream, const XmlElDiffXSecRecord & rec_xsec); 

private:

  XmlElDiffXSecRecord(); 
  XmlElDiffXSecRecord(const XmlElDiffXSecRecord & rec_xsec);
  virtual ~XmlElDiffXSecRecord();
  
  static XmlElDiffXSecRecord * _instance;
};

} // nuvld namespace
} // genie namespace

#endif // _E_DIFF_XSEC_RECORD_H_
