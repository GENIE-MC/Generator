//_____________________________________________________________________________
/*!

\class    genie::nuvld::vXSecRecord

\brief    XML record structure for neutrino cross section data points

\author   Costas Andreopoulos (Rutherford Lab.)  <costas.andreopoulos \at stfc.ac.uk>

\created  August 2003
*/
//_____________________________________________________________________________

#ifndef _V_XSEC_RECORD_H_
#define _V_XSEC_RECORD_H_

#include <iostream>

#include "RecordStructure.h"

using std::ostream;

namespace genie {
namespace nuvld {
  
class vXSecRecord : public RecordStructure {

public:
  static vXSecRecord * Instance(void);

  friend ostream & operator<<(ostream & stream, const vXSecRecord & rec_xsec); 

private:

  vXSecRecord(); 
  vXSecRecord(const vXSecRecord & rec_xsec);
  virtual ~vXSecRecord();
  
  static vXSecRecord * _instance;

};

} // nuvld namespace
} // genie namespace

#endif // _V_XSEC_RECORD_H_
