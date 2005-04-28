//_____________________________________________________________________________
/*!

\class    genie::nuvld::SFRecord

\brief    XML record structure for structure function data points

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>

\created  August 2003
*/
//_____________________________________________________________________________

#ifndef _SF_RECORD_H_
#define _SF_RECORD_H_

#include <iostream>

#include "RecordStructure.h"

using std::ostream;

namespace genie {
namespace nuvld {
  
class SFRecord : public RecordStructure {

public:
  static SFRecord * Instance(void);

  friend ostream & operator<<(ostream & stream, const SFRecord & rec_xsec); 

private:

  SFRecord(); 
  SFRecord(const SFRecord & rec_xsec);
  virtual ~SFRecord();
  
  static SFRecord * _instance;

};

} // nuvld namespace
} // genie namespace

#endif // _SF_RECORD_H_
