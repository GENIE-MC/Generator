//____________________________________________________________________________
/*!

\class    genie::nuvld::XmlSFRecord

\brief    XML record structure for structure function data points

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  Aug, 2003

\cpright  Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

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
