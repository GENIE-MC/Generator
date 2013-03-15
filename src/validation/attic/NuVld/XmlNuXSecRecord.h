//____________________________________________________________________________
/*!

\class    genie::nuvld::XmlNuXSecRecord

\brief    XML record structure for neutrino cross section data points                          

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  Aug, 2003

\cpright  Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

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
