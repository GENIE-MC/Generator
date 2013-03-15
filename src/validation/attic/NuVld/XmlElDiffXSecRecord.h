//____________________________________________________________________________
/*!

\class    genie::nuvld::XmlElDiffXSecRecord

\brief

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  Aug, 2003

\cpright  Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

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
