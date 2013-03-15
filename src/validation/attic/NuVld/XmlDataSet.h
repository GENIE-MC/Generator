//____________________________________________________________________________
/*!

\class    genie::nuvld::XmlDataSet

\brief    Holds the contents of a parsed NuValidator XML data file

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  Aug, 2003

\cpright  Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _XML_DATA_SET_H_
#define _XML_DATA_SET_H_

#include <string>
#include <map>

#include "ValidationTools/NuVld/XmlExperimentMeasurements.h"

using std::string;
using std::map;

namespace genie {
namespace nuvld {
  
class XmlDataSet {

public:

  XmlDataSet();
 ~XmlDataSet();

  void Add(string str_unique_id, XmlExperimentMeasurements * e_meas);

  const map<string, XmlExperimentMeasurements *> & Get(void) const;

private:

  map<string, XmlExperimentMeasurements *> *  _data;
};

} // nuvld namespace
} // genie namespace

#endif // _DATA_SET_H_
