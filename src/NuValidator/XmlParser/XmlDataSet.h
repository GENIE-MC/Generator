//_____________________________________________________________________________
/*!

\class    genie::nuvld::XmlDataSet

\brief    Holds the contents of a parsed NuValidator XML data file

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>

\created  August 2003
*/
//_____________________________________________________________________________

#ifndef _XML_DATA_SET_H_
#define _XML_DATA_SET_H_

#include <string>
#include <map>

#include "ExperimentMeasurements.h"

using std::string;
using std::map;

namespace genie {
namespace nuvld {
  
class XmlDataSet {

public:

  XmlDataSet();
 ~XmlDataSet();

  void Add(string str_unique_id, ExperimentMeasurements * e_meas);

  const map<string, ExperimentMeasurements *> & Get(void) const;

private:

  map<string, ExperimentMeasurements *> *  _data;
};

} // nuvld namespace
} // genie namespace

#endif // _DATA_SET_H_
