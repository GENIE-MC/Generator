//_____________________________________________________________________________
/*!

\class    genie::nuvld::MeasurementIdList

\brief

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>

\created  January 2004
*/
//_____________________________________________________________________________

#ifndef _MEASUREMENT_ID_LIST_H_
#define _MEASUREMENT_ID_LIST_H_

#include <string>
#include <vector>
#include <iostream>

#include "MeasurementId.h"

using std::string;
using std::vector;
using std::ostream;

namespace genie {
namespace nuvld {

class MeasurementIdList 
{
public:

  friend ostream & operator << (ostream & stream, const MeasurementIdList & list);

  MeasurementIdList();
  MeasurementIdList(const MeasurementIdList * mlist);
  ~MeasurementIdList();

  void  AddId     (MeasurementId * id);
  bool  IdExists  (MeasurementId * id) const;

  const MeasurementId * GetId(string experiment, string measurement_tag) const;
  const MeasurementId * GetId(unsigned int ipos) const;

  unsigned int NIds  (void)              const;
  void         Print  (ostream & stream)  const;

private:

  vector<MeasurementId *> _id_list;
};

} // nuvld namespace
} // genie namespace

#endif
