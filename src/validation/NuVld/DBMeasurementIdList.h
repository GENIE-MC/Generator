//_____________________________________________________________________________
/*!

\class    genie::nuvld::DBMeasurementIdList

\brief

\author   Costas Andreopoulos (Rutherford Lab.)  <costas.andreopoulos \at stfc.ac.uk>

\created  January 2004
*/
//_____________________________________________________________________________

#ifndef _MEASUREMENT_ID_LIST_H_
#define _MEASUREMENT_ID_LIST_H_

#include <string>
#include <vector>
#include <iostream>

#include "DBMeasurementId.h"

using std::string;
using std::vector;
using std::ostream;

namespace genie {
namespace nuvld {

class DBMeasurementIdList 
{
public:

  friend ostream & operator << (ostream & stream, const DBMeasurementIdList & list);

  DBMeasurementIdList();
  DBMeasurementIdList(const DBMeasurementIdList * mlist);
  ~DBMeasurementIdList();

  void  AddId     (DBMeasurementId * id);
  bool  IdExists  (DBMeasurementId * id) const;

  const DBMeasurementId * GetId(string experiment, string measurement_tag) const;
  const DBMeasurementId * GetId(unsigned int ipos) const;

  unsigned int NIds  (void)              const;
  void         Print  (ostream & stream)  const;

private:

  vector<DBMeasurementId *> _id_list;
};

} // nuvld namespace
} // genie namespace

#endif
