//_____________________________________________________________________________
/*!

\class    genie::nuvld::RecordBase

\brief

\author   Costas Andreopoulos (Rutherford Lab.)  <costas.andreopoulos \at stfc.ac.uk>

\created  August 2003         
*/
//_____________________________________________________________________________

#include "RecordBase.h"

using std::endl;

namespace genie {
namespace nuvld {
  
//_____________________________________________________________________________
ostream & operator << (ostream & stream, const RecordBase & rec_base)
{
  map<string, string>::const_iterator rec_iter;

  stream << "----------------------------------------printing record:" << endl;
  for(rec_iter = rec_base._rec.begin(); 
                      rec_iter != rec_base._rec.end(); ++rec_iter) 
                                    stream << (*rec_iter).first << " :....... " 
                                                 << (*rec_iter).second << endl;
  return stream;
}
//_____________________________________________________________________________
RecordBase::RecordBase()
{

}
//_____________________________________________________________________________
RecordBase::~RecordBase()
{

}
//_____________________________________________________________________________

} // nuvld namespace
} // genie namespace
