//_____________________________________________________________________________
/*!

\class    genie::nuvld::MeasurementIdList

\brief

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>

\created  January 2004                    
*/
//_____________________________________________________________________________

#include "DBUtils/MeasurementIdList.h"
#include "Messenger/Messenger.h"

using std::endl;
using std::cout;

using namespace genie::nuvld;

//______________________________________________________________________________
namespace genie {
 namespace nuvld {
    ostream & operator << (ostream & stream, const MeasurementIdList & list)
    {
      list.Print(stream);
      return stream;
    }
 }
}    
//______________________________________________________________________________
MeasurementIdList::MeasurementIdList()
{

}
//______________________________________________________________________________
MeasurementIdList::MeasurementIdList(const MeasurementIdList * mlist)
{
  vector<MeasurementId *>::const_iterator mid_iter;

  for( mid_iter = mlist->_id_list.begin();
                               mid_iter != mlist->_id_list.end(); ++mid_iter) {

      MeasurementId * mid = new MeasurementId( *mid_iter );

      _id_list.push_back(mid);
  }
}
//______________________________________________________________________________
MeasurementIdList::~MeasurementIdList()
{

}
//______________________________________________________________________________
unsigned int MeasurementIdList::NIds(void) const
{
  return _id_list.size();
}
//______________________________________________________________________________
void MeasurementIdList::AddId(MeasurementId * id)
{
  if(! IdExists(id) ) _id_list.push_back(id);
}
//______________________________________________________________________________
const MeasurementId * MeasurementIdList::GetId(
                                string experiment, string measurement_tag) const
{
  vector<MeasurementId *>::const_iterator list_iter;

  for(list_iter = _id_list.begin(); list_iter != _id_list.end(); ++list_iter) {

     if( experiment.compare( (*list_iter)->Experiment() ) == 0 &&
         measurement_tag.compare( (*list_iter)->MeasurementTag() ) == 0 )
                                                             return *list_iter;
  }

  return 0;
}
//______________________________________________________________________________
const MeasurementId * MeasurementIdList::GetId(unsigned int ipos) const
{
  SLOG("NuVld", pDEBUG) << "IPOS = " << ipos << " and NIDS = " << NIds();

  if(ipos < NIds()) return _id_list[ipos];
  else return 0;
}
//______________________________________________________________________________
bool MeasurementIdList::IdExists(MeasurementId * id) const
{
  vector<MeasurementId *>::const_iterator list_iter;

  for(list_iter = _id_list.begin(); list_iter != _id_list.end(); ++list_iter) {

     if( id->Experiment().compare( (*list_iter)->Experiment() ) == 0 &&
         id->MeasurementTag().compare( (*list_iter)->MeasurementTag() ) == 0 )
                                                                   return true;
  }

  return false;
}
//______________________________________________________________________________
void MeasurementIdList::Print(ostream & stream) const
{
  vector<MeasurementId *>::const_iterator list_iter;

  int id = 0;
  for(list_iter = _id_list.begin(); list_iter != _id_list.end(); ++list_iter) {

     stream << "--------- MeasurementId: " << id++ << " ------------" << endl;
     stream << *(*list_iter);
  }  
}
//______________________________________________________________________________


