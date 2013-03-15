//____________________________________________________________________________
/*
 Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory - Jan 12, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Dec 08, 2008 - CA
   NuValidator package refurbishment. Removed the neugen3 dependency. 
   Moved all sources to $GENIE/src/ValidationTools/NuVld.
   Some clases have been renamed.
*/
//____________________________________________________________________________ 

#include "ValidationTools/NuVld/DBMeasurementIdList.h"
#include "Messenger/Messenger.h"

using std::endl;
using std::cout;

using namespace genie::nuvld;

//______________________________________________________________________________
namespace genie {
 namespace nuvld {
    ostream & operator << (ostream & stream, const DBMeasurementIdList & list)
    {
      list.Print(stream);
      return stream;
    }
 }
}    
//______________________________________________________________________________
DBMeasurementIdList::DBMeasurementIdList()
{

}
//______________________________________________________________________________
DBMeasurementIdList::DBMeasurementIdList(const DBMeasurementIdList * mlist)
{
  vector<DBMeasurementId *>::const_iterator mid_iter;

  for( mid_iter = mlist->_id_list.begin();
                               mid_iter != mlist->_id_list.end(); ++mid_iter) {

      DBMeasurementId * mid = new DBMeasurementId( *mid_iter );

      _id_list.push_back(mid);
  }
}
//______________________________________________________________________________
DBMeasurementIdList::~DBMeasurementIdList()
{

}
//______________________________________________________________________________
unsigned int DBMeasurementIdList::NIds(void) const
{
  return _id_list.size();
}
//______________________________________________________________________________
void DBMeasurementIdList::AddId(DBMeasurementId * id)
{
  if(! IdExists(id) ) _id_list.push_back(id);
}
//______________________________________________________________________________
const DBMeasurementId * DBMeasurementIdList::GetId(
                                string experiment, string measurement_tag) const
{
  vector<DBMeasurementId *>::const_iterator list_iter;

  for(list_iter = _id_list.begin(); list_iter != _id_list.end(); ++list_iter) {

     if( experiment.compare( (*list_iter)->Experiment() ) == 0 &&
         measurement_tag.compare( (*list_iter)->XmlMeasurementTag() ) == 0 )
                                                             return *list_iter;
  }

  return 0;
}
//______________________________________________________________________________
const DBMeasurementId * DBMeasurementIdList::GetId(unsigned int ipos) const
{
  SLOG("NuVld", pDEBUG) << "IPOS = " << ipos << " and NIDS = " << NIds();

  if(ipos < NIds()) return _id_list[ipos];
  else return 0;
}
//______________________________________________________________________________
bool DBMeasurementIdList::IdExists(DBMeasurementId * id) const
{
  vector<DBMeasurementId *>::const_iterator list_iter;

  for(list_iter = _id_list.begin(); list_iter != _id_list.end(); ++list_iter) {

     if( id->Experiment().compare( (*list_iter)->Experiment() ) == 0 &&
         id->XmlMeasurementTag().compare( (*list_iter)->XmlMeasurementTag() ) == 0 )
                                                                   return true;
  }

  return false;
}
//______________________________________________________________________________
void DBMeasurementIdList::Print(ostream & stream) const
{
  vector<DBMeasurementId *>::const_iterator list_iter;

  int id = 0;
  for(list_iter = _id_list.begin(); list_iter != _id_list.end(); ++list_iter) {

     stream << "--------- DBMeasurementId: " << id++ << " ------------" << endl;
     stream << *(*list_iter);
  }  
}
//______________________________________________________________________________


