//____________________________________________________________________________
/*!

\class    genie::FunctionMap

\brief

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  July 01, 2004

*/
//____________________________________________________________________________

#include "Messenger/Messenger.h"
#include "Numerical/FunctionMap.h"

using namespace genie;

//____________________________________________________________________________
FunctionMap::FunctionMap()
{

}
//____________________________________________________________________________
FunctionMap::FunctionMap(const UnifGrid & grid) 
{
  fGrid = new UnifGrid(grid);
}
//____________________________________________________________________________
FunctionMap::~FunctionMap()
{
  delete fGrid;
}
//____________________________________________________________________________
void FunctionMap::AddPoint(double value, const vector<int> & pos)
{
  long int uid = fGrid->GridPoint2UId(pos);

  fFuncMap.insert( map<long int, double>::value_type(uid, value) );
}  
//____________________________________________________________________________
void FunctionMap::AddPoint(double value, int iposx, int iposy, int iposz)
{
  vector<int> pos(3);

  pos[0] = iposx;
  pos[1] = iposy;
  pos[2] = iposz;

  long int uid = fGrid->GridPoint2UId(pos);

  fFuncMap.insert( map<long int, double>::value_type(uid, value) );
}
//____________________________________________________________________________
void FunctionMap::AddPoint(double value, int iposx, int iposy)
{
  vector<int> pos(2);

  pos[0] = iposx;
  pos[1] = iposy;

  long int uid = fGrid->GridPoint2UId(pos);

  fFuncMap.insert( map<long int, double>::value_type(uid, value) );
}
//____________________________________________________________________________
void FunctionMap::AddPoint(double value, int iposx)
{
  vector<int> pos(1);

  pos[0] = iposx;

  long int uid = fGrid->GridPoint2UId(pos);

  fFuncMap.insert( map<long int, double>::value_type(uid, value) );
}
//____________________________________________________________________________
double FunctionMap::Func (const vector<int> & pos) const
{
  long int uid = fGrid->GridPoint2UId(pos);

  if( fFuncMap.count(uid) == 1 ) {
      map<long int, double>::const_iterator iter = fFuncMap.find(uid);
      return iter->second;
  } else
      LOG("Numerical", pWARN) << "No point with UId = " << uid;

  return 0;
}
//____________________________________________________________________________
double FunctionMap::Func (int iposx, int iposy, int iposz) const
{
  vector<int> pos;

  pos.push_back(iposx);
  pos.push_back(iposy);
  pos.push_back(iposz);

  return Func(pos);
}
//____________________________________________________________________________
double FunctionMap::Func (int iposx, int iposy) const
{
  vector<int> pos;

  pos.push_back(iposx);
  pos.push_back(iposy);

  return Func(pos);
}
//____________________________________________________________________________
double FunctionMap::Func (int iposx) const
{
  vector<int> pos;

  pos.push_back(iposx);

  return Func(pos);
}
//____________________________________________________________________________
