//____________________________________________________________________________
/*!

\class    genie::FunctionMap

\brief    Utility class to hold the values of a scalar function estimated on
          the points of a multi-dimensional uniform grid.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  July 01, 2004

*/
//____________________________________________________________________________

#include "Messenger/Messenger.h"
#include "Numerical/FunctionMap.h"
#include "Numerical/UnifGrid.h"
#include "Numerical/UnifGridDimension.h"

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
void FunctionMap::SetValue(double value, const vector<unsigned int> & pos)
{
  long int uid = fGrid->UId(pos);
  this->SetValue(value, uid);
}
//____________________________________________________________________________
void FunctionMap::SetValue(double value, const vector<double> & pos)
{
  long int uid = fGrid->UId(pos);
  this->SetValue(value, uid);
}
//____________________________________________________________________________
void FunctionMap::SetValue(double value, long int uid)
{
  fFuncMap.insert( map<long int, double>::value_type(uid, value) );
  fIsSet.insert( map<long int, bool>::value_type(uid, true) );
}
//____________________________________________________________________________
double FunctionMap::Value(const vector<unsigned int> & pos) const
{
  long int uid = fGrid->UId(pos);
  return this->Value(uid);
}
//____________________________________________________________________________
double FunctionMap::Value(const vector<double> & pos) const
{
  long int uid = fGrid->UId(pos);
  return this->Value(uid);
}
//____________________________________________________________________________
double FunctionMap::Value(long int uid) const
{
  if( fFuncMap.count(uid) == 1 ) {
      map<long int, double>::const_iterator iter = fFuncMap.find(uid);
      return iter->second;
  } else {
      LOG("FMap", pWARN) << "No point with UId = " << uid;
  }
  return 0;
}
//____________________________________________________________________________
bool FunctionMap::ValueIsSet(const vector<unsigned int> & pos) const
{
  long int uid = fGrid->UId(pos);
  return this->ValueIsSet(uid);
}
//____________________________________________________________________________
bool FunctionMap::ValueIsSet(const vector<double> & pos) const
{
  long int uid = fGrid->UId(pos);
  return this->ValueIsSet(uid);
}
//____________________________________________________________________________
bool FunctionMap::ValueIsSet(long int uid) const
{
  if( fIsSet.count(uid) == 1 ) {
      map<long int, bool>::const_iterator iter = fIsSet.find(uid);
      return iter->second;
  } else
      return false;
}
//____________________________________________________________________________
void FunctionMap::IncreaseGridDensity(unsigned int np)
{
  LOG("FMap", pINFO) << "Increasing grid density to " << np << " points/dim";

  unsigned int nd = fGrid->GetNDimensions();

  // copy of old grid
  UnifGrid * old_grid = new UnifGrid(*fGrid);

  // increase grid density
  for(unsigned int i = 0; i < nd; i++) {
    fGrid->GridDimension(i)->ReBin(np);
  }
  // new grid
  UnifGrid * new_grid = fGrid;

  // initialize grid position vector
  vector<unsigned int> pos(nd);
  for(unsigned int id = 0; id < nd; id++) {
    pos[id] = 0;
  }

  // number of points in the old grid - where the function was evaluated
  unsigned int ntot = old_grid->NPoints();

  LOG("FMap", pINFO) << "Number of grid points = " << ntot;

  // loop over all old grid points and translate their UIds
  map<long int, double> func_map;
  for(unsigned int i=0; i<ntot; i++) {

    old_grid->Point(i, pos); // current point
    long int old_uid = old_grid->UId(pos); // old grid UId
    long int new_uid = new_grid->UId(pos); // new grid UId

    // get previously computed value
    double val = fFuncMap[old_uid];
    // and associate it with its new grid uid
    func_map.insert(map<long int, double>::value_type(new_uid, val));
  }

  // delete old grid
  delete old_grid;
  fFuncMap.clear();
  fIsSet.clear();

  // refresh fFuncMap and fIsSet
  map<long int, double>::iterator iter;
  for(iter = func_map.begin(); iter != func_map.end(); ++iter) {
    long int uid = iter->first;
    double   val = iter->second;
    fFuncMap.insert( map<long int, double>::value_type(uid, val)  );
    fIsSet.insert( map<long int, double>::value_type(uid, true)  );
  }
}
//____________________________________________________________________________


