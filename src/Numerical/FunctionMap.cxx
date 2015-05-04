//____________________________________________________________________________
/*
 Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory - July 01, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

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
  //fIsSet.insert( map<long int, bool>::value_type(uid, true) );
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
  map<long int, double>::const_iterator iter = fFuncMap.find(uid);
  if( iter == fFuncMap.end() ) {
    LOG("FMap", pWARN) << "No point with UId = " << uid;
    return 0.;
  } else return iter->second;
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
  //  map<long int, bool>::const_iterator iter = fIsSet.find(uid);
  //  if( iter == fIsSet.end() ) return false;
  //  else return iter->second;

  map<long int, double>::const_iterator iter = fFuncMap.find(uid);
  if( iter == fFuncMap.end() ) return false;
  else return true;
}
//____________________________________________________________________________
void FunctionMap::IncreaseGridDensity(unsigned int np, int in_dim)
{
// np must be 2**N+1
//

  unsigned int nd = fGrid->GetNDimensions();

  // copy of old grid
  UnifGrid * old_grid = new UnifGrid(*fGrid);

  // increase grid density
  if(in_dim==-1) {
    LOG("FMap", pINFO) << "Increasing density: " << np << " points / dim";
    for(unsigned int i = 0; i < nd; i++) {
      fGrid->GridDimension(i)->ReBin(np);
    }
  } else {
    LOG("FMap", pINFO)
     << "Increasing density: " << np << " points @ " << in_dim << "th dim";
    unsigned int i = (unsigned int) in_dim;
    fGrid->GridDimension(i)->ReBin(np);
  }

  // new grid
  UnifGrid * new_grid = fGrid;

  LOG("FMap", pDEBUG) << "Old Grid: " << *old_grid;
  LOG("FMap", pDEBUG) << "New Grid: " << *new_grid;

  if(fFuncMap.size() == 0) {
    delete old_grid;
    return;
  }

  // initialize grid position vector
  vector<unsigned int> old_grid_pos(nd), new_grid_pos(nd);
  for(unsigned int id = 0; id < nd; id++) {
    old_grid_pos[id] = 0;
    new_grid_pos[id] = 0;
  }

  // number of points in the old grid - where the function was evaluated
  unsigned int ntot = old_grid->NPoints();
  LOG("FMap", pINFO) << "Number of grid points (old grid) = " << ntot;

  LOG("FMap", pNOTICE)
   << "Moving to new grid: func values associated with old grid points";

  // loop over all old grid points and translate their UIds
  map<long int, double> func_map;
  for(unsigned int i=0; i<ntot; i++) {

    old_grid->Point(i, old_grid_pos); // current point
    long int old_uid = old_grid->UId(old_grid_pos); // old grid UId

    if(in_dim==-1) {
      for(unsigned int i = 0; i < nd; i++) {
	int nmult = (np-1) / ((*old_grid)[i]->NPoints() - 1);
         new_grid_pos[i] = old_grid_pos[i] * nmult;
      }
    } else {
      int nmult = (np-1) / ((*old_grid)[in_dim]->NPoints()-1);
      new_grid_pos[in_dim] = old_grid_pos[in_dim] * nmult;
    }
    long int new_uid = new_grid->UId(new_grid_pos); // new grid UId

    LOG("FMap", pDEBUG) << "UId: " << old_uid << " -> " << new_uid;

    // get previously computed value
    double val = fFuncMap[old_uid];
    // and associate it with its new grid uid
    func_map.insert(map<long int, double>::value_type(new_uid, val));
  }

  // delete old grid
  delete old_grid;
  fFuncMap.clear();

  // refresh fFuncMap and fIsSet
  map<long int, double>::iterator iter;
  for(iter = func_map.begin(); iter != func_map.end(); ++iter) {
    long int uid = iter->first;
    double   val = iter->second;
    this->SetValue(val,uid);
  }
}
//____________________________________________________________________________


