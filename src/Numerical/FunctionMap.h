//____________________________________________________________________________
/*!

\class    genie::FunctionMap

\brief

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  July 01, 2004

*/
//____________________________________________________________________________

#ifndef _FUNCTION_MAP_H_
#define _FUNCTION_MAP_H_

#include <vector>
#include <map>

#include "Numerical/UnifGrid.h"

using std::vector;
using std::map;

namespace genie {

class FunctionMap
{
public:

  FunctionMap();
  FunctionMap(const UnifGrid & grid);
  virtual ~FunctionMap();

  //-- add points to the Function map
  
  void AddPoint (double value, const vector<int> & pos);
  
  //-- get value at a given point

  double Func   (const vector<int> & pos) const;
  
  //-- simpler interface for lower dimensions
  
  void   AddPoint (double value, int iposx, int iposy, int iposz);
  void   AddPoint (double value, int iposx, int iposy);
  void   AddPoint (double value, int iposx);
  
  double Func     (int iposx, int iposy, int iposz) const;
  double Func     (int iposx, int iposy) const;
  double Func     (int iposx) const;
  
  //-- access the grid

  const UnifGrid & GetGrid(void) const { return *fGrid; }
  
private:

  UnifGrid * fGrid;
  map<long int, double> fFuncMap;
};

}        // namespace

#endif   // _FUNCTION_MAP_H_
