//____________________________________________________________________________
/*!

\class    genie::FunctionMap

\brief    Utility class to hold the values of a scalar function estimated on
          the points of a multi-dimensional uniform grid.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  July 01, 2004

\cpright  Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _FUNCTION_MAP_H_
#define _FUNCTION_MAP_H_

#include <vector>
#include <map>

using std::vector;
using std::map;

namespace genie {

class UnifGrid;

class FunctionMap
{
public:
  FunctionMap();
  FunctionMap(const UnifGrid & grid);
  virtual ~FunctionMap();

  //-- access & modify the grid
  const UnifGrid & GetGrid (void) const { return *fGrid; }
  void  IncreaseGridDensity (unsigned int np, int in_dim=-1);

  //-- set the function value for a point at a given grid position
  //   or coordinate vector
  void  SetValue (double value, const vector<unsigned int> & pos);
  void  SetValue (double value, const vector<double> & pos);

  //-- get value at the input grid position or coordinate vector
  double Value (const vector<unsigned int> & pos) const;
  double Value (const vector<double> & pos) const;

  //-- checks whether the function value for the point or the input
  //   coordinate vector is already set
  bool   ValueIsSet (const vector<unsigned int> & pos) const;
  bool   ValueIsSet (const vector<double> & pos) const;

private:

  void   SetValue   (double value, long int uid);
  double Value      (long int uid) const;
  bool   ValueIsSet (long int uid) const;

  UnifGrid *            fGrid;
  map<long int, double> fFuncMap;
  //map<long int, bool>   fIsSet;
};

}        // namespace

#endif   // _FUNCTION_MAP_H_
