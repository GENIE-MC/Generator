//____________________________________________________________________________
/*!

\class    genie::Simpson1D

\brief    

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 03, 2004

*/ 
//____________________________________________________________________________

#include "Numerical/Simpson1D.h"
#include "Messenger/Messenger.h"

using namespace genie;

//____________________________________________________________________________
Simpson1D::Simpson1D():
IntegratorI()
{
  fName     = "genie::Simpson1D";
  fParamSet = "NoConfig";
}
//____________________________________________________________________________
Simpson1D::~Simpson1D()
{

}
//____________________________________________________________________________
double Simpson1D::Integrate(FunctionMap & func_map) const
{
  const UnifGrid & grid = func_map.GetGrid();
  
  assert(grid.GetNDimensions() == 1);
  
  int    N    = grid[0]->npoints;
  double step = grid[0]->step;
  
  LOG("Simpson1D", pDEBUG) << "N-Points = " << N << ", Step-Size = " << step;

  double sum = (func_map.Func(0) + func_map.Func(N-1)) / 2.;

  for(int i = 0; i < N-1; i++)  sum += ( func_map.Func(i) * (i%2 + 1) );

  sum *= (2.*step/3.);

  return sum;
}
//____________________________________________________________________________

