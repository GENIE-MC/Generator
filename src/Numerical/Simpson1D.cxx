//____________________________________________________________________________
/*!

\class    genie::Simpson1D

\brief    

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 03, 2004

*/ 
//____________________________________________________________________________

#include <TMath.h>

#include "Numerical/Simpson1D.h"
#include "Messenger/Messenger.h"

using namespace genie;

//____________________________________________________________________________
Simpson1D::Simpson1D():
IntegratorI("genie::Simpson1D")
{

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
double Simpson1D::EvalError(FunctionMap & func_map) const
{
// If f(x) is continuous in [a,b], then the error in Simpson's rule is no
// larger than err = ((b-a)^5 / 180*n^4)*|fmax|

  const UnifGrid & grid = func_map.GetGrid();  
  assert(grid.GetNDimensions() == 1);

  int    N = grid[0]->npoints;
  double L = grid[0]->max - grid[0]->min;

  double fmax = 0;
  for(int i = 0; i < N; i++) 
      fmax = TMath::Max(fmax, TMath::Abs(func_map.Func(i))); 

  double err = TMath::Power(L,5) * fmax / (180*TMath::Power(N,4));
  return err;
}
//____________________________________________________________________________

