//____________________________________________________________________________
/*!

\program gtestBLI2DUnifGrid

\brief   Program used for testing / debugging GENIE's BLI2DUnifGrid

\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

\created May 29, 2009

\cpright Copyright (c) 2003-2016, GENIE Neutrino MC Generator Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#include <TMath.h>
#include <TFile.h>
#include <TNtuple.h>

#include "Messenger/Messenger.h"
#include "Numerical/BLI2D.h"
#include "Numerical/RandomGen.h"

using namespace genie;

double func(double x, double y);

int main(int /*argc*/, char ** /*argv*/)
{
  int npoints=10000;

  int nx = 100;
  int ny = 100;

  double xmin = -5;
  double xmax =  5;
  double ymin = -5;
  double ymax =  5;

  double dx = (xmax-xmin)/(nx-1);
  double dy = (ymax-ymin)/(ny-1);

  BLI2DUnifGrid biln(nx,xmin,xmax,ny,ymin,ymax);

  for(int ix=0; ix<nx; ix++) {
      double x = xmin + ix * dx;
      for(int iy=0; iy<ny; iy++) {
         double y = ymin + iy * dy;      
         double z = func(x,y);
         biln.AddPoint(x,y,z);
      }
  }

  RandomGen * rnd = RandomGen::Instance();

  TNtuple * nt = new TNtuple("nt","billinear interpolation validation","x:y:ztrue:zeval");

  for(int ip=0; ip<npoints; ip++) {
   double rx = rnd->RndGen().Uniform();
   double ry = rnd->RndGen().Uniform();
   double x  = xmin + (xmax-xmin)*rx;
   double y  = ymin + (ymax-ymin)*ry;
   double zt = func(x,y);
   double ze = biln.Evaluate(x,y);

   nt->Fill(x,y,zt,ze);
  }

  TFile f("./bli2dug.root","recreate");
  nt->Write();
  f.Close();

  LOG("test", pINFO)  << "Done!";
  return 0;
}


double func(double x, double y)
{
  //return x*y+10*x-7*y+1;
  return TMath::Sin(x)/x * TMath::Sin(y)/y;
}
