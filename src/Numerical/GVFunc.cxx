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

#include <cassert>

#include "Numerical/GVFunc.h"

using namespace genie;

//____________________________________________________________________________
GVFunc::GVFunc(unsigned int ndim) :
GFunc(ndim)
{
  assert(ndim>=1);
  fOutV.resize(ndim);
}
//____________________________________________________________________________
GVFunc::~GVFunc()
{

}
//____________________________________________________________________________

