//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab - May 06, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include "Framework/Utils/Range1.h"

using namespace genie;

//____________________________________________________________________________
Range1F_t::Range1F_t(void) :
min(0.),
max(0.)
{

}
//____________________________________________________________________________
Range1F_t::Range1F_t(float _min, float _max) :
min(_min),
max(_max)
{

}
//____________________________________________________________________________
Range1F_t::Range1F_t(const Range1F_t & r) :
min(r.min),
max(r.max)
{

}
//____________________________________________________________________________
Range1F_t::~Range1F_t(void)
{

}
//____________________________________________________________________________
void Range1F_t::Copy(const Range1F_t & r)
{
  min = r.min;
  max = r.max;
}
//____________________________________________________________________________
Range1D_t::Range1D_t(void) :
min(0.),
max(0.)
{

}
//____________________________________________________________________________
Range1D_t::Range1D_t(double _min, double _max) :
min(_min),
max(_max)
{

}
//____________________________________________________________________________
Range1D_t::Range1D_t(const Range1D_t & r) :
min(r.min),
max(r.max)
{

}
//____________________________________________________________________________
Range1D_t::~Range1D_t(void)
{

}
//____________________________________________________________________________
void Range1D_t::Copy(const Range1D_t & r)
{
  min = r.min;
  max = r.max;
}
//____________________________________________________________________________
Range1I_t::Range1I_t(void) :
min(0),
max(0)
{

}
//____________________________________________________________________________
Range1I_t::Range1I_t(int _min, int _max) :
min(_min),
max(_max)
{

}
//____________________________________________________________________________
Range1I_t::Range1I_t(const Range1I_t & r) :
min(r.min),
max(r.max)
{

}
//____________________________________________________________________________
Range1I_t::~Range1I_t(void)
{

}
//____________________________________________________________________________
void Range1I_t::Copy(const Range1I_t & r)
{
  min = r.min;
  max = r.max;
}
//____________________________________________________________________________


