//_________________________________________________________
/*!

\class    genie::Range1F_t
\brief    A simple [min,max] interval for floats.

\class    genie::Range1D_t
\brief    A simple [min,max] interval for doubles.

\class    genie::Range1I_t
\brief    A simple [min,max] interval for integers.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 06, 2004

*/
//_________________________________________________________

#include "Utils/Range1.h"

using namespace genie;

//_________________________________________________________
Range1F_t::Range1F_t(void) :
min(0.),
max(0.)
{

}
//_________________________________________________________
Range1F_t::Range1F_t(float _min, float _max) :
min(_min),
max(_max)
{

}
//_________________________________________________________
Range1F_t::Range1F_t(const Range1F_t & r) :
min(r.min),
max(r.max)
{

}
//_________________________________________________________
Range1F_t::~Range1F_t(void)
{

}
//_________________________________________________________
void Range1F_t::Copy(const Range1F_t & r)
{
  min = r.min;
  max = r.max;
}
//_________________________________________________________
Range1D_t::Range1D_t(void) :
min(0.),
max(0.)
{

}
//_________________________________________________________
Range1D_t::Range1D_t(double _min, double _max) :
min(_min),
max(_max)
{

}
//_________________________________________________________
Range1D_t::Range1D_t(const Range1D_t & r) :
min(r.min),
max(r.max)
{

}
//_________________________________________________________
Range1D_t::~Range1D_t(void)
{

}
//_________________________________________________________
void Range1D_t::Copy(const Range1D_t & r)
{
  min = r.min;
  max = r.max;
}
//_________________________________________________________
Range1I_t::Range1I_t(void) :
min(0),
max(0)
{

}
//_________________________________________________________
Range1I_t::Range1I_t(int _min, int _max) :
min(_min),
max(_max)
{

}
//_________________________________________________________
Range1I_t::Range1I_t(const Range1I_t & r) :
min(r.min),
max(r.max)
{

}
//_________________________________________________________
Range1I_t::~Range1I_t(void)
{

}
//_________________________________________________________
void Range1I_t::Copy(const Range1I_t & r)
{
  min = r.min;
  max = r.max;
}
//_________________________________________________________


