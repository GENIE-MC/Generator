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

#ifndef _RANGE_1_T_H_
#define _RANGE_1_T_H_

namespace genie {

class Range1F_t
{
public:
  float min;
  float max;
};

class Range1D_t
{
public:
  double min;
  double max;
};

class Range1I_t
{
public:
  int min;
  int max;
};



}      // genie namespace

#endif // _RANGE_1_T_H_

