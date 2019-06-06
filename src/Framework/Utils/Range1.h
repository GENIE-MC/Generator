//__________________________________________________________________________
/*!

\class    genie::Range1F_t
\brief    A simple [min,max] interval for floats.

\class    genie::Range1D_t
\brief    A simple [min,max] interval for doubles.

\class    genie::Range1I_t
\brief    A simple [min,max] interval for integers.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  May 06, 2004

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//__________________________________________________________________________

#ifndef _RANGE_1_T_H_
#define _RANGE_1_T_H_

namespace genie {

class Range1F_t
{
public:
  Range1F_t  (void);
  Range1F_t  (float _min, float _max);
  Range1F_t  (const Range1F_t & r);
  ~Range1F_t (void);

  void Copy  (const Range1F_t & r);

  float min;
  float max;
};

class Range1D_t
{
public:
  Range1D_t  (void);
  Range1D_t  (double _min, double _max);
  Range1D_t  (const Range1D_t & r);
  ~Range1D_t (void);

  void Copy  (const Range1D_t & r);

  double min;
  double max;
};

class Range1I_t
{
public:
  Range1I_t  (void);
  Range1I_t  (int _min, int _max);
  Range1I_t  (const Range1I_t & r);
  ~Range1I_t (void);

  void Copy  (const Range1I_t & r);

  int min;
  int max;
};



}      // genie namespace

#endif // _RANGE_1_T_H_

