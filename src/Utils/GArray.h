//____________________________________________________________________________
/*!

\class    genie::GArray1D, genie::GArray2D

\brief    C++ arrays wrapped in an class deriving from TNamed.
          A simple workaround for adding C++ arrays in TObjArrays and ROOT
          files in a straightforward way.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Laboratory

\created  March 25, 2015

\cpright  Copyright (c) 2003-2016, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _GENIE_ARRAYS_H_
#define _GENIE_ARRAYS_H_

#include <iostream>
#include <string>
#include <TNamed.h>

class TRootIOCtor;
class TH1D;
class TH2D;

using std::string;

namespace genie {

class GArray1D;
class GArray2D;

// .......................................
// 1-dimensional array of doubles
// .......................................

std::ostream & operator << (std::ostream & stream, const GArray1D & arr);

class GArray1D: public TNamed
{
public:
  GArray1D();
  GArray1D(int n, double xinit=0.);
  GArray1D(const GArray1D & arr);
  GArray1D(const GArray1D & arr, string new_name);
  GArray1D(TRootIOCtor*);
 ~GArray1D();

  using TNamed::Print; //suppress clang 'hides overloaded virtual function [-Woverloaded-virtual]' warnings

  bool   InRange  (int i) const;
  bool   InRange2 (int i) const;
  double Get      (int i) const;
  void   Add      (int i, double x);
  void   Set      (int i, double x);
  void   SetAll   (double x);
  void   Scale    (double s);
  void   Print    (std::ostream & stream) const;

  double GetMaxBin(void) const;

  double operator () (int i);

  int      Size;  
  double * Data; //[Size]

ClassDef(GArray1D,1)
};

// .......................................
// 2-dimensional array of doubles
// .......................................

std::ostream & operator << (std::ostream & stream, const GArray2D & arr);

class GArray2D: public TNamed
{
public:
  GArray2D();
  GArray2D(int n0, int n1, double xinit=0.);
  GArray2D(const GArray2D & arr);
  GArray2D(const GArray2D & arr, string new_name);
  GArray2D(TRootIOCtor*);
 ~GArray2D();

  using TNamed::Print; //suppress clang 'hides overloaded virtual function [-Woverloaded-virtual]' warnings

  bool   InRange  (int i, int j) const;
  bool   InRange2 (int i, int j) const;
  double Get      (int i, int j) const;
  void   Add      (int i, int j, double x);
  void   Set      (int i, int j, double x);
  void   SetAll   (double x);
  void   Scale    (double s);
  void   Print    (std::ostream & stream) const;

  double operator () (int i, int j);

  int      Size;  
  int      Size0;  
  int      Size1;  
  double*  Data; //[Size]

ClassDef(GArray2D,1)
};

namespace utils {
namespace arr {

  // Convert input `central_values' GArray1D to a TH1D (variable-size binning supported)
  TH1D * GArray2TH1D(
    string name, const genie::GArray1D* binning, const genie::GArray1D* central_values, const genie::GArray1D* abs_err, Option_t * opt);

  // Get a slice of the input `central_values' GArray2D and convert it to a TH1D (variable-size binning supported)
  // For example:
  // - to create 1D slice g(i) from 2D array f(i,j) for j=2, use fix_idx = 2, fix_dim = 1
  // - to create 1D slice g(j) from 2D array f(i,j) for i=3, use fix_idx = 3, fix_dim = 0
  TH1D * GArray2TH1D(
    string name, const genie::GArray1D* binning, const genie::GArray2D* central_values, const genie::GArray2D* abs_err, int fix_idx, int fix_dim, Option_t * opt);

  // Convert input `central_values' GArray2D to a TH2D (variable-size binning supported)
  TH2D * GArray2TH2D(
    string name, const genie::GArray1D* binning_x, const genie::GArray1D* binning_y, const genie::GArray2D* central_values, const genie::GArray2D* abs_err, Option_t * opt);

  // In cases where the array describes bin boundaries, find the bin (range: [0, N-1]) that contains the input value.
  // Return -1 for under/overlows.
  int  BinIdx (double x, const genie::GArray1D* binning);
  // In cases where the array describes bin boundaries, return the size of the bith that includes x, or the size of the i^th bin
  // Return 0 for under/overlows.
  double  BinSize (double x, const genie::GArray1D* binning);
  double  BinSize (int idx,  const genie::GArray1D* binning);
  // Fill the input 1-D or 2-D array 
  // If wght_with_bin_sz is true, the weight is scaled by the bin size (wght/bin_size is added in the histogram).
  bool Fill(double x, double wght, bool wght_with_bin_sz, const genie::GArray1D* binning, genie::GArray1D* contents);
  bool Fill(double x, double y, double wght, bool wght_with_bin_sz, const genie::GArray1D* binning_x, const genie::GArray1D* binning_y, genie::GArray2D* contents);
  bool Fill(double x, double y, double wght, bool wght_with_bin_sz, const genie::GArray1D* xlow, const genie::GArray1D* xhigh, const genie::GArray1D* ylow, const genie::GArray1D* yhigh, genie::GArray1D* contents);

}      // arr   namespace
}      // utils namespace
}      // genie namespace

#endif // _GENIE_ARRAYS_H_

