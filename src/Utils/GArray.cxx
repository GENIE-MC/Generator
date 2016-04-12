//____________________________________________________________________________
/*
 Copyright (c) 2003-2016, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Laboratory

 For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________

//#define __skip_range_check__

#include <TRootIOCtor.h>
#include <TMath.h>
#include <TH1D.h>
#include <TH2D.h>

#include "Utils/GArray.h"
#include "Messenger/Messenger.h"

using namespace genie;
using std::endl;
using std::ostream;

ClassImp(GArray1D);
ClassImp(GArray2D);

namespace genie
{
  ostream & operator << (ostream & stream, const GArray1D & arr)
  {
     arr.Print(stream);
     return stream;
  }
  ostream & operator << (ostream & stream, const GArray2D & arr)
  {
     arr.Print(stream);
     return stream;
  }
}

//____________________________________________________________________________
//                               GArray1D
//____________________________________________________________________________
GArray1D::GArray1D() :
TNamed()
{
  this->Size = 0;
  this->Data = 0;
}
//____________________________________________________________________________
GArray1D::GArray1D(int n, double xinit) :
TNamed()
{  
  this->Size = 0;
  this->Data = 0;

  if(n>0) 
  {
    this->Size = n;
    this->Data = new double[n];
  }
  this->SetAll(xinit);
}
//____________________________________________________________________________
GArray1D::GArray1D(const GArray1D & arr) :
TNamed()
{
  this->Size = arr.Size;
  this->Data = new double[arr.Size];

  for(int idx=0; idx < arr.Size; idx++) {
    this->Data[idx] = arr.Data[idx];
  }

  this->SetName(arr.GetName());
}
//____________________________________________________________________________
GArray1D::GArray1D(const GArray1D & arr, string new_name) :
TNamed()
{
  this->Size = arr.Size;
  this->Data = new double[arr.Size];

  for(int idx=0; idx < arr.Size; idx++) {
    this->Data[idx] = arr.Data[idx];
  }

  this->SetName(new_name.c_str());
}
//____________________________________________________________________________
GArray1D::GArray1D(TRootIOCtor*)
{
  this->Size = 0;
  this->Data = 0;
}
//____________________________________________________________________________
GArray1D::~GArray1D()
{
  if(this->Data) { delete [] this->Data; }
}
//____________________________________________________________________________
bool GArray1D::InRange(int idx) const
{
  if(idx >= 0 && idx < this->Size) 
  {
     return true;
  } 
  LOG("GArray", pERROR)
     << "Accessing element " << idx 
     << " in array with elements in range [0, " << this->Size << ")";
  return false;
}
//____________________________________________________________________________
bool GArray1D::InRange2(int idx) const
{
  if(idx >= 0 && idx < this->Size) 
  {
     return true;
  } 
  return false;
}
//____________________________________________________________________________
double GArray1D::Get(int idx) const
{
#ifndef __skip_range_check__
  if(!(this->InRange(idx))) return -9999999;
#endif

  return this->Data[idx];
}
//____________________________________________________________________________
void GArray1D::Add(int idx, double x)
{
  double x_old = this->Data[idx];
  double x_new = x_old + x;
  this->Data[idx] = x_new;
}
//____________________________________________________________________________
void GArray1D::Set(int idx, double x)
{
#ifndef __skip_range_check__
  if(!(this->InRange(idx))) return;
#endif

  this->Data[idx] = x;
}
//____________________________________________________________________________
void GArray1D::SetAll(double x)
{
  for (int idx = 0; idx < this->Size; idx++) {
    this->Set(idx, x);
  }
}
//____________________________________________________________________________
void GArray1D::Scale(double s)
{
  for (int idx = 0; idx < this->Size; idx++) {
    double x = this->Get(idx);
    this->Set(idx, s*x);
  }
}
//____________________________________________________________________________
double GArray1D::GetMaxBin(void) const // needed?
{
  return this->Get(Size-1);
}
//____________________________________________________________________________
void genie::GArray1D::Print(ostream & stream) const
{
  stream << endl;
  for (int idx = 0; idx < this->Size; idx++) {
    double x = this->Get(idx);
    stream << idx << " --> " << x << endl;
  }
}
//____________________________________________________________________________
double genie::GArray1D::operator () (int idx) 
{
  return this->Get(idx);
}
//____________________________________________________________________________
//                               GArray2D
//____________________________________________________________________________
GArray2D::GArray2D() :
TNamed()
{
  this->Size  = 0;
  this->Size0 = 0;
  this->Size1 = 0;
  this->Data  = 0;
}
//____________________________________________________________________________
GArray2D::GArray2D(int n0, int n1, double xinit) :
TNamed()
{  
  this->Size  = 0;
  this->Size0 = 0;
  this->Size1 = 0;
  this->Data  = 0;

  if(n0>0 && n1>0) 
  {
    this->Size  = n0*n1;
    this->Size0 = n0;
    this->Size1 = n1;
    this->Data  = new double[this->Size];
  }
  this->SetAll(xinit);
}
//____________________________________________________________________________
GArray2D::GArray2D(const GArray2D & arr) :
TNamed()
{  
  this->Size  = arr.Size;
  this->Size0 = arr.Size0;
  this->Size1 = arr.Size1;
  this->Data  = new double[arr.Size];

  for(int idx=0; idx < arr.Size; idx++) {
    this->Data[idx] = arr.Data[idx];
  }

  this->SetName(arr.GetName());
}
//____________________________________________________________________________
GArray2D::GArray2D(const GArray2D & arr, string new_name) :
TNamed()
{  
  this->Size  = arr.Size;
  this->Size0 = arr.Size0;
  this->Size1 = arr.Size1;
  this->Data  = new double[arr.Size];

  for(int idx=0; idx < arr.Size; idx++) {
    this->Data[idx] = arr.Data[idx];
  }

  this->SetName(new_name.c_str());
}
//____________________________________________________________________________
GArray2D::GArray2D(TRootIOCtor*)
{
  this->Size  = 0;
  this->Size0 = 0;
  this->Size1 = 0;
  this->Data  = 0;
} 
//____________________________________________________________________________ 
GArray2D::~GArray2D() 
{
  if(this->Data) { delete [] this->Data; }
}
//____________________________________________________________________________
bool GArray2D::InRange(int i, int j) const
{
  if(i >= 0 && i < this->Size0) 
  {
    if(j >= 0 && j < this->Size1) 
    {
      return true;
    }
  } 
  LOG("GArray", pERROR)
    << "Accessing element i,j: " << i << "," << j 
    << " in 2-D array with elements in range "
    << "[0, " <<  this->Size0 << "), [0, " << this->Size1 << ")";
  return false;
}
//____________________________________________________________________________
bool GArray2D::InRange2(int i, int j) const
{
  if(i >= 0 && i < this->Size0) 
  {
    if(j >= 0 && j < this->Size1) 
    {
      return true;
    }
  } 
  return false;
}
//____________________________________________________________________________
double GArray2D::Get(int i, int j) const
{
#ifndef __skip_range_check__
  if(!(this->InRange(i,j))) return -9999999;
#endif

  int idx = i * this->Size1 + j;
  return this->Data[idx];
}
//____________________________________________________________________________
void GArray2D::Add(int i, int j, double x)
{
  int idx = i * this->Size1 + j;
  double x_old = this->Data[idx];
  double x_new = x_old + x;
  this->Data[idx] = x_new;
}
//____________________________________________________________________________
void GArray2D::Set(int i, int j, double x)
{
#ifndef __skip_range_check__
  if(!(this->InRange(i,j))) return;
#endif

  int idx = i * this->Size1 + j;
  this->Data[idx] = x;
}
//____________________________________________________________________________
void GArray2D::SetAll(double x)
{
  for (int idx = 0; idx < this->Size; idx++) {
     this->Data[idx] = x;
  }
}
//____________________________________________________________________________
void GArray2D::Scale(double s)
{
  for (int idx = 0; idx < this->Size; idx++) {
     double x = this->Data[idx];
     this->Data[idx] = s*x;
  }
}
//____________________________________________________________________________
void genie::GArray2D::Print(ostream & /*stream*/) const
{

}
//____________________________________________________________________________
double genie::GArray2D::operator () (int i, int j) 
{
  return this->Get(i,j);
}
//____________________________________________________________________________
//____________________________________________________________________________
//____________________________________________________________________________
TH1D * genie::utils::arr::GArray2TH1D(
  string name, const genie::GArray1D* binning, 
  const genie::GArray1D* central_values, const genie::GArray1D* abs_err, Option_t * opt)
{
  if(!binning) return 0;
  if(!central_values) return 0;

  const int nbins = binning->Size - 1;

  string op = (opt) ? opt : "";
  bool pos_only = (op.find("suppress-nonpositive") != string::npos);

  TH1D * h = new TH1D(name.c_str(),"", nbins, binning->Data);

  for(int ibin = 1; ibin <= nbins; ibin++) 
  { 
    int idx = ibin - 1;
    double val = central_values->Get(idx);
    if(val <= 0. && pos_only) continue;
    h->SetBinContent(ibin, val);   
    if(abs_err) {
       double err = abs_err->Get(idx);
       h->SetBinError(ibin, err);   
    }
  }
  return h;
}
//____________________________________________________________________________
TH1D * genie::utils::arr::GArray2TH1D(
   string name, const genie::GArray1D* binning, 
   const genie::GArray2D* central_values, const genie::GArray2D* abs_err,
   int fix_idx, int fix_dim, Option_t * opt)
{
  if(!binning) return 0;
  if(!central_values) return 0;

  if(fix_dim < 0 || fix_dim >= 2) return 0;

  const int nbins = binning->Size - 1;

  string op = (opt) ? opt : "";
  bool pos_only = (op.find("suppress-nonpositive") != string::npos);

  TH1D * h = new TH1D(name.c_str(),"", nbins, binning->Data);

  for(int ibin = 1; ibin <= nbins; ibin++) 
  { 
    int idx = ibin - 1;
    double val = (fix_dim == 0) ? 
      central_values->Get(fix_idx, idx) : 
      central_values->Get(idx, fix_idx);
    if(val <= 0. && pos_only) continue;
    h->SetBinContent(ibin, val);   
    if(abs_err) {
      double err = (fix_dim == 0) ?
         abs_err->Get(fix_idx, idx) : 
         abs_err->Get(idx, fix_idx);
      h->SetBinError(ibin, err);   
    }
  }
  return h;
}
//____________________________________________________________________________
TH2D * genie::utils::arr::GArray2TH2D(
   string name, const genie::GArray1D* binning_x, const genie::GArray1D* binning_y, 
   const genie::GArray2D* central_values, const genie::GArray2D* abs_err, Option_t * opt)
{
  if(!binning_x) return 0;
  if(!binning_y) return 0;
  if(!central_values) return 0;

  const int nbins_x = binning_x->Size - 1;
  const int nbins_y = binning_y->Size - 1;

  string op = (opt) ? opt : "";
  bool pos_only = (op.find("suppress-nonpositive") != string::npos);

  TH2D * h = new TH2D(
    name.c_str(),"", nbins_x, binning_x->Data, nbins_y, binning_y->Data);

  for(int ibin_x = 1; ibin_x <= nbins_x; ibin_x++) 
  {
    int idx_x = ibin_x - 1;
    for(int ibin_y = 1; ibin_y <= nbins_y; ibin_y++) 
    { 
      int idx_y = ibin_y - 1;
      double val = central_values->Get(idx_x, idx_y);
      if(val <= 0. && pos_only) continue;
      h->SetBinContent(ibin_x, ibin_y, val);   
      if(abs_err) {
         double err = abs_err->Get(idx_x, idx_y);
         h->SetBinError(ibin_x, ibin_y, err);   
      }
    }
  }
  return h;
}
//____________________________________________________________________________
int genie::utils::arr::BinIdx(double x, const genie::GArray1D* binning)
{
  if(!binning) return -1;
  const int n = binning->Size;
  if (x < binning->Get(0) || x > binning->Get(n-1)) return -1;
  return TMath::BinarySearch(n, binning->Data, x); 
}
//____________________________________________________________________________
double genie::utils::arr::BinSize (double x, const genie::GArray1D* binning)
{
  int idx = genie::utils::arr::BinIdx(x, binning);
  return genie::utils::arr::BinSize(idx, binning);
}
//____________________________________________________________________________
double genie::utils::arr::BinSize (int idx,  const genie::GArray1D* binning)
{
  if(!binning) return 0;
  if(idx < 0 ) return 0;  
  const int n = binning->Size;
  if(idx + 1 > n - 1) return 0.;
  double size = binning->Get(idx+1) - binning->Get(idx); 
  return size;
}
//____________________________________________________________________________
bool genie::utils::arr::Fill(
  double x, double wght, bool wght_with_bin_sz,
  const genie::GArray1D* binning, genie::GArray1D* contents)
{
  int idx = genie::utils::arr::BinIdx(x, binning);
  //LOG("GArray", pDEBUG) << "Bin edges: " << *binning;
  //LOG("GArray", pDEBUG) << "Entry " << x << " will be stored at element = " << idx;
  if(idx == -1) return false;
  if(!contents) return false;
  double scale = 1;
  if(wght_with_bin_sz) {
    double sz = genie::utils::arr::BinSize(idx,binning);
    if(sz > 0) {
       scale = 1/sz;
    }
  }
  contents->Add(idx,scale*wght);
  return true;
}
//____________________________________________________________________________
bool genie::utils::arr::Fill(
  double x, double y, double wght, bool wght_with_bin_sz,
  const genie::GArray1D* binning_x, const genie::GArray1D* binning_y, 
  genie::GArray2D* contents)
{
  int idx_x = genie::utils::arr::BinIdx(x, binning_x);
  int idx_y = genie::utils::arr::BinIdx(y, binning_y);
  if(idx_x == -1 || idx_y == -1) return false;
  if(!contents) return false;
  double scale = 1;
  if(wght_with_bin_sz) {
    double sz_x = genie::utils::arr::BinSize(idx_x,binning_x);
    double sz_y = genie::utils::arr::BinSize(idx_y,binning_y);
    double sz   = sz_x * sz_y;
    if(sz > 0) {
       scale = 1/sz;
    }
  }
  contents->Add(idx_x,idx_y,scale*wght);
  return true;
}
//____________________________________________________________________________

