//____________________________________________________________________________
/*!

\class    genie::mc_vs_data::NuXSecData

\brief    Utility class to handle neutrino cross-section data archive

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  Apr 17, 2012

\cpright  Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _NU_XSEC_DATA_H_
#define _NU_XSEC_DATA_H_

#include <string>
#include <vector>

#include <TFile.h>
#include <TTree.h>
#include <TGraphAsymmErrors.h>

using std::string;
using std::vector;

namespace genie {
namespace mc_vs_data {

const int buffer_size = 100;

const int kNMaxDataSets = 25; // max number of datasets in single plot

const int kDataPointStyle[kNMaxDataSets] = 
{ 
  20,      20,    20,       20,    20,
  21,      21,    21,       21,    21,
  24,      24,    24,       24,    24,
  25,      25,    25,       25,    25,
  29,      29,    29,       29,    29
};
const int kDataPointColor[kNMaxDataSets] = 
{
  kBlack,  kRed,  kGreen+1, kBlue, kMagenta+1, 
  kBlack,  kRed,  kGreen+1, kBlue, kMagenta+1, 
  kBlack,  kRed,  kGreen+1, kBlue, kMagenta+1, 
  kBlack,  kRed,  kGreen+1, kBlue, kMagenta+1, 
  kBlack,  kRed,  kGreen+1, kBlue, kMagenta+1
};

class NuXSecData
{
public:
  NuXSecData();
 ~NuXSecData();

  bool OpenArchive(string data_archive_file_name);

  vector<TGraphAsymmErrors *> 
     Retrieve(string keys, double Emin=0., double Emax=200., bool scale_with_E=false);

private:
  
  void Init    (void);
  void CleanUp (void);

  TFile * fNuXSecDataFile;
  TTree * fNuXSecDataTree;

  // tree branches
  char   fDataset  [buffer_size];
  char   fCitation [buffer_size];
  double fE;
  double fEmin;
  double fEmax;
  double fXSecDatum;
  double fXSecDatumErrP;
  double fXSecDatumErrM;
};

} // mc_vs_data namespace
} // genie namepace

#endif  // _NU_XSEC_DATA_H_
