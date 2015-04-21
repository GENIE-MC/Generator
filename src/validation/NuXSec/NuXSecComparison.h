//____________________________________________________________________________
/*!

\class    genie::mc_vs_data::NuXSecComparison

\brief    Utility class to hold the information needed so as to setup a
          neutrino cross-section data/MC comparison.
          This info includes:
          - which datasets to extract from the archive,
          - how to compute the corresponding GENIE cross-section prediction,
          - how to format the generated comparison plot

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  Apr 17, 2012

\cpright  Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _NU_XSEC_COMPARISON_H_
#define _NU_XSEC_COMPARISON_H_

#include <string>

#include "validation/NuXSec/NuXSecFunc.h"

using std::string;

namespace genie {
namespace mc_vs_data {

class NuXSecComparison
{
public:
  NuXSecComparison(
    string id, string label, 
    string dataset_keys, NuXSecFunc * xsec_func,
    double Emin,  double Emax, 
    bool in_logx, bool in_logy,  bool scale_with_E = false);
 ~NuXSecComparison();

  string       ID          (void) const { return fID;            }
  string       Label       (void) const { return fLabel;         }
  string       DataSetKeys (void) const { return fDataSetKeys;   }
  NuXSecFunc * XSecFunc    (void) const { return fXSecFunc;      }
  double       Emin        (void) const { return fEmin;          }
  double       Emax        (void) const { return fEmax;          }
  bool         InLogX      (void) const { return fInLogX;        }
  bool         InLogY      (void) const { return fInLogY;        }
  bool         ScaleWithE  (void) const { return fScaleWithE;    }

  bool SameMCPrediction(const NuXSecComparison * c);

private:

  string       fID;               // A data/MC comparison ID
  string       fLabel;            // (Roo)TeX label for data/MC comparison plot 
  string       fDataSetKeys;      // Semicolon-separated list of keys for all datasets included in comparison 
  NuXSecFunc * fXSecFunc;         // Cross-section calculator
  double       fEmin;             // Minimum energy plotted 
  double       fEmax;             // Maximum energy plotted 
  bool         fInLogX;           // Decide whether to show the x axis in linear or log scale 
  bool         fInLogY;           // Decide whether to show the y axis in linear or log scale 
  bool         fScaleWithE;       // Decide whether to plot sigma or sigma/E 
  
};

} // mc_vs_data namespace
} // genie namepace

#endif  // _NU_XSEC_COMPARISON
