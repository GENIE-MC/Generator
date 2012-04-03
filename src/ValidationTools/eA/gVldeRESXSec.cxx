//____________________________________________________________________________
/*!

\program gvld_e_res_xsec

\brief   Compares GENIE with resonance electro-production cross-section data.

         The data come from the JLab Hall C archive maintained at:
         https://hallcweb.jlab.org/resdata/database/
         A local copy may be found in:
         $GENIE/data/validation/eA/xsec/differential/res/
         The archive contains ~12k data points.

         Syntax:
           gvld_e_res_xsec -m model [-d data_archive_location]

         Options:

           [] Denotes an optional argument.

           -m Specify GENIE resonance electro-production cross-section model.
              If not set, only data -no GENIE predictions- will be displayed.
           -d Full path to the electron scattering archive.
              By default, will pick the one at:
              $GENIE/data/validation/eA/xsec/differential/res/eRES.root

         Example:

            % gvld_e_res_xsec -m genie::ReinSeghalRESPXSec/Default

\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory

\created Oct 16, 2009 

\cpright Copyright (c) 2003-2011, GENIE Neutrino MC Generator Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#include <cstdlib>
#include <cassert>
#include <sstream>
#include <string>

#include <TSystem.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TPostScript.h>
#include <TLatex.h>
#include <TH1D.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TPavesText.h>
#include <TText.h>
#include <TLegend.h>
#include <TBox.h>

#include "Algorithm/AlgFactory.h"
#include "Base/XSecAlgorithmI.h"
#include "BaryonResonance/BaryonResonance.h"
#include "BaryonResonance/BaryonResUtils.h"
//#include "Conventions/GBuild.h"
#include "Conventions/Constants.h"
#include "Conventions/Units.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGUtils.h"
#include "PDG/PDGCodes.h"
#include "Utils/CmdLnArgParser.h"
#include "Utils/StringUtils.h"
#include "Utils/SystemUtils.h"
#include "Utils/Style.h"

using std::ostringstream;
using std::string;

using namespace genie;
using namespace genie::constants;

//_________________________________________________________________________________
// Utility class to hold info on plotted datasets
class eResDataSetDescription_t
{
public:
  eResDataSetDescription_t(
     int tgtpdg, string expt, double E, double theta, bool show) :
    fTgtPdg (tgtpdg), 
    fExpt   (expt), 
    fE      (E), 
    fTheta  (theta),
    fShow   (show)
  {   
  }
  eResDataSetDescription_t() 
  {   
  }
  int    TgtPdg   (void) const { return fTgtPdg; }
  int    TgtZ     (void) const { return pdg::IonPdgCodeToZ(fTgtPdg); }
  int    TgtA     (void) const { return pdg::IonPdgCodeToA(fTgtPdg); }
  string TgtName  (void) const { 
    // all data are either in Hydrogen or Deuterium
    if(fTgtPdg == 1000010010) return "Hydrogen";
    if(fTgtPdg == 1000010020) return "Deuterium";
    return "Other";
  }
  string Expt     (void) const { return fExpt; }
  double E        (void) const { return fE; }
  double Theta    (void) const { return fTheta; }
  bool   Show     (void) const { return fShow; }
  string LabelTeX (void) const { 
    ostringstream label;
    label << fExpt << " (" << this->TgtName() << "), ";
    label << "E = " << fE << " GeV, ";
    label << "#theta = " << fTheta << "^{o}";
    return label.str();
  }
private:
  int    fTgtPdg;  //
  string fExpt;    //
  double fE;       //
  double fTheta;   //
  bool   fShow;    //
};
//_________________________________________________________________________________

/* 
..............................................................................
ELECTRON-NUCLEUS RESONANCE SCATTERING DATASET DESCRIPTIONS
..............................................................................
*/
const int kNumOfDataSets = 286;

eResDataSetDescription_t * kDataSet[kNumOfDataSets] =
{
   new eResDataSetDescription_t(1000010010, "jlab_e00_002",   2.2375,  33.130, true), // 0
   new eResDataSetDescription_t(1000010010, "jlab_e00_002",   2.2375,  41.320, true), // 1
   new eResDataSetDescription_t(1000010010, "jlab_e00_002",   2.2375,  50.610, true), // 2
   new eResDataSetDescription_t(1000010010, "jlab_e00_002",   3.0440,  12.005, true), // 3
   new eResDataSetDescription_t(1000010010, "jlab_e00_002",   3.0440,  19.830, true), // 4
   new eResDataSetDescription_t(1000010010, "jlab_e00_002",   3.0440,  22.710, true), // 5
   new eResDataSetDescription_t(1000010010, "jlab_e00_002",   3.0440,  24.900, true), // 6
   new eResDataSetDescription_t(1000010010, "jlab_e00_002",   3.0440,  32.300, true), // 7
   new eResDataSetDescription_t(1000010010, "jlab_e00_002",   3.0440,  36.190, true), // 8
   new eResDataSetDescription_t(1000010010, "jlab_e00_002",   4.4125,  13.510, true), // 9
   new eResDataSetDescription_t(1000010010, "jlab_e00_002",   4.4125,  11.210, true), // 10
   new eResDataSetDescription_t(1000010010, "jlab_e00_002",   4.4125,  16.500, true), // 11
   new eResDataSetDescription_t(1000010010, "jlab_e00_002",   5.5005,  11.210, true), // 12

   new eResDataSetDescription_t(1000010010, "jlab_e94_110",   1.150,   47.950, true), // 13
   new eResDataSetDescription_t(1000010010, "jlab_e94_110",   1.150,   59.970, true), // 14
   new eResDataSetDescription_t(1000010010, "jlab_e94_110",   1.884,   33.930, true), // 15
   new eResDataSetDescription_t(1000010010, "jlab_e94_110",   1.884,   47.940, true), // 16
   new eResDataSetDescription_t(1000010010, "jlab_e94_110",   2.238,   21.950, true), // 17
   new eResDataSetDescription_t(1000010010, "jlab_e94_110",   2.238,   31.930, true), // 18
   new eResDataSetDescription_t(1000010010, "jlab_e94_110",   2.238,   42.950, true), // 19
   new eResDataSetDescription_t(1000010010, "jlab_e94_110",   2.238,   58.950, true), // 20
   new eResDataSetDescription_t(1000010010, "jlab_e94_110",   2.238,   79.950, true), // 21
   new eResDataSetDescription_t(1000010010, "jlab_e94_110",   3.118,   12.450, true), // 22
   new eResDataSetDescription_t(1000010010, "jlab_e94_110",   3.118,   15.950, true), // 23
   new eResDataSetDescription_t(1000010010, "jlab_e94_110",   3.118,   19.440, true), // 24
   new eResDataSetDescription_t(1000010010, "jlab_e94_110",   3.118,   22.950, true), // 25
   new eResDataSetDescription_t(1000010010, "jlab_e94_110",   3.118,   32.950, true), // 26
   new eResDataSetDescription_t(1000010010, "jlab_e94_110",   3.118,   40.950, true), // 27
   new eResDataSetDescription_t(1000010010, "jlab_e94_110",   3.118,   49.950, true), // 28
   new eResDataSetDescription_t(1000010010, "jlab_e94_110",   3.118,   61.950, true), // 29
   new eResDataSetDescription_t(1000010010, "jlab_e94_110",   3.118,   77.950, true), // 30
   new eResDataSetDescription_t(1000010010, "jlab_e94_110",   4.110,   38.950, true), // 31
   new eResDataSetDescription_t(1000010010, "jlab_e94_110",   4.412,   38.940, true), // 32
   new eResDataSetDescription_t(1000010010, "jlab_e94_110",   4.412,   44.960, true), // 33
   new eResDataSetDescription_t(1000010010, "jlab_e94_110",   4.412,   50.950, true), // 34
   new eResDataSetDescription_t(1000010010, "jlab_e94_110",   5.498,   12.970, true), // 35
   new eResDataSetDescription_t(1000010010, "jlab_e94_110",   5.498,   15.460, true), // 36
   new eResDataSetDescription_t(1000010010, "jlab_e94_110",   5.498,   17.940, true), // 37
   new eResDataSetDescription_t(1000010010, "jlab_e94_110",   5.498,   20.450, true), // 38
   new eResDataSetDescription_t(1000010010, "jlab_e94_110",   5.498,   22.950, true), // 39
   new eResDataSetDescription_t(1000010010, "jlab_e94_110",   5.498,   25.450, true), // 40
   new eResDataSetDescription_t(1000010010, "jlab_e94_110",   5.498,   27.950, true), // 41
   new eResDataSetDescription_t(1000010010, "jlab_e94_110",   5.498,   30.450, true), // 42
   new eResDataSetDescription_t(1000010010, "jlab_e94_110",   5.498,   32.970, true), // 43
   new eResDataSetDescription_t(1000010010, "jlab_e94_110",   5.498,   35.460, true), // 44

   new eResDataSetDescription_t(1000010010, "jlab_ioanna",    2.445,   20.000, true), // 45
   new eResDataSetDescription_t(1000010010, "jlab_ioanna",    2.445,   30.000, true), // 46
   new eResDataSetDescription_t(1000010010, "jlab_ioanna",    2.445,   38.500, true), // 47
   new eResDataSetDescription_t(1000010010, "jlab_ioanna",    2.445,   70.010, true), // 48
   new eResDataSetDescription_t(1000010010, "jlab_ioanna",    3.245,   26.980, true), // 49
   new eResDataSetDescription_t(1000010010, "jlab_ioanna",    4.045,   47.990, true), // 50
   new eResDataSetDescription_t(1000010010, "jlab_ioanna",    4.054,   24.030, true), // 51
   new eResDataSetDescription_t(1000010010, "jlab_ioanna",    4.054,   30.000, true), // 52
   new eResDataSetDescription_t(1000010010, "jlab_ioanna",    4.054,   39.990, true), // 53

   new eResDataSetDescription_t(1000010010, "slac_e133",      9.772,   10.000, true), // 54
   new eResDataSetDescription_t(1000010010, "slac_e133",     12.604,   10.000, true), // 55
   new eResDataSetDescription_t(1000010010, "slac_e133",     15.740,   10.000, true), // 56
   new eResDataSetDescription_t(1000010010, "slac_e133",     18.520,   10.000, true), // 57
   new eResDataSetDescription_t(1000010010, "slac_e133",     20.998,   10.000, true), // 58

   new eResDataSetDescription_t(1000010010, "slac_e140",      4.502,   34.009, true), // 59
   new eResDataSetDescription_t(1000010010, "slac_e140",      4.504,   18.020, true), // 60
   new eResDataSetDescription_t(1000010010, "slac_e140",      4.506,   26.015, true), // 61
   new eResDataSetDescription_t(1000010010, "slac_e140",      5.808,   34.009, true), // 62
   new eResDataSetDescription_t(1000010010, "slac_e140",      6.500,   59.999, true), // 63
   new eResDataSetDescription_t(1000010010, "slac_e140",      6.500,   18.000, true), // 64
   new eResDataSetDescription_t(1000010010, "slac_e140",      6.509,   18.020, true), // 65
   new eResDataSetDescription_t(1000010010, "slac_e140",      6.711,   26.015, true), // 66
   new eResDataSetDescription_t(1000010010, "slac_e140",      7.000,   49.998, true), // 67
   new eResDataSetDescription_t(1000010010, "slac_e140",      7.019,   10.000, true), // 68
   new eResDataSetDescription_t(1000010010, "slac_e140",      7.912,   34.009, true), // 69
   new eResDataSetDescription_t(1000010010, "slac_e140",      8.614,   18.020, true), // 70
   new eResDataSetDescription_t(1000010010, "slac_e140",      8.700,   25.993, true), // 71
   new eResDataSetDescription_t(1000010010, "slac_e140",      8.713,   26.015, true), // 72
   new eResDataSetDescription_t(1000010010, "slac_e140",      9.022,   10.000, true), // 73
   new eResDataSetDescription_t(1000010010, "slac_e140",      9.999,   14.991, true), // 74
   new eResDataSetDescription_t(1000010010, "slac_e140",      9.999,   33.992, true), // 75
   new eResDataSetDescription_t(1000010010, "slac_e140",     10.014,   34.009, true), // 76
   new eResDataSetDescription_t(1000010010, "slac_e140",     10.027,    5.988, true), // 77
   new eResDataSetDescription_t(1000010010, "slac_e140",     10.392,   18.020, true), // 78
   new eResDataSetDescription_t(1000010010, "slac_e140",     10.400,   60.000, true), // 79
   new eResDataSetDescription_t(1000010010, "slac_e140",     10.400,   18.000, true), // 80
   new eResDataSetDescription_t(1000010010, "slac_e140",     10.998,   10.000, true), // 81
   new eResDataSetDescription_t(1000010010, "slac_e140",     11.884,   25.993, true), // 82
   new eResDataSetDescription_t(1000010010, "slac_e140",     11.900,   26.015, true), // 83
   new eResDataSetDescription_t(1000010010, "slac_e140",     12.497,   14.991, true), // 84
   new eResDataSetDescription_t(1000010010, "slac_e140",     12.498,   18.996, true), // 85
   new eResDataSetDescription_t(1000010010, "slac_e140",     13.000,    4.000, true), // 86
   new eResDataSetDescription_t(1000010010, "slac_e140",     13.290,   60.000, true), // 87
   new eResDataSetDescription_t(1000010010, "slac_e140",     13.300,   15.000, true), // 88
   new eResDataSetDescription_t(1000010010, "slac_e140",     13.300,   18.000, true), // 89
   new eResDataSetDescription_t(1000010010, "slac_e140",     13.320,   18.020, true), // 90
   new eResDataSetDescription_t(1000010010, "slac_e140",     13.500,    6.000, true), // 91
   new eResDataSetDescription_t(1000010010, "slac_e140",     13.545,   10.000, true), // 92
   new eResDataSetDescription_t(1000010010, "slac_e140",     13.549,    5.988, true), // 93
   new eResDataSetDescription_t(1000010010, "slac_e140",     15.004,   14.991, true), // 94
   new eResDataSetDescription_t(1000010010, "slac_e140",     15.004,   18.996, true), // 95
   new eResDataSetDescription_t(1000010010, "slac_e140",     15.005,   25.993, true), // 96
   new eResDataSetDescription_t(1000010010, "slac_e140",     15.022,   26.015, true), // 97
   new eResDataSetDescription_t(1000010010, "slac_e140",     15.204,   10.000, true), // 98
   new eResDataSetDescription_t(1000010010, "slac_e140",     16.000,    4.000, true), // 99
   new eResDataSetDescription_t(1000010010, "slac_e140",     16.000,    6.000, true), // 100
   new eResDataSetDescription_t(1000010010, "slac_e140",     16.000,   15.000, true), // 101
   new eResDataSetDescription_t(1000010010, "slac_e140",     16.000,   18.000, true), // 102
   new eResDataSetDescription_t(1000010010, "slac_e140",     16.075,    5.988, true), // 103
   new eResDataSetDescription_t(1000010010, "slac_e140",     17.706,   10.000, true), // 104
   new eResDataSetDescription_t(1000010010, "slac_e140",     18.010,    4.000, true), // 105
   new eResDataSetDescription_t(1000010010, "slac_e140",     18.018,   18.996, true), // 106
   new eResDataSetDescription_t(1000010010, "slac_e140",     18.018,   25.993, true), // 107
   new eResDataSetDescription_t(1000010010, "slac_e140",     19.350,   10.000, true), // 108
   new eResDataSetDescription_t(1000010010, "slac_e140",     19.500,    6.000, true), // 109
   new eResDataSetDescription_t(1000010010, "slac_e140",     19.500,   10.000, true), // 110
   new eResDataSetDescription_t(1000010010, "slac_e140",     19.500,   15.000, true), // 111
   new eResDataSetDescription_t(1000010010, "slac_e140",     19.500,   18.000, true), // 112
   new eResDataSetDescription_t(1000010010, "slac_e140",     19.500,   20.600, true), // 113
   new eResDataSetDescription_t(1000010010, "slac_e140",     19.544,    5.988, true), // 114
   new eResDataSetDescription_t(1000010010, "slac_e140",     20.001,   18.996, true), // 115
   new eResDataSetDescription_t(1000010010, "slac_e140",     20.005,    4.000, true), // 116

   new eResDataSetDescription_t(1000010010, "slac_e140x",     1.997,   46.160, true), // 117
   new eResDataSetDescription_t(1000010010, "slac_e140x",     3.301,   50.000, true), // 118
   new eResDataSetDescription_t(1000010010, "slac_e140x",     4.014,   50.700, true), // 119
   new eResDataSetDescription_t(1000010010, "slac_e140x",     4.951,   56.010, true), // 120
   new eResDataSetDescription_t(1000010010, "slac_e140x",     5.600,   58.000, true), // 121
   new eResDataSetDescription_t(1000010010, "slac_e140x",     6.453,   59.510, true), // 122

   new eResDataSetDescription_t(1000010010, "slac_e49a10",    4.892,   10.000, true), // 123
   new eResDataSetDescription_t(1000010010, "slac_e49a10",    7.019,   10.000, true), // 124
   new eResDataSetDescription_t(1000010010, "slac_e49a10",    9.022,   10.000, true), // 125
   new eResDataSetDescription_t(1000010010, "slac_e49a10",   10.998,   10.000, true), // 126
   new eResDataSetDescription_t(1000010010, "slac_e49a10",   13.545,   10.000, true), // 127
   new eResDataSetDescription_t(1000010010, "slac_e49a10",   15.204,   10.000, true), // 128
   new eResDataSetDescription_t(1000010010, "slac_e49a10",   17.706,   10.000, true), // 129
   new eResDataSetDescription_t(1000010010, "slac_e49a10",   19.350,   10.000, true), // 130

   new eResDataSetDescription_t(1000010010, "slac_e49a6",     4.511,    5.988, true), // 131
   new eResDataSetDescription_t(1000010010, "slac_e49a6",     7.014,    5.988, true), // 132
   new eResDataSetDescription_t(1000010010, "slac_e49a6",    10.027,    5.988, true), // 133
   new eResDataSetDescription_t(1000010010, "slac_e49a6",    13.549,    5.988, true), // 134
   new eResDataSetDescription_t(1000010010, "slac_e49a6",    16.075,    5.988, true), // 135
   new eResDataSetDescription_t(1000010010, "slac_e49a6",    19.544,    5.988, true), // 136

   new eResDataSetDescription_t(1000010010, "slac_e49b",      4.502,   34.009, true), // 137
   new eResDataSetDescription_t(1000010010, "slac_e49b",      4.504,   18.020, true), // 138
   new eResDataSetDescription_t(1000010010, "slac_e49b",      4.506,   26.015, true), // 139
   new eResDataSetDescription_t(1000010010, "slac_e49b",      5.808,   34.009, true), // 140
   new eResDataSetDescription_t(1000010010, "slac_e49b",      6.509,   18.020, true), // 141
   new eResDataSetDescription_t(1000010010, "slac_e49b",      6.711,   26.015, true), // 142
   new eResDataSetDescription_t(1000010010, "slac_e49b",      7.912,   34.009, true), // 143
   new eResDataSetDescription_t(1000010010, "slac_e49b",      8.614,   18.020, true), // 144
   new eResDataSetDescription_t(1000010010, "slac_e49b",      8.713,   26.015, true), // 145
   new eResDataSetDescription_t(1000010010, "slac_e49b",     10.014,   34.009, true), // 146
   new eResDataSetDescription_t(1000010010, "slac_e49b",     10.392,   18.020, true), // 147
   new eResDataSetDescription_t(1000010010, "slac_e49b",     11.900,   26.015, true), // 148
   new eResDataSetDescription_t(1000010010, "slac_e49b",     12.518,   18.020, true), // 149
   new eResDataSetDescription_t(1000010010, "slac_e49b",     12.518,   34.009, true), // 150
   new eResDataSetDescription_t(1000010010, "slac_e49b",     13.320,   18.020, true), // 151
   new eResDataSetDescription_t(1000010010, "slac_e49b",     15.022,   26.015, true), // 152
   new eResDataSetDescription_t(1000010010, "slac_e49b",     17.027,   18.020, true), // 153

   new eResDataSetDescription_t(1000010010, "slac_e61",       4.499,    4.000, true), // 154
   new eResDataSetDescription_t(1000010010, "slac_e61",       7.000,    4.000, true), // 155
   new eResDataSetDescription_t(1000010010, "slac_e61",       9.993,    4.000, true), // 156
   new eResDataSetDescription_t(1000010010, "slac_e61",      13.000,    4.000, true), // 157
   new eResDataSetDescription_t(1000010010, "slac_e61",      16.000,    4.000, true), // 158
   new eResDataSetDescription_t(1000010010, "slac_e61",      18.010,    4.000, true), // 159
   new eResDataSetDescription_t(1000010010, "slac_e61",      20.005,    4.000, true), // 160

   new eResDataSetDescription_t(1000010010, "slac_e87",       4.500,   25.994, true), // 161
   new eResDataSetDescription_t(1000010010, "slac_e87",       4.501,   33.993, true), // 162
   new eResDataSetDescription_t(1000010010, "slac_e87",       6.700,   25.994, true), // 163
   new eResDataSetDescription_t(1000010010, "slac_e87",       8.700,   25.994, true), // 164
   new eResDataSetDescription_t(1000010010, "slac_e87",       9.999,   14.991, true), // 165
   new eResDataSetDescription_t(1000010010, "slac_e87",       9.999,   33.993, true), // 166
   new eResDataSetDescription_t(1000010010, "slac_e87",      11.883,   25.994, true), // 167
   new eResDataSetDescription_t(1000010010, "slac_e87",      12.497,   14.991, true), // 168
   new eResDataSetDescription_t(1000010010, "slac_e87",      12.498,   18.997, true), // 169
   new eResDataSetDescription_t(1000010010, "slac_e87",      15.004,   14.991, true), // 170
   new eResDataSetDescription_t(1000010010, "slac_e87",      15.004,   18.997, true), // 171
   new eResDataSetDescription_t(1000010010, "slac_e87",      15.005,   25.994, true), // 172
   new eResDataSetDescription_t(1000010010, "slac_e87",      18.018,   18.997, true), // 173
   new eResDataSetDescription_t(1000010010, "slac_e87",      18.018,   25.994, true), // 174
   new eResDataSetDescription_t(1000010010, "slac_e87",      20.001,   18.997, true), // 175

   new eResDataSetDescription_t(1000010010, "slac_e891",      6.500,   59.999, true), // 176
   new eResDataSetDescription_t(1000010010, "slac_e891",      7.000,   49.998, true), // 177
   new eResDataSetDescription_t(1000010010, "slac_e891",     10.400,   60.000, true), // 178
   new eResDataSetDescription_t(1000010010, "slac_e891",     13.290,   60.000, true), // 179
   new eResDataSetDescription_t(1000010010, "slac_e891",     16.002,   59.999, true), // 180
   new eResDataSetDescription_t(1000010010, "slac_e891",     19.505,   60.000, true), // 181

   new eResDataSetDescription_t(1000010010, "slac_e8920",     6.500,   18.000, true), // 182
   new eResDataSetDescription_t(1000010010, "slac_e8920",     7.000,    6.000, true), // 183
   new eResDataSetDescription_t(1000010010, "slac_e8920",    10.400,   18.000, true), // 184
   new eResDataSetDescription_t(1000010010, "slac_e8920",    13.300,   15.000, true), // 185
   new eResDataSetDescription_t(1000010010, "slac_e8920",    13.300,   18.000, true), // 186
   new eResDataSetDescription_t(1000010010, "slac_e8920",    13.500,    6.000, true), // 187
   new eResDataSetDescription_t(1000010010, "slac_e8920",    16.000,    6.000, true), // 188
   new eResDataSetDescription_t(1000010010, "slac_e8920",    16.000,   15.000, true), // 189
   new eResDataSetDescription_t(1000010010, "slac_e8920",    16.000,   18.000, true), // 190
   new eResDataSetDescription_t(1000010010, "slac_e8920",    19.500,    6.000, true), // 191
   new eResDataSetDescription_t(1000010010, "slac_e8920",    19.500,   10.000, true), // 192
   new eResDataSetDescription_t(1000010010, "slac_e8920",    19.500,   15.000, true), // 193
   new eResDataSetDescription_t(1000010010, "slac_e8920",    19.500,   18.000, true), // 194
   new eResDataSetDescription_t(1000010010, "slac_e8920",    19.500,   20.600, true), // 195

   new eResDataSetDescription_t(1000010010, "slac_ne11",      5.507,   15.146, true), // 196
   new eResDataSetDescription_t(1000010010, "slac_ne11",      5.507,   18.981, true), // 197
   new eResDataSetDescription_t(1000010010, "slac_ne11",      5.507,   22.805, true), // 198
   new eResDataSetDescription_t(1000010010, "slac_ne11",      5.507,   26.823, true), // 199
   new eResDataSetDescription_t(1000010010, "slac_ne11",      9.800,   13.248, true), // 200
   new eResDataSetDescription_t(1000010010, "slac_ne11",      9.800,   15.367, true), // 201
   new eResDataSetDescription_t(1000010010, "slac_ne11",      9.800,   17.516, true), // 202
   new eResDataSetDescription_t(1000010010, "slac_ne11",      9.800,   19.753, true), // 203

   new eResDataSetDescription_t(1000010010, "slac_onen1haf",  5.000,    1.482, true), // 204
   new eResDataSetDescription_t(1000010010, "slac_onen1haf",  7.103,    1.482, true), // 205
   new eResDataSetDescription_t(1000010010, "slac_onen1haf",  9.301,    1.482, true), // 206
   new eResDataSetDescription_t(1000010010, "slac_onen1haf", 11.799,    1.482, true), // 207
   new eResDataSetDescription_t(1000010010, "slac_onen1haf", 13.805,    1.482, true), // 208
   new eResDataSetDescription_t(1000010010, "slac_onen1haf", 15.604,    1.482, true), // 209
   new eResDataSetDescription_t(1000010010, "slac_onen1haf", 17.298,    1.482, true), // 210
   new eResDataSetDescription_t(1000010010, "slac_onen1haf", 18.703,    1.482, true), // 211
   new eResDataSetDescription_t(1000010010, "slac_onen1haf", 19.995,    1.482, true), // 212

   new eResDataSetDescription_t(1000010020, "jlab_ioanna",    2.445,   20.000, true), // 213
   new eResDataSetDescription_t(1000010020, "jlab_ioanna",    2.445,   30.000, true), // 214
   new eResDataSetDescription_t(1000010020, "jlab_ioanna",    2.445,   70.010, true), // 215
   new eResDataSetDescription_t(1000010020, "jlab_ioanna",    3.245,   26.980, true), // 216
   new eResDataSetDescription_t(1000010020, "jlab_ioanna",    4.045,   30.000, true), // 217
   new eResDataSetDescription_t(1000010020, "jlab_ioanna",    4.045,   39.990, true), // 218
   new eResDataSetDescription_t(1000010020, "jlab_ioanna",    4.045,   48.000, true), // 219
   new eResDataSetDescription_t(1000010020, "jlab_ioanna",    4.054,   24.030, true), // 220

   new eResDataSetDescription_t(1000010020, "slac_e140x",     3.301,   50.000, true), // 221
   new eResDataSetDescription_t(1000010020, "slac_e140x",     4.014,   50.700, true), // 222
   new eResDataSetDescription_t(1000010020, "slac_e140x",     4.951,   56.010, true), // 223
   new eResDataSetDescription_t(1000010020, "slac_e140x",     5.600,   58.000, true), // 224
   new eResDataSetDescription_t(1000010020, "slac_e140x",     6.453,   59.510, true), // 225

   new eResDataSetDescription_t(1000010020, "slac_e49a10",    4.892,   10.000, true), // 226
   new eResDataSetDescription_t(1000010020, "slac_e49a10",    7.019,   10.000, true), // 227
   new eResDataSetDescription_t(1000010020, "slac_e49a10",    9.022,   10.000, true), // 228
   new eResDataSetDescription_t(1000010020, "slac_e49a10",   10.998,   10.000, true), // 229
   new eResDataSetDescription_t(1000010020, "slac_e49a10",   13.545,   10.000, true), // 230
   new eResDataSetDescription_t(1000010020, "slac_e49a10",   15.204,   10.000, true), // 231
   new eResDataSetDescription_t(1000010020, "slac_e49a10",   17.706,   10.000, true), // 232
   new eResDataSetDescription_t(1000010020, "slac_e49a10",   19.350,   10.000, true), // 233

   new eResDataSetDescription_t(1000010020, "slac_e49a6",     4.511,    5.988, true), // 234
   new eResDataSetDescription_t(1000010020, "slac_e49a6",     7.014,    5.988, true), // 235
   new eResDataSetDescription_t(1000010020, "slac_e49a6",    10.027,    5.988, true), // 236
   new eResDataSetDescription_t(1000010020, "slac_e49a6",    13.549,    5.988, true), // 237
   new eResDataSetDescription_t(1000010020, "slac_e49a6",    16.075,    5.988, true), // 238
   new eResDataSetDescription_t(1000010020, "slac_e49a6",    19.544,    5.988, true), // 239

   new eResDataSetDescription_t(1000010020, "slac_e49b",      4.502,   34.009, true), // 240
   new eResDataSetDescription_t(1000010020, "slac_e49b",      4.504,   18.020, true), // 241
   new eResDataSetDescription_t(1000010020, "slac_e49b",      4.506,   26.015, true), // 242
   new eResDataSetDescription_t(1000010020, "slac_e49b",      5.808,   34.009, true), // 243
   new eResDataSetDescription_t(1000010020, "slac_e49b",      6.509,   18.020, true), // 244
   new eResDataSetDescription_t(1000010020, "slac_e49b",      6.711,   26.015, true), // 245
   new eResDataSetDescription_t(1000010020, "slac_e49b",      7.912,   34.009, true), // 246
   new eResDataSetDescription_t(1000010020, "slac_e49b",      8.614,   18.020, true), // 247
   new eResDataSetDescription_t(1000010020, "slac_e49b",      8.713,   26.015, true), // 248
   new eResDataSetDescription_t(1000010020, "slac_e49b",     10.014,   34.009, true), // 249
   new eResDataSetDescription_t(1000010020, "slac_e49b",     10.392,   18.020, true), // 250
   new eResDataSetDescription_t(1000010020, "slac_e49b",     11.900,   26.015, true), // 251
   new eResDataSetDescription_t(1000010020, "slac_e49b",     12.518,   34.009, true), // 252
   new eResDataSetDescription_t(1000010020, "slac_e49b",     13.320,   18.020, true), // 253
   new eResDataSetDescription_t(1000010020, "slac_e49b",     15.022,   26.015, true), // 254
   new eResDataSetDescription_t(1000010020, "slac_e49b",     17.027,   18.020, true), // 255

   new eResDataSetDescription_t(1000010020, "slac_e61",       4.499,    4.000, true), // 256
   new eResDataSetDescription_t(1000010020, "slac_e61",       6.999,    4.000, true), // 257
   new eResDataSetDescription_t(1000010020, "slac_e61",       9.993,    4.000, true), // 258
   new eResDataSetDescription_t(1000010020, "slac_e61",      12.987,    4.000, true), // 259
   new eResDataSetDescription_t(1000010020, "slac_e61",      16.000,    4.000, true), // 260
   new eResDataSetDescription_t(1000010020, "slac_e61",      18.010,    4.000, true), // 261
   new eResDataSetDescription_t(1000010020, "slac_e61",      20.005,    4.000, true), // 262

   new eResDataSetDescription_t(1000010020, "slac_e891",      6.500,   59.999, true), // 263
   new eResDataSetDescription_t(1000010020, "slac_e891",      7.000,   49.998, true), // 264
   new eResDataSetDescription_t(1000010020, "slac_e891",     10.400,   60.000, true), // 265
   new eResDataSetDescription_t(1000010020, "slac_e891",     13.290,   60.000, true), // 266
   new eResDataSetDescription_t(1000010020, "slac_e891",     16.002,   59.999, true), // 267
   new eResDataSetDescription_t(1000010020, "slac_e891",     19.505,   60.000, true), // 268

   new eResDataSetDescription_t(1000010020, "slac_e8920",     6.500,   18.000, true), // 269
   new eResDataSetDescription_t(1000010020, "slac_e8920",     7.000,    6.000, true), // 270
   new eResDataSetDescription_t(1000010020, "slac_e8920",    10.400,   18.000, true), // 271
   new eResDataSetDescription_t(1000010020, "slac_e8920",    13.300,   15.000, true), // 272
   new eResDataSetDescription_t(1000010020, "slac_e8920",    13.300,   18.000, true), // 273
   new eResDataSetDescription_t(1000010020, "slac_e8920",    13.500,    6.000, true), // 274
   new eResDataSetDescription_t(1000010020, "slac_e8920",    16.000,    6.000, true), // 275
   new eResDataSetDescription_t(1000010020, "slac_e8920",    16.000,   15.000, true), // 276
   new eResDataSetDescription_t(1000010020, "slac_e8920",    16.000,   18.000, true), // 277
   new eResDataSetDescription_t(1000010020, "slac_e8920",    19.500,    6.000, true), // 278
   new eResDataSetDescription_t(1000010020, "slac_e8920",    19.500,    5.000, true), // 279
   new eResDataSetDescription_t(1000010020, "slac_e8920",    19.500,   18.000, true), // 280
   new eResDataSetDescription_t(1000010020, "slac_e8920",    19.500,   20.600, true), // 281

   new eResDataSetDescription_t(1000010020, "slac_ne11",      5.507,   15.146, true), // 282
   new eResDataSetDescription_t(1000010020, "slac_ne11",      5.507,   18.981, true), // 283
   new eResDataSetDescription_t(1000010020, "slac_ne11",      5.507,   22.805, true), // 284
   new eResDataSetDescription_t(1000010020, "slac_ne11",      5.507,   26.823, true)  // 285
};

// function prototypes
void           Init               (void);
void           Plot               (void);
void           End                (void);
void           AddCoverPage       (void);
TGraphErrors * Data               (int iset);
TGraph *       Model              (int iset, int imodel);
void           Draw               (int iset);
void           GetCommandLineArgs (int argc, char ** argv);
void           PrintSyntax        (void);

// command-line arguments
string gOptGModelName   = "";
string gOptGModelConf   = "";
string gOptDataFilename = "";

// dbase information
const char * kDefDataFile = "data/validation/eA/xsec/differential/res/eRES.root";  

// globals
TFile *        gResDataFile  = 0;
TTree *        gResDataTree  = 0;
TPostScript *  gPS           = 0;
TCanvas *      gC            = 0;
bool           gShowModel    = false;

// consts

const int kNCx = 2; // number of columns in TCanvas::Divide()
const int kNCy = 2; // number of rows    in TCanvas::Divide()

const int kNRes=18;  
Resonance_t kResId[kNRes] = {
   kP33_1232, kS11_1535, kD13_1520, kS11_1650,
   kD13_1700, kD15_1675, kS31_1620, kD33_1700,
   kP11_1440, kP33_1600, kP13_1720, kF15_1680,
   kP31_1910, kP33_1920, kF35_1905, kF37_1950,
   kP11_1710, kF17_1970 
};

// current program draws predictions only for the explicit resonance-production
// model at W<Wcut
const bool kDrawHatchcedScalingRegion = true; 

const double kWcut = 1.7; // Wcut from UserPhysicsOptions.xml

//_________________________________________________________________________________
int main(int argc, char ** argv)
{
  GetCommandLineArgs (argc,argv);

  Init();
  Plot();
  End();

  LOG("gvldtest", pINFO)  << "Done!";
  return 0;
}
//_________________________________________________________________________________
void Plot(void)
{
  // loop over data sets
  for(int iset = 0; iset < kNumOfDataSets; iset++) 
  {
    if(!kDataSet[iset]->Show()) continue;
    LOG("gvldtest", pNOTICE) 
      << "Producing plots for: " << kDataSet[iset]->LabelTeX();
    Draw(iset);
  }
}
//_________________________________________________________________________________
void Init(void)
{
  LOG("gvldtest", pNOTICE) << "Initializing...";

  // get TTree with electron-scattering data
  if( ! utils::system::FileExists(gOptDataFilename) ) {
      LOG("gvldtest", pFATAL) 
         << "Can not find file: " << gOptDataFilename;
      gAbortingInErr = true;
      exit(1);
  }
  gResDataFile = new TFile(gOptDataFilename.c_str(),"read");  
  gResDataTree = (TTree *) gResDataFile->Get("resnt");
  if(!gResDataTree) {
      LOG("gvldtest", pFATAL) 
         << "Can not find TTree `resnt' in file: " << gOptDataFilename;
      gAbortingInErr = true;
      exit(1);
  }

  // genie style
  utils::style::SetDefaultStyle();

  // canvas
  gC = new TCanvas("c","",20,20,500,650);
  gC->SetBorderMode(0);
  gC->SetFillColor(0);
  gC->SetGridx();
  gC->SetGridy();

  // output postscript file
  gPS = new TPostScript("eRES.genie_vs_data.ps", 111);

  // cover page
  AddCoverPage();
}
//_________________________________________________________________________________
void AddCoverPage(void)
{
  // header
  gPS->NewPage();
  gC->Range(0,0,100,100);
  TPavesText hdr(10,40,90,70,3,"tr");
  hdr.AddText(" ");
  hdr.AddText("Resonance Electro-Production: GENIE vs data");
  hdr.AddText(" ");
  hdr.Draw();
  gC->Update();
}
//_________________________________________________________________________________
void End(void)
{
  LOG("gvldtest", pNOTICE) << "Cleaning up...";

  gPS->Close();

  delete gC;
  delete gPS;

  gResDataFile->Close();
}
//_________________________________________________________________________________
// Corresponding GENIE prediction for the `iset' data set 
//.................................................................................
TGraph * Model(int iset, int imodel)
{
  if(!gShowModel) return 0;

  LOG("gvldtest", pNOTICE) 
    << "Getting GENIE prediction (model ID = " 
    << imodel << ", data set ID = " << iset << ")";

  AlgFactory * algf = AlgFactory::Instance();
  const XSecAlgorithmI * xsec_alg = 
     dynamic_cast<const XSecAlgorithmI *> (
          algf->GetAlgorithm(gOptGModelName, gOptGModelConf));
  if(!xsec_alg) return 0;

  double M  = kNucleonMass;
  double M2 = M*M;

  double E     = kDataSet[iset]->E();
  double theta = kDataSet[iset]->Theta();
  double costh = TMath::Cos(2*kPi*theta/360.);

  LOG("gvldtest", pINFO) 
     << " E = " << E 
     << ", theta = " << theta << " (cos(theta) = " << costh << ")";

  //int target_pdgc = kDataSet[iset]->TgtPdg();
  int Z = kDataSet[iset]->TgtZ();
  int A = kDataSet[iset]->TgtA();
  int N = A-Z;
  bool   tgt_has_p = (Z>0);
  bool   tgt_has_n = (N>0);
  double frac_p    = (double) Z / (double) A;
  double frac_n    = (double) N / (double) A;

  Interaction * ep_res = 0;
  Interaction * en_res = 0;

  if(tgt_has_p) {
     ep_res = Interaction::RESEM(1000010010, kPdgProton,  kPdgElectron, E);
  }
  if(tgt_has_n) {
     en_res = Interaction::RESEM(1000010010, kPdgNeutron, kPdgElectron, E);
  } //is this right kPdgNeutron???  SAD

  const int n = 500;
  double d2sig_dEpdOmega_array[n];
  double W2_array[n];

  double Epmin = 0.01;
  double Epmax = E;
  double dEp = (Epmax-Epmin)/(n-1);

  int nn0 = 0; // non-zero xsec values

  for(int i=0; i<n; i++) {
     double Ep = Epmin + i*dEp;
     double Q2 = 2*E*Ep*(1-costh);
     double W2 = M2 + 2*M*(E-Ep)-Q2;
     double W  = TMath::Sqrt( TMath::Max(0.,W2) );

     if(W2 <= 0) {
        LOG("gvldtest", pDEBUG) 
          << "Ep = " << Ep << ", Q2 = " << Q2 << ", W2 = " << W2
          << "... Skipping point";
        d2sig_dEpdOmega_array[i] = 0.;
        W2_array[i] = 0.;
        continue;
     }

     if(tgt_has_p) {
       ep_res->KinePtr()->SetW (W);
       ep_res->KinePtr()->SetQ2(Q2);
     }
     if(tgt_has_n) {
       en_res->KinePtr()->SetW (W);
       en_res->KinePtr()->SetQ2(Q2);
     }

     double d2sig_dWdQ2 = 0;
    
     for(int ires=0; ires<kNRes; ires++) {

        double d2sig_dWdQ2_res_p = 0.;
        double d2sig_dWdQ2_res_n = 0.;

        if(tgt_has_p) {
          ep_res->ExclTagPtr()->SetResonance(kResId[ires]);
          d2sig_dWdQ2_res_p = xsec_alg->XSec(ep_res,kPSWQ2fE) / units::nb;
          d2sig_dWdQ2_res_p = TMath::Max(0., d2sig_dWdQ2_res_p);
          LOG("gvldtest", pINFO) 
             << "d2xsec_dWdQ2(ep; " << utils::res::AsString(kResId[ires])
             << "; E = " << E << " GeV, W = " << W << " GeV, Q2 = " << Q2 << " GeV^2"
             << "; Ep = " << Ep << " GeV, theta = " << theta << " deg) = " 
             << d2sig_dWdQ2_res_p << " nbarn/GeV^3";
        }
        if(tgt_has_n) {
          en_res->ExclTagPtr()->SetResonance(kResId[ires]);
          d2sig_dWdQ2_res_n = xsec_alg->XSec(en_res,kPSWQ2fE) / units::nb;
          d2sig_dWdQ2_res_n = TMath::Max(0., d2sig_dWdQ2_res_n);
          LOG("gvldtest", pINFO) 
             << "d2xsec_dWdQ2(en; " << utils::res::AsString(kResId[ires])
             << "; E = " << E << " GeV, W = " << W << " GeV, Q2 = " << Q2 << " GeV^2"
             << "; Ep = " << Ep << " GeV, theta = " << theta << " deg) = " 
             << d2sig_dWdQ2_res_n << " nbarn/GeV^3";
        }

        d2sig_dWdQ2 += (frac_p*d2sig_dWdQ2_res_p + frac_n*d2sig_dWdQ2_res_n);
     }

     // d^2 sigma / dW dQ^2 --> d^2sigma / dE' dOmega
     double jacobian = (E*Ep)*(M+2*E*(1-costh))/(kPi*W);
     double d2sig_dEpdOmega = jacobian * d2sig_dWdQ2;

     if(TMath::IsNaN(d2sig_dEpdOmega)) {
        LOG("gvldtest", pWARN) << "Got a NaN!";
        d2sig_dEpdOmega = 0;
     }

     if(d2sig_dEpdOmega>0) nn0++;

     d2sig_dEpdOmega_array[i] = TMath::Max(0., d2sig_dEpdOmega);
     W2_array[i] = W2;
  }

  if(tgt_has_p) { delete ep_res; }
  if(tgt_has_n) { delete en_res; }

  LOG("gvldtest", pNOTICE) << "Computed " << nn0 << " non-zero xsec values";
  
  TGraph * gr = new TGraph(n,W2_array,d2sig_dEpdOmega_array);

  return gr;
}
//_________________________________________________________________________________
TGraphErrors * Data(int iset)
{
  const double dE      = 1.0E-3;
  const double dtheta  = 2.5E-2;

  double E     = kDataSet[iset]->E();
  double theta = kDataSet[iset]->Theta();

  int Z = kDataSet[iset]->TgtZ();
  int A = kDataSet[iset]->TgtA();

  const char * selection = 
    Form("E > %f && E < %f && theta > %f && theta < %f && Z == %d && A == %d",
         E     - dE,
         E     + dE,
         theta - dtheta,
         theta + dtheta,
         Z,A);

  gResDataTree->Draw("W2:xsec:xsec_err", selection, "goff");
  
  int n = gResDataTree->GetSelectedRows();

  LOG("gvldtest", pNOTICE) 
    << "Found " << n << " data points in the xsec archive";

  if(n == 0) return 0; // return null graph

  // Data returned by TTree::Draw() are not necessarily ordered in W
  // Do the ordering here before building the graph
  int    *  idx = new int   [n];
  double *  xv  = new double[n];
  double *  yv  = new double[n];
  double *  dyv = new double[n];

  TMath::Sort(n,gResDataTree->GetV1(),idx,false);

  for(int i=0; i<n; i++) {
     int ii = idx[i];
     xv [i] = (gResDataTree->GetV1())[ii];
     yv [i] = (gResDataTree->GetV2())[ii];
     dyv[i] = (gResDataTree->GetV3())[ii];
  }

  TGraphErrors * gr = new TGraphErrors(n,xv,yv,0,dyv);

  delete [] idx;
  delete [] xv;
  delete [] yv;
  delete [] dyv;

  return gr;
}
//_________________________________________________________________________________
void Draw(int iset)
{
  // get all measurements for the current channel from the NuValidator MySQL dbase
  TGraphErrors * data = Data(iset);

  // get the corresponding GENIE model prediction
  TGraph * model = Model(iset,0);

  if(!model && !data) return;

  int plots_per_page = kNCx * kNCy;
  int iplot = 1 + iset % plots_per_page;

  if(iplot == 1) {
     gPS->NewPage();
     gC -> Clear();
     gC -> Divide(kNCx,kNCy);
  }

  gC -> GetPad(iplot) -> Range(0,0,100,100);
  gC -> GetPad(iplot) -> SetFillColor(0);
  gC -> GetPad(iplot) -> SetBorderMode(0);
  gC -> GetPad(iplot) -> cd();

  double xmin = 0.0, scale_xmin = 0.5;
  double xmax = 0.0, scale_xmax = 1.2;
  double ymin = 0.0, scale_ymin = 0.4;
  double ymax = 0.0, scale_ymax = 1.2;

  TH1F * hframe = 0;
  bool have_frame = false;

  if(data) {
    xmin  = ( data->GetX() )[TMath::LocMin(data->GetN(),data->GetX())];
    xmax  = ( data->GetX() )[TMath::LocMax(data->GetN(),data->GetX())];
    ymin  = ( data->GetY() )[TMath::LocMin(data->GetN(),data->GetY())];
    ymax  = ( data->GetY() )[TMath::LocMax(data->GetN(),data->GetY())];
    xmin  = TMath::Max(xmin, 0.5); // some data go very low 
    if(model) {
       ymin  = TMath::Min(
        ymin, ( model->GetY() )[TMath::LocMin(model->GetN(),model->GetY())]);
       ymax  = TMath::Max(
        ymax, ( model->GetY() )[TMath::LocMax(model->GetN(),model->GetY())]);
    }
    hframe = (TH1F*) gC->GetPad(iplot)->DrawFrame(
        scale_xmin*xmin, scale_ymin*ymin, scale_xmax*xmax, scale_ymax*ymax);
    have_frame = true;
    utils::style::Format(data, 1,1,1,1,8,0.8);
    data->Draw("P");
  }//data?

  if(model) {
    if(!have_frame) {
       xmin  = ( model->GetX() )[TMath::LocMin(model->GetN(),model->GetX())];
       xmax  = ( model->GetX() )[TMath::LocMax(model->GetN(),model->GetX())];
       ymin  = ( model->GetY() )[TMath::LocMin(model->GetN(),model->GetY())];
       ymax  = ( model->GetY() )[TMath::LocMax(model->GetN(),model->GetY())];
       hframe = (TH1F*) gC->GetPad(iplot)->DrawFrame(
         scale_xmin*xmin, scale_ymin*ymin, scale_xmax*xmax, scale_ymax*ymax);
    }
    utils::style::Format(model, 1,1,1,1,1,1);
    model->Draw("L");
  }

  //hframe->Draw();
  hframe->GetXaxis()->SetTitle("W^{2} (GeV^{2})");
  hframe->GetYaxis()->SetTitle("d^{2}#sigma / d#Omega dE (nb/sr/GeV)");

  // scaling region
  TBox * scaling_region = 0;
  if(kDrawHatchcedScalingRegion) {
    double W2c = kWcut*kWcut;
    if(W2c > scale_xmin*xmin && W2c < scale_xmax*xmax) {
       scaling_region = new TBox(
           W2c, scale_ymin*ymin, scale_xmax*xmax, scale_ymax*ymax);
       scaling_region->SetFillColor(kRed);
       scaling_region->SetFillStyle(3005);
       scaling_region->Draw();
    }
  }

  // some data show the elastic peak - mark the are to avoid confusion
  if(xmin < 1) {
    double Wm2 = 1.21; // between the QE and Delta peaks
    TBox * qe_peak = new TBox(
       scale_xmin*xmin, scale_ymin*ymin, Wm2, scale_ymax*ymax);
     qe_peak->SetFillColor(kBlue);
     qe_peak->SetFillStyle(3005);
     qe_peak->Draw();
  }

  // title
  TLatex * title = new TLatex(
     scale_xmin*xmin + 0.2*(scale_xmax*xmax-scale_xmin*xmin),
    1.01*scale_ymax*ymax,kDataSet[iset]->LabelTeX().c_str());
  title->SetTextSize(0.027);
  title->Draw();

  gC->GetPad(iplot)->Update();
  gC->Update();
}
//_________________________________________________________________________________
// Parsing command-line arguments, check/form filenames, etc
//.................................................................................
void GetCommandLineArgs(int argc, char ** argv)
{
  LOG("gvldtest", pNOTICE) << "*** Parsing command line arguments";

  CmdLnArgParser parser(argc,argv);

  // get the name and configu
  if(parser.OptionExists('d')){
     string filename = parser.ArgAsString('d');
     gOptDataFilename = filename;
  } else {
     if(gSystem->Getenv("GENIE")) {
        string base_dir = string( gSystem->Getenv("GENIE") );
        string filename = base_dir + kDefDataFile;
        gOptDataFilename = filename;
     } else { 
        LOG("gvldtest", pFATAL) 
          << "\n Please make sure that $GENIE is defined, or use the -d option"
          << "\n You didn't specify a data file and I can not pick the default one either";
        gAbortingInErr = true;
        exit(1);
     }
  }

  // get the name and configuration for the GENIE model to be tested
  if(parser.OptionExists('m')){
     string model = parser.ArgAsString('m');
     vector<string> modelv = utils::str::Split(model,"/");
     assert(modelv.size()==2);
     gOptGModelName = modelv[0];
     gOptGModelConf = modelv[1];
     gShowModel = true;
  } else {
     gShowModel = false;
  }
}
//_________________________________________________________________________________
void PrintSyntax(void)
{
  LOG("gvldtest", pNOTICE)
    << "\n\n" << "Syntax:" << "\n"
    << "  gvld_e_res_xsec -m model [-d data_archive_location]\n";
}
//_________________________________________________________________________________
