//____________________________________________________________________________
/*!

\program gvld_e_qel_xsec

\brief   Compares GENIE with quasi-elastic electron scattering data.

         The data come from the Quasielastic Electron Nucleus Scattering Data
         Archive maintained by Donal Day:
         http://http://faculty.virginia.edu/qes-archive/
         See also O. Benhar, D. Day and I. Sick, Rev.Mod.Phys.80, 189, 2008
         A local copy of the data may be found in:
         $GENIE/data/validation/eA/xsec/differential/qe
         The archive contains ~19k data points.

         Syntax:
           gvld_e_qel_xsec [-g input_file_list] [-d data_archive_location]

         Options:

           [] Denotes an optional argument.

           -d Full path to the electron scattering archive.
              By default, will pick the one distributed with GENIE.
           -g An XML file with GENIE inputs (cross sections and event samples).
              If not set, only data -no GENIE predictions- will be displayed.
              Multiple models can be included in the input file, each identified by 
              a "name" (all model predictions will be overlayed).

              <?xml version="1.0" encoding="ISO-8859-1"?>
              <vld_inputs>
                 <model name="a_model_name">
                   <xsec_file>             /path/model_1/xsec.root </xsec_file>
                   <evt_file format="gst"> /path/model_1/ev0.root  </evt_file>
                   <evt_file format="gst"> /path/model_1/ev1.root  </evt_file>
                   <evt_file format="gst"> /path/model_1/ev2.root  </evt_file>
                   ...
                 </model>

                 <model name="another_model_name">
                   <xsec_file>             /path/model_2/xsec.root </xsec_file>
                   <evt_file format="gst"> /path/model_2/ev0.root  </evt_file>
                   <evt_file format="gst"> /path/model_2/ev1.root  </evt_file>
                   <evt_file format="gst"> /path/model_2/ev2.root  </evt_file>
                   ...
                 </model>
                 ...
              </vld_inputs>

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
#include <TH2D.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TPavesText.h>
#include <TText.h>
#include <TLegend.h>
#include <TBox.h>

#include "Conventions/Constants.h"
#include "Conventions/Units.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGUtils.h"
#include "PDG/PDGCodes.h"
#include "Utils/CmdLnArgParser.h"
#include "Utils/VldTestInputs.h"
#include "Utils/SystemUtils.h"
#include "Utils/Style.h"

using std::ostringstream;
using std::string;

using namespace genie;
using namespace genie::utils::vld;
using namespace genie::constants;

//_________________________________________________________________________________
// Utility class to hold info on plotted datasets
class eQEDataSetDescription_t
{
public:
  eQEDataSetDescription_t(
     int tgtpdg, string citation, double E, double theta, bool show) :
    fTgtPdg   (tgtpdg), 
    fCitation (citation), 
    fE        (E), 
    fTheta    (theta),
    fShow     (show)
  {   
  }
  eQEDataSetDescription_t() 
  {   
  }
  int    TgtPdg   (void) const { return fTgtPdg; }
  int    TgtZ     (void) const { return pdg::IonPdgCodeToZ(fTgtPdg); }
  int    TgtA     (void) const { return pdg::IonPdgCodeToA(fTgtPdg); }
  string TgtName  (void) const { 
    if(fTgtPdg == 1000010010) return "1H";
    if(fTgtPdg == 1000010020) return "2D";
    if(fTgtPdg == 1000010030) return "3H";
    if(fTgtPdg == 1000020030) return "3He";
    if(fTgtPdg == 1000020040) return "4He";
    if(fTgtPdg == 1000060120) return "12C";
    if(fTgtPdg == 1000080160) return "16O";
    if(fTgtPdg == 1000130270) return "27Al";
    if(fTgtPdg == 1000200400) return "40Ca";
    if(fTgtPdg == 1000200480) return "48Ca";
    if(fTgtPdg == 1000260560) return "56Fe";
    if(fTgtPdg == 1000791970) return "197Au";
    if(fTgtPdg == 1000822080) return "208Pb";
    if(fTgtPdg == 1000922380) return "238U";
    return "Other";
  }
  string Citation (void) const { return fCitation; }
  double E        (void) const { return fE; }
  double Theta    (void) const { return fTheta; }
  bool   Show     (void) const { return fShow; }
  string LabelTeX (void) const { 
    ostringstream label;
    label << this->TgtName() << "(e,e^{'}) ";
    label << "E = " << fE << " GeV, ";
    label << "#theta = " << fTheta << "^{o}";
    label << " [" << fCitation << "]";
    return label.str();
  }
private:
  int    fTgtPdg;   //
  string fCitation; //
  double fE;        //
  double fTheta;    //
  bool   fShow;     //
};
//_________________________________________________________________________________

/* 
..............................................................................
QUASIELASTIC ELECTRON-NUCLEUS DATASET DESCRIPTIONS
..............................................................................
*/
const int kNumOfDataSets = 551;

eQEDataSetDescription_t * kDataSet[kNumOfDataSets] =
{
   //2H
   new eQEDataSetDescription_t(1000010020, "Parker:1986",          0.220,  180.000,  true),
   new eQEDataSetDescription_t(1000010020, "Parker:1986",          0.269,  180.000,  true),
   new eQEDataSetDescription_t(1000010020, "Parker:1986",          0.320,  180.000,  true),
   new eQEDataSetDescription_t(1000010020, "Arnold:1988us",        0.843,  180.000,  true),
   new eQEDataSetDescription_t(1000010020, "Arnold:1988us",        1.020,  180.000,  true),
   new eQEDataSetDescription_t(1000010020, "Arnold:1988us",        1.189,  180.000,  true),
   new eQEDataSetDescription_t(1000010020, "Arnold:1988us",        1.281,  180.000,  true),
   new eQEDataSetDescription_t(1000010020, "Lung-thesis:1992",     1.511,   89.970,  true),
   new eQEDataSetDescription_t(1000010020, "Lung-thesis:1992",     1.511,   90.070,  true),
   new eQEDataSetDescription_t(1000010020, "Lung-thesis:1992",     1.968,   55.210,  true),
   new eQEDataSetDescription_t(1000010020, "Lung-thesis:1992",     1.968,   89.950,  true),
   new eQEDataSetDescription_t(1000010020, "Lung-thesis:1992",     1.968,   89.970,  true),
   new eQEDataSetDescription_t(1000010020, "Arrington:1995hs",     2.015,   38.840,  true),
   new eQEDataSetDescription_t(1000010020, "Lung-thesis:1992",     2.407,   41.110,  true),
   new eQEDataSetDescription_t(1000010020, "Lung-thesis:1992",     2.407,   58.880,  true),
   new eQEDataSetDescription_t(1000010020, "Lung-thesis:1992",     2.407,   89.970,  true),
   new eQEDataSetDescription_t(1000010020, "Lung-thesis:1992",     2.837,   44.990,  true),
   new eQEDataSetDescription_t(1000010020, "Lung-thesis:1992",     2.837,   61.210,  true),
   new eQEDataSetDescription_t(1000010020, "Lung-thesis:1992",     2.837,   89.970,  true),
   new eQEDataSetDescription_t(1000010020, "Arrington:1995hs",     3.188,   47.680,  true),
   new eQEDataSetDescription_t(1000010020, "Arrington:1998ps",     4.045,   15.000,  true),
   new eQEDataSetDescription_t(1000010020, "Arrington:1998ps",     4.045,   23.000,  true),
   new eQEDataSetDescription_t(1000010020, "Arrington:1998ps",     4.045,   30.000,  true),
   new eQEDataSetDescription_t(1000010020, "Arrington:1998ps",     4.045,   37.000,  true),
   new eQEDataSetDescription_t(1000010020, "Arrington:1998ps",     4.045,   45.000,  true),
   new eQEDataSetDescription_t(1000010020, "Arrington:1998ps",     4.045,   55.000,  true),
   new eQEDataSetDescription_t(1000010020, "Arrington:1995hs",     4.212,   53.390,  true),
   new eQEDataSetDescription_t(1000010020, "Arrington:1995hs",     5.120,   56.640,  true),
   new eQEDataSetDescription_t(1000010020, "Lung-thesis:1992",     5.507,   15.150,  true),
   new eQEDataSetDescription_t(1000010020, "Lung-thesis:1992",     5.507,   18.980,  true),
   new eQEDataSetDescription_t(1000010020, "Lung-thesis:1992",     5.507,   22.800,  true),
   new eQEDataSetDescription_t(1000010020, "Lung-thesis:1992",     5.507,   26.820,  true),
   new eQEDataSetDescription_t(1000010020, "Schutz:1976he",        6.519,    8.000,  true),
   new eQEDataSetDescription_t(1000010020, "Schutz:1976he",        7.302,    8.000,  true),
   new eQEDataSetDescription_t(1000010020, "Schutz:1976he",        8.981,    8.000,  true),
   new eQEDataSetDescription_t(1000010020, "Schutz:1976he",        9.718,    8.000,  true),
   new eQEDataSetDescription_t(1000010020, "Rock:1991jy",          9.744,   10.000,  true),
   new eQEDataSetDescription_t(1000010020, "Schutz:1976he",       10.407,    8.000,  true),
   new eQEDataSetDescription_t(1000010020, "Schutz:1976he",       11.671,    8.000,  true),
   new eQEDataSetDescription_t(1000010020, "Rock:1991jy",         12.565,   10.000,  true),
   new eQEDataSetDescription_t(1000010020, "Schutz:1976he",       12.821,    8.000,  true),
   new eQEDataSetDescription_t(1000010020, "Schutz:1976he",       14.878,    8.000,  true),
   new eQEDataSetDescription_t(1000010020, "Rock:1991jy",         15.730,   10.000,  true),
   new eQEDataSetDescription_t(1000010020, "Rock:1991jy",         17.301,   10.000,  true),
   new eQEDataSetDescription_t(1000010020, "Schutz:1976he",       18.375,    8.000,  true),
   new eQEDataSetDescription_t(1000010020, "Rock:1991jy",         18.476,   10.000,  true),
   new eQEDataSetDescription_t(1000010020, "Rock:1991jy",         20.999,   10.000,  true),

   //3H
   new eQEDataSetDescription_t(1000010030, "Dow:1988rk",           0.066,   54.000, true),
   new eQEDataSetDescription_t(1000010030, "Dow:1988rk",           0.099,  134.500, true),
   new eQEDataSetDescription_t(1000010030, "Dow:1988rk",           0.110,   54.000, true),
   new eQEDataSetDescription_t(1000010030, "Dow:1988rk",           0.145,  134.500, true),
   new eQEDataSetDescription_t(1000010030, "Dow:1988rk",           0.187,  134.500, true),
   new eQEDataSetDescription_t(1000010030, "Dow:1988rk",           0.198,   54.000, true),
   new eQEDataSetDescription_t(1000010030, "Dow:1988rk",           0.287,   54.000, true),
   new eQEDataSetDescription_t(1000010030, "Dow:1988rk",           0.289,  134.500, true),
   new eQEDataSetDescription_t(1000010030, "Dow:1988rk",           0.367,  134.500, true),
   new eQEDataSetDescription_t(1000010030, "Dow:1988rk",           0.368,   54.000, true),
   new eQEDataSetDescription_t(1000010030, "Dow:1988rk",           0.469,   54.000, true),
   new eQEDataSetDescription_t(1000010030, "Dow:1988rk",           0.507,  134.500, true),
   new eQEDataSetDescription_t(1000010030, "Dow:1988rk",           0.558,   54.000, true),
   new eQEDataSetDescription_t(1000010030, "Dow:1988rk",           0.652,   54.000, true),
   new eQEDataSetDescription_t(1000010030, "Dow:1988rk",           0.790,   54.000, true),
   new eQEDataSetDescription_t(1000010030, "Dow:1988rk",           0.790,   88.000, true),
   new eQEDataSetDescription_t(1000010030, "Dow:1988rk",           0.790,  108.000, true),

   //3He
   new eQEDataSetDescription_t(1000020030, "Dow:1988rk",           0.099,  134.500, true),
   new eQEDataSetDescription_t(1000020030, "Dow:1988rk",           0.110,   54.000, true),
   new eQEDataSetDescription_t(1000020030, "Dow:1988rk",           0.145,  134.500, true),
   new eQEDataSetDescription_t(1000020030, "Marchand:1985us",      0.169,   90.000, true),
   new eQEDataSetDescription_t(1000020030, "Dow:1988rk",           0.188,  134.500, true),
   new eQEDataSetDescription_t(1000020030, "Dow:1988rk",           0.198,   54.000, true),
   new eQEDataSetDescription_t(1000020030, "Marchand:1985us",      0.219,   36.000, true),
   new eQEDataSetDescription_t(1000020030, "Marchand:1985us",      0.219,   60.000, true),
   new eQEDataSetDescription_t(1000020030, "Marchand:1985us",      0.219,   90.000, true),
   new eQEDataSetDescription_t(1000020030, "Marchand:1985us",      0.269,   36.000, true),
   new eQEDataSetDescription_t(1000020030, "Marchand:1985us",      0.269,   60.000, true),
   new eQEDataSetDescription_t(1000020030, "Marchand:1985us",      0.269,   90.000, true),
   new eQEDataSetDescription_t(1000020030, "Marchand:1985us",      0.269,  144.500, true),
   new eQEDataSetDescription_t(1000020030, "Dow:1988rk",           0.287,   54.000, true),
   new eQEDataSetDescription_t(1000020030, "Dow:1988rk",           0.289,  134.500, true),
   new eQEDataSetDescription_t(1000020030, "Marchand:1985us",      0.319,   36.000, true),
   new eQEDataSetDescription_t(1000020030, "Marchand:1985us",      0.319,   60.000, true),
   new eQEDataSetDescription_t(1000020030, "Marchand:1985us",      0.319,   90.000, true),
   new eQEDataSetDescription_t(1000020030, "Marchand:1985us",      0.319,  144.500, true),
   new eQEDataSetDescription_t(1000020030, "Dow:1988rk",           0.367,  134.500, true),
   new eQEDataSetDescription_t(1000020030, "Dow:1988rk",           0.368,   54.000, true),
   new eQEDataSetDescription_t(1000020030, "Marchand:1985us",      0.379,   36.000, true),
   new eQEDataSetDescription_t(1000020030, "Marchand:1985us",      0.379,   60.000, true),
   new eQEDataSetDescription_t(1000020030, "Marchand:1985us",      0.379,   90.000, true),
   new eQEDataSetDescription_t(1000020030, "Marchand:1985us",      0.379,  144.500, true),
   new eQEDataSetDescription_t(1000020030, "Marchand:1985us",      0.439,   36.000, true),
   new eQEDataSetDescription_t(1000020030, "Marchand:1985us",      0.439,   60.000, true),
   new eQEDataSetDescription_t(1000020030, "Marchand:1985us",      0.439,   90.000, true),
   new eQEDataSetDescription_t(1000020030, "Marchand:1985us",      0.439,  144.500, true),
   new eQEDataSetDescription_t(1000020030, "Dow:1988rk",           0.469,   54.000, true),
   new eQEDataSetDescription_t(1000020030, "Marchand:1985us",      0.500,   36.000, true),
   new eQEDataSetDescription_t(1000020030, "Marchand:1985us",      0.500,   60.000, true),
   new eQEDataSetDescription_t(1000020030, "Marchand:1985us",      0.500,   90.000, true),
   new eQEDataSetDescription_t(1000020030, "Marchand:1985us",      0.500,  144.500, true),
   new eQEDataSetDescription_t(1000020030, "Mccarthy:1976re",      0.500,   60.000, true),
   new eQEDataSetDescription_t(1000020030, "Dow:1988rk",           0.507,  134.500, true),
   new eQEDataSetDescription_t(1000020030, "Dow:1988rk",           0.558,   54.000, true),
   new eQEDataSetDescription_t(1000020030, "Marchand:1985us",      0.560,   36.000, true),
   new eQEDataSetDescription_t(1000020030, "Marchand:1985us",      0.560,   60.000, true),
   new eQEDataSetDescription_t(1000020030, "Marchand:1985us",      0.560,   90.000, true),
   new eQEDataSetDescription_t(1000020030, "Marchand:1985us",      0.560,  144.500, true),
   new eQEDataSetDescription_t(1000020030, "Marchand:1985us",      0.620,   36.000, true),
   new eQEDataSetDescription_t(1000020030, "Marchand:1985us",      0.620,   60.000, true),
   new eQEDataSetDescription_t(1000020030, "Marchand:1985us",      0.620,   90.000, true),
   new eQEDataSetDescription_t(1000020030, "Marchand:1985us",      0.620,  144.500, true),
   new eQEDataSetDescription_t(1000020030, "Dow:1988rk",           0.652,   54.000, true),
   new eQEDataSetDescription_t(1000020030, "Marchand:1985us",      0.667,   36.000, true),
   new eQEDataSetDescription_t(1000020030, "Marchand:1985us",      0.667,   60.000, true),
   new eQEDataSetDescription_t(1000020030, "Marchand:1985us",      0.667,   90.000, true),
   new eQEDataSetDescription_t(1000020030, "Marchand:1985us",      0.667,  144.500, true),
   new eQEDataSetDescription_t(1000020030, "Dow:1988rk",           0.790,   54.000, true),
   new eQEDataSetDescription_t(1000020030, "Dow:1988rk",           0.790,   88.000, true),
   new eQEDataSetDescription_t(1000020030, "Dow:1988rk",           0.790,  108.000, true),
   new eQEDataSetDescription_t(1000020030, "Meziani:1992xr",       0.900,   85.000, true),
   new eQEDataSetDescription_t(1000020030, "Meziani:1992xr",       1.100,   85.000, true),
   new eQEDataSetDescription_t(1000020030, "Meziani:1992xr",       2.700,   15.000, true),
   new eQEDataSetDescription_t(1000020030, "Day:1979bx",           2.814,    8.000, true),
   new eQEDataSetDescription_t(1000020030, "Day:1979bx",           3.258,    8.000, true),
   new eQEDataSetDescription_t(1000020030, "Meziani:1992xr",       3.300,   15.000, true),
   new eQEDataSetDescription_t(1000020030, "Meziani:1992xr",       3.600,   15.000, true),
   new eQEDataSetDescription_t(1000020030, "Day:1979bx",           3.651,    8.000, true),
   new eQEDataSetDescription_t(1000020030, "Meziani:1992xr",       3.900,   15.000, true),
   new eQEDataSetDescription_t(1000020030, "Meziani:1992xr",       4.300,   15.000, true),
   new eQEDataSetDescription_t(1000020030, "Day:1979bx",           6.483,    8.000, true),
   new eQEDataSetDescription_t(1000020030, "Day:1979bx",           7.257,    8.000, true),
   new eQEDataSetDescription_t(1000020030, "Day:1979bx",           7.959,    8.000, true),
   new eQEDataSetDescription_t(1000020030, "Day:1979bx",           8.606,    8.000, true),
   new eQEDataSetDescription_t(1000020030, "Day:1979bx",           8.607,    8.000, true),
   new eQEDataSetDescription_t(1000020030, "Day:1979bx",           9.210,    8.000, true),
   new eQEDataSetDescription_t(1000020030, "Day:1979bx",           9.778,    8.000, true),
   new eQEDataSetDescription_t(1000020030, "Day:1979bx",          10.316,    8.000, true),
   new eQEDataSetDescription_t(1000020030, "Day:1979bx",          10.950,    8.000, true),
   new eQEDataSetDescription_t(1000020030, "Day:1979bx",          10.954,    8.000, true),
   new eQEDataSetDescription_t(1000020030, "Day:1979bx",          11.558,    8.000, true),
   new eQEDataSetDescription_t(1000020030, "Day:1979bx",          12.685,    8.000, true),
   new eQEDataSetDescription_t(1000020030, "Day:1979bx",          14.696,    8.000, true),
   new eQEDataSetDescription_t(1000020030, "Day:1979bx",          14.700,    8.000, true),

   //4He
   new eQEDataSetDescription_t(1000020040, "Zghiche:1993xg",       0.150,   60.000, true), 
   new eQEDataSetDescription_t(1000020040, "Zghiche:1993xg",       0.150,   75.000, true), 
   new eQEDataSetDescription_t(1000020040, "Zghiche:1993xg",       0.150,   90.000, true), 
   new eQEDataSetDescription_t(1000020040, "Zghiche:1993xg",       0.150,  145.000, true), 
   new eQEDataSetDescription_t(1000020040, "Zghiche:1993xg",       0.200,   60.000, true), 
   new eQEDataSetDescription_t(1000020040, "Zghiche:1993xg",       0.200,   75.000, true), 
   new eQEDataSetDescription_t(1000020040, "Zghiche:1993xg",       0.200,   90.000, true), 
   new eQEDataSetDescription_t(1000020040, "Zghiche:1993xg",       0.200,  145.000, true), 
   new eQEDataSetDescription_t(1000020040, "Zghiche:1993xg",       0.250,   34.000, true), 
   new eQEDataSetDescription_t(1000020040, "Zghiche:1993xg",       0.250,   60.000, true), 
   new eQEDataSetDescription_t(1000020040, "Zghiche:1993xg",       0.250,   75.000, true), 
   new eQEDataSetDescription_t(1000020040, "Zghiche:1993xg",       0.250,   90.000, true), 
   new eQEDataSetDescription_t(1000020040, "Zghiche:1993xg",       0.250,  145.000, true), 
   new eQEDataSetDescription_t(1000020040, "vonReden:1990ah",      0.279,  134.500, true), 
   new eQEDataSetDescription_t(1000020040, "Zghiche:1993xg",       0.300,   34.000, true), 
   new eQEDataSetDescription_t(1000020040, "Zghiche:1993xg",       0.300,   60.000, true), 
   new eQEDataSetDescription_t(1000020040, "Zghiche:1993xg",       0.300,   75.000, true), 
   new eQEDataSetDescription_t(1000020040, "Zghiche:1993xg",       0.300,   90.000, true), 
   new eQEDataSetDescription_t(1000020040, "Zghiche:1993xg",       0.300,  145.000, true), 
   new eQEDataSetDescription_t(1000020040, "vonReden:1990ah",      0.328,  134.500, true), 
   new eQEDataSetDescription_t(1000020040, "Zghiche:1993xg",       0.350,   34.000, true), 
   new eQEDataSetDescription_t(1000020040, "Zghiche:1993xg",       0.350,   60.000, true), 
   new eQEDataSetDescription_t(1000020040, "Zghiche:1993xg",       0.350,   75.000, true), 
   new eQEDataSetDescription_t(1000020040, "Zghiche:1993xg",       0.350,   90.000, true), 
   new eQEDataSetDescription_t(1000020040, "Zghiche:1993xg",       0.350,  145.000, true), 
   new eQEDataSetDescription_t(1000020040, "vonReden:1990ah",      0.364,   60.000, true), 
   new eQEDataSetDescription_t(1000020040, "vonReden:1990ah",      0.367,  134.500, true), 
   new eQEDataSetDescription_t(1000020040, "Zghiche:1993xg",       0.400,   34.000, true), 
   new eQEDataSetDescription_t(1000020040, "Zghiche:1993xg",       0.400,   60.000, true), 
   new eQEDataSetDescription_t(1000020040, "Zghiche:1993xg",       0.400,   75.000, true), 
   new eQEDataSetDescription_t(1000020040, "Zghiche:1993xg",       0.400,   90.000, true), 
   new eQEDataSetDescription_t(1000020040, "Zghiche:1993xg",       0.400,  145.000, true), 
   new eQEDataSetDescription_t(1000020040, "Zghiche:1993xg",       0.450,   34.000, true), 
   new eQEDataSetDescription_t(1000020040, "Zghiche:1993xg",       0.450,   60.000, true), 
   new eQEDataSetDescription_t(1000020040, "Zghiche:1993xg",       0.450,   75.000, true), 
   new eQEDataSetDescription_t(1000020040, "Zghiche:1993xg",       0.450,   90.000, true), 
   new eQEDataSetDescription_t(1000020040, "vonReden:1990ah",      0.476,   60.000, true), 
   new eQEDataSetDescription_t(1000020040, "Zghiche:1993xg",       0.500,   34.000, true), 
   new eQEDataSetDescription_t(1000020040, "Mccarthy:1976re",      0.500,   60.000, true), 
   new eQEDataSetDescription_t(1000020040, "Zghiche:1993xg",       0.500,   60.000, true), 
   new eQEDataSetDescription_t(1000020040, "Zghiche:1993xg",       0.500,   75.000, true), 
   new eQEDataSetDescription_t(1000020040, "Zghiche:1993xg",       0.500,  145.000, true), 
   new eQEDataSetDescription_t(1000020040, "Zghiche:1993xg",       0.550,   34.000, true), 
   new eQEDataSetDescription_t(1000020040, "Zghiche:1993xg",       0.550,   60.000, true), 
   new eQEDataSetDescription_t(1000020040, "Zghiche:1993xg",       0.550,   75.000, true), 
   new eQEDataSetDescription_t(1000020040, "Zghiche:1993xg",       0.550,  145.000, true), 
   new eQEDataSetDescription_t(1000020040, "vonReden:1990ah",      0.589,   60.000, true), 
   new eQEDataSetDescription_t(1000020040, "Zghiche:1993xg",       0.600,   34.000, true), 
   new eQEDataSetDescription_t(1000020040, "Zghiche:1993xg",       0.600,   60.000, true), 
   new eQEDataSetDescription_t(1000020040, "Zghiche:1993xg",       0.600,   75.000, true), 
   new eQEDataSetDescription_t(1000020040, "Zghiche:1993xg",       0.640,   34.000, true), 
   new eQEDataSetDescription_t(1000020040, "Zghiche:1993xg",       0.640,   60.000, true), 
   new eQEDataSetDescription_t(1000020040, "Zghiche:1993xg",       0.640,   75.000, true), 
   new eQEDataSetDescription_t(1000020040, "vonReden:1990ah",      0.724,   60.000, true), 
   new eQEDataSetDescription_t(1000020040, "O'Connell:1987ag",     0.730,   37.100, true),
   new eQEDataSetDescription_t(1000020040, "Meziani:1992xr",       0.900,   85.000, true), 
   new eQEDataSetDescription_t(1000020040, "Sealock:1989nx",       0.961,   37.500, true), 
   new eQEDataSetDescription_t(1000020040, "Meziani:1992xr",       1.100,   85.000, true), 
   new eQEDataSetDescription_t(1000020040, "Sealock:1989nx",       1.108,   37.500, true), 
   new eQEDataSetDescription_t(1000020040, "Sealock:1989nx",       1.299,   37.500, true), 
   new eQEDataSetDescription_t(1000020040, "Day:1993md",           2.020,   15.022, true), 
   new eQEDataSetDescription_t(1000020040, "Day:1993md",           2.020,   20.016, true), 
   new eQEDataSetDescription_t(1000020040, "Meziani:1992xr",       2.700,   15.000, true), 
   new eQEDataSetDescription_t(1000020040, "Meziani:1992xr",       3.300,   15.000, true), 
   new eQEDataSetDescription_t(1000020040, "Day:1993md",           3.595,   16.020, true), 
   new eQEDataSetDescription_t(1000020040, "Day:1993md",           3.595,   20.016, true), 
   new eQEDataSetDescription_t(1000020040, "Day:1993md",           3.595,   25.012, true), 
   new eQEDataSetDescription_t(1000020040, "Day:1993md",           3.595,   30.010, true), 
   new eQEDataSetDescription_t(1000020040, "Meziani:1992xr",       3.600,   15.000, true), 
   new eQEDataSetDescription_t(1000020040, "Meziani:1992xr",       3.900,   15.000, true), 
   new eQEDataSetDescription_t(1000020040, "Meziani:1992xr",       4.300,   15.000, true), 
   new eQEDataSetDescription_t(1000020040, "Rock:1981aa",          6.465,    8.000, true), 
   new eQEDataSetDescription_t(1000020040, "Rock:1981aa",          7.235,    8.000, true), 
   new eQEDataSetDescription_t(1000020040, "Rock:1981aa",          7.933,    8.000, true), 
   new eQEDataSetDescription_t(1000020040, "Rock:1981aa",          8.576,    8.000, true), 
   new eQEDataSetDescription_t(1000020040, "Rock:1981aa",          9.175,    8.000, true), 
   new eQEDataSetDescription_t(1000020040, "Rock:1981aa",          9.738,    8.000, true), 
   new eQEDataSetDescription_t(1000020040, "Rock:1981aa",         11.267,    8.000, true),

   //12C
   new eQEDataSetDescription_t(1000060120, "Barreau:1983ht",       0.160,   36.000, true), 
   new eQEDataSetDescription_t(1000060120, "Barreau:1983ht",       0.161,   60.000, true), 
   new eQEDataSetDescription_t(1000060120, "Barreau:1983ht",       0.200,   36.000, true), 
   new eQEDataSetDescription_t(1000060120, "Barreau:1983ht",       0.200,   60.000, true), 
   new eQEDataSetDescription_t(1000060120, "Barreau:1983ht",       0.240,   36.000, true), 
   new eQEDataSetDescription_t(1000060120, "Barreau:1983ht",       0.240,   60.000, true), 
   new eQEDataSetDescription_t(1000060120, "Barreau:1983ht",       0.280,   60.000, true), 
   new eQEDataSetDescription_t(1000060120, "Barreau:1983ht",       0.320,   36.000, true), 
   new eQEDataSetDescription_t(1000060120, "Barreau:1983ht",       0.320,   60.000, true), 
   new eQEDataSetDescription_t(1000060120, "Barreau:1983ht",       0.361,   60.000, true), 
   new eQEDataSetDescription_t(1000060120, "Barreau:1983ht",       0.400,   36.000, true), 
   new eQEDataSetDescription_t(1000060120, "Barreau:1983ht",       0.401,   60.000, true), 
   new eQEDataSetDescription_t(1000060120, "Barreau:1983ht",       0.440,   60.000, true), 
   new eQEDataSetDescription_t(1000060120, "Barreau:1983ht",       0.480,   36.000, true), 
   new eQEDataSetDescription_t(1000060120, "Barreau:1983ht",       0.480,   60.000, true), 
   new eQEDataSetDescription_t(1000060120, "Whitney:1974hr",       0.500,   60.000, true), 
   new eQEDataSetDescription_t(1000060120, "Barreau:1983ht",       0.519,   60.000, true), 
   new eQEDataSetDescription_t(1000060120, "Barreau:1983ht",       0.560,   36.000, true), 
   new eQEDataSetDescription_t(1000060120, "Barreau:1983ht",       0.560,   60.000, true), 
   new eQEDataSetDescription_t(1000060120, "Barreau:1983ht",       0.560,  145.000, true), 
   new eQEDataSetDescription_t(1000060120, "Barreau:1983ht",       0.620,   36.000, true), 
   new eQEDataSetDescription_t(1000060120, "Barreau:1983ht",       0.620,   60.000, true), 
   new eQEDataSetDescription_t(1000060120, "Barreau:1983ht",       0.680,   36.000, true), 
   new eQEDataSetDescription_t(1000060120, "Barreau:1983ht",       0.680,   60.000, true), 
   new eQEDataSetDescription_t(1000060120, "O'Connell:1987ag",     0.730,   37.100, true), 
   new eQEDataSetDescription_t(1000060120, "Sealock:1989nx",       0.961,   37.500, true), 
   new eQEDataSetDescription_t(1000060120, "Sealock:1989nx",       1.108,   37.500, true), 
   new eQEDataSetDescription_t(1000060120, "Sealock:1989nx",       1.299,   37.500, true), 
   new eQEDataSetDescription_t(1000060120, "Baran:1988tw",         1.300,   11.950, true), 
   new eQEDataSetDescription_t(1000060120, "Baran:1988tw",         1.300,   13.540, true), 
   new eQEDataSetDescription_t(1000060120, "Baran:1988tw",         1.500,   11.950, true), 
   new eQEDataSetDescription_t(1000060120, "Baran:1988tw",         1.500,   13.540, true), 
   new eQEDataSetDescription_t(1000060120, "Sealock:1989nx",       1.501,   37.500, true), 
   new eQEDataSetDescription_t(1000060120, "Baran:1988tw",         1.650,   11.950, true), 
   new eQEDataSetDescription_t(1000060120, "Baran:1988tw",         1.650,   13.540, true), 
   new eQEDataSetDescription_t(1000060120, "Bagdasaryan:1988hp",   1.930,   16.000, true), 
   new eQEDataSetDescription_t(1000060120, "Bagdasaryan:1988hp",   1.930,   18.000, true), 
   new eQEDataSetDescription_t(1000060120, "Zeller:1973ge",        2.000,   15.000, true), 
   new eQEDataSetDescription_t(1000060120, "Arrington:1995hs",     2.015,   35.510, true), 
   new eQEDataSetDescription_t(1000060120, "Day:1993md",           2.020,   15.022, true), 
   new eQEDataSetDescription_t(1000060120, "Day:1993md",           2.020,   20.016, true), 
   new eQEDataSetDescription_t(1000060120, "Bagdasaryan:1988hp",   2.130,   16.000, true), 
   new eQEDataSetDescription_t(1000060120, "Bagdasaryan:1988hp",   2.130,   18.000, true), 
   new eQEDataSetDescription_t(1000060120, "Zeller:1973ge",        2.500,   15.000, true), 
   new eQEDataSetDescription_t(1000060120, "Zeller:1973ge",        2.700,   15.000, true), 
   new eQEDataSetDescription_t(1000060120, "Arrington:1995hs",     3.188,   47.680, true), 
   new eQEDataSetDescription_t(1000060120, "Day:1993md",           3.595,   16.020, true), 
   new eQEDataSetDescription_t(1000060120, "Day:1993md",           3.595,   20.016, true), 
   new eQEDataSetDescription_t(1000060120, "Day:1993md",           3.595,   25.012, true), 
   new eQEDataSetDescription_t(1000060120, "Day:1993md",           3.595,   30.010, true), 
   new eQEDataSetDescription_t(1000060120, "Day:1993md",           3.605,   16.020, true), 
   new eQEDataSetDescription_t(1000060120, "Arrington:1998ps",     4.045,   15.000, true), 
   new eQEDataSetDescription_t(1000060120, "Arrington:1998ps",     4.045,   23.000, true), 
   new eQEDataSetDescription_t(1000060120, "Arrington:1998ps",     4.045,   30.000, true), 
   new eQEDataSetDescription_t(1000060120, "Arrington:1998ps",     4.045,   37.000, true), 
   new eQEDataSetDescription_t(1000060120, "Arrington:1998ps",     4.045,   45.000, true), 
   new eQEDataSetDescription_t(1000060120, "Arrington:1998ps",     4.045,   55.000, true), 
   new eQEDataSetDescription_t(1000060120, "Arrington:1998ps",     4.045,   74.000, true), 
   new eQEDataSetDescription_t(1000060120, "Arrington:1995hs",     4.212,   53.389, true), 
   new eQEDataSetDescription_t(1000060120, "Arrington:1995hs",     5.120,   56.639, true), 

   //16O
   new eQEDataSetDescription_t(1000080160, "Anghinolfi:1996vm",    0.700,   32.000, true),
   new eQEDataSetDescription_t(1000080160, "O\'Connell:1987ag",    0.737,   37.100, true),
   new eQEDataSetDescription_t(1000080160, "Anghinolfi:1996vm",    0.880,   32.000, true),
   new eQEDataSetDescription_t(1000080160, "Anghinolfi:1996vm",    1.080,   32.000, true),
   new eQEDataSetDescription_t(1000080160, "Anghinolfi:1996vm",    1.200,   32.000, true),
   new eQEDataSetDescription_t(1000080160, "Anghinolfi:1996vm",    1.500,   32.000, true),

   //27Al
   new eQEDataSetDescription_t(1000130270, "Bosted:1992fy",        1.968,   55.210, true), 
   new eQEDataSetDescription_t(1000130270, "Day:1993md",           2.020,   15.022, true), 
   new eQEDataSetDescription_t(1000130270, "Day:1993md",           2.020,   20.016, true), 
   new eQEDataSetDescription_t(1000130270, "Bosted:1992fy",        2.407,   41.110, true), 
   new eQEDataSetDescription_t(1000130270, "Bosted:1992fy",        2.407,   58.880, true), 
   new eQEDataSetDescription_t(1000130270, "Rock-pc",              2.814,    8.000, true), 
   new eQEDataSetDescription_t(1000130270, "Bosted:1992fy",        2.837,   44.990, true), 
   new eQEDataSetDescription_t(1000130270, "Bosted:1992fy",        2.837,   61.200, true), 
   new eQEDataSetDescription_t(1000130270, "Rock-pc",              3.258,    8.000, true), 
   new eQEDataSetDescription_t(1000130270, "Bosted:1992fy",        3.400,   34.690, true), 
   new eQEDataSetDescription_t(1000130270, "Bosted:1992fy",        3.400,   44.480, true), 
   new eQEDataSetDescription_t(1000130270, "Bosted:1992fy",        3.400,   57.570, true), 
   new eQEDataSetDescription_t(1000130270, "Day:1993md",           3.595,   16.020, true), 
   new eQEDataSetDescription_t(1000130270, "Day:1993md",           3.595,   20.016, true), 
   new eQEDataSetDescription_t(1000130270, "Day:1993md",           3.595,   25.012, true), 
   new eQEDataSetDescription_t(1000130270, "Day:1993md",           3.595,   30.010, true), 
   new eQEDataSetDescription_t(1000130270, "Rock-pc",              3.651,    8.000, true), 
   new eQEDataSetDescription_t(1000130270, "Bosted:1992fy",        3.956,   28.410, true), 
   new eQEDataSetDescription_t(1000130270, "Bosted:1992fy",        3.956,   35.380, true), 
   new eQEDataSetDescription_t(1000130270, "Bosted:1992fy",        3.956,   43.710, true), 
   new eQEDataSetDescription_t(1000130270, "Bosted:1992fy",        3.956,   59.290, true), 
   new eQEDataSetDescription_t(1000130270, "Bosted:1992fy",        4.507,   35.590, true), 
   new eQEDataSetDescription_t(1000130270, "Bosted:1992fy",        4.507,   45.660, true), 
   new eQEDataSetDescription_t(1000130270, "Bosted:1992fy",        5.507,   15.140, true), 
   new eQEDataSetDescription_t(1000130270, "Bosted:1992fy",        5.507,   18.980, true), 
   new eQEDataSetDescription_t(1000130270, "Bosted:1992fy",        5.507,   22.800, true), 
   new eQEDataSetDescription_t(1000130270, "Bosted:1992fy",        5.507,   26.820, true), 
   new eQEDataSetDescription_t(1000130270, "Bosted:1992fy",        5.507,   32.830, true), 
   new eQEDataSetDescription_t(1000130270, "Rock-pc",              7.257,    8.000, true), 
   new eQEDataSetDescription_t(1000130270, "Bosted:1992fy",        9.800,   13.250, true), 
   new eQEDataSetDescription_t(1000130270, "Bosted:1992fy",        9.800,   15.370, true), 
   new eQEDataSetDescription_t(1000130270, "Bosted:1992fy",        9.800,   17.520, true), 
   new eQEDataSetDescription_t(1000130270, "Bosted:1992fy",        9.800,   19.750, true), 
   new eQEDataSetDescription_t(1000130270, "Rock-pc",             10.950,    8.000, true), 
   new eQEDataSetDescription_t(1000130270, "Rock-pc",             14.700,    8.000, true), 

   //40Ca
   new eQEDataSetDescription_t(1000200400, "Meziani:1984is",       0.120,   90.000, true), 
   new eQEDataSetDescription_t(1000200400, "Meziani:1984is",       0.120,  140.000, true), 
   new eQEDataSetDescription_t(1000200400, "Williamson:1997",      0.130,  140.000, true), 
   new eQEDataSetDescription_t(1000200400, "Williamson:1997",      0.150,   90.000, true), 
   new eQEDataSetDescription_t(1000200400, "Meziani:1984is",       0.160,   60.000, true), 
   new eQEDataSetDescription_t(1000200400, "Meziani:1984is",       0.160,   90.000, true), 
   new eQEDataSetDescription_t(1000200400, "Meziani:1984is",       0.160,  140.000, true), 
   new eQEDataSetDescription_t(1000200400, "Williamson:1997",      0.160,  140.000, true), 
   new eQEDataSetDescription_t(1000200400, "Williamson:1997",      0.189,  140.000, true), 
   new eQEDataSetDescription_t(1000200400, "Meziani:1984is",       0.200,   60.000, true), 
   new eQEDataSetDescription_t(1000200400, "Meziani:1984is",       0.200,   90.000, true), 
   new eQEDataSetDescription_t(1000200400, "Williamson:1997",      0.200,   90.000, true), 
   new eQEDataSetDescription_t(1000200400, "Meziani:1984is",       0.200,  140.000, true), 
   new eQEDataSetDescription_t(1000200400, "Williamson:1997",      0.219,  140.000, true), 
   new eQEDataSetDescription_t(1000200400, "Meziani:1984is",       0.240,   60.000, true), 
   new eQEDataSetDescription_t(1000200400, "Meziani:1984is",       0.240,   90.000, true), 
   new eQEDataSetDescription_t(1000200400, "Meziani:1984is",       0.240,  140.000, true), 
   new eQEDataSetDescription_t(1000200400, "Williamson:1997",      0.248,   90.000, true), 
   new eQEDataSetDescription_t(1000200400, "Williamson:1997",      0.249,  140.000, true), 
   new eQEDataSetDescription_t(1000200400, "Meziani:1984is",       0.280,   60.000, true), 
   new eQEDataSetDescription_t(1000200400, "Meziani:1984is",       0.280,   90.000, true), 
   new eQEDataSetDescription_t(1000200400, "Meziani:1984is",       0.280,  140.000, true), 
   new eQEDataSetDescription_t(1000200400, "Williamson:1997",      0.288,  140.000, true), 
   new eQEDataSetDescription_t(1000200400, "Williamson:1997",      0.298,   90.000, true), 
   new eQEDataSetDescription_t(1000200400, "Meziani:1984is",       0.320,   90.000, true), 
   new eQEDataSetDescription_t(1000200400, "Meziani:1984is",       0.320,  140.000, true), 
   new eQEDataSetDescription_t(1000200400, "Williamson:1997",      0.327,  140.000, true), 
   new eQEDataSetDescription_t(1000200400, "Meziani:1984is",       0.340,   60.000, true), 
   new eQEDataSetDescription_t(1000200400, "Williamson:1997",      0.347,   90.000, true), 
   new eQEDataSetDescription_t(1000200400, "Williamson:1997",      0.348,   45.500, true), 
   new eQEDataSetDescription_t(1000200400, "Meziani:1984is",       0.360,   90.000, true), 
   new eQEDataSetDescription_t(1000200400, "Meziani:1984is",       0.360,  140.000, true), 
   new eQEDataSetDescription_t(1000200400, "Williamson:1997",      0.367,  140.000, true), 
   new eQEDataSetDescription_t(1000200400, "Williamson:1997",      0.372,   90.000, true), 
   new eQEDataSetDescription_t(1000200400, "Meziani:1984is",       0.400,   60.000, true), 
   new eQEDataSetDescription_t(1000200400, "Meziani:1984is",       0.400,   90.000, true), 
   new eQEDataSetDescription_t(1000200400, "Meziani:1984is",       0.400,  140.000, true), 
   new eQEDataSetDescription_t(1000200400, "Williamson:1997",      0.408,   45.500, true), 
   new eQEDataSetDescription_t(1000200400, "Meziani:1984is",       0.440,  140.000, true), 
   new eQEDataSetDescription_t(1000200400, "Williamson:1997",      0.471,   45.500, true), 
   new eQEDataSetDescription_t(1000200400, "Meziani:1984is",       0.480,   60.000, true), 
   new eQEDataSetDescription_t(1000200400, "Meziani:1984is",       0.480,  140.000, true), 
   new eQEDataSetDescription_t(1000200400, "Whitney:1974hr",       0.500,   60.000, true), 
   new eQEDataSetDescription_t(1000200400, "Williamson:1997",      0.545,   45.500, true), 
   new eQEDataSetDescription_t(1000200400, "Meziani:1984is",       0.560,   60.000, true), 
   new eQEDataSetDescription_t(1000200400, "Williamson:1997",      0.628,   45.500, true), 
   new eQEDataSetDescription_t(1000200400, "Williamson:1997",      0.681,   45.500, true), 
   new eQEDataSetDescription_t(1000200400, "Williamson:1997",      0.739,   45.500, true), 
   new eQEDataSetDescription_t(1000200400, "Williamson:1997",      0.841,   45.500, true),

   //48Ca
   new eQEDataSetDescription_t(1000200480, "Meziani:1984is",       0.120,   90.000, true), 
   new eQEDataSetDescription_t(1000200480, "Meziani:1984is",       0.120,  140.000, true), 
   new eQEDataSetDescription_t(1000200480, "Meziani:1984is",       0.160,   60.000, true), 
   new eQEDataSetDescription_t(1000200480, "Meziani:1984is",       0.160,   90.000, true), 
   new eQEDataSetDescription_t(1000200480, "Meziani:1984is",       0.160,  140.000, true), 
   new eQEDataSetDescription_t(1000200480, "Meziani:1984is",       0.200,   60.000, true), 
   new eQEDataSetDescription_t(1000200480, "Meziani:1984is",       0.200,   90.000, true), 
   new eQEDataSetDescription_t(1000200480, "Meziani:1984is",       0.200,  140.000, true), 
   new eQEDataSetDescription_t(1000200480, "Meziani:1984is",       0.240,   60.000, true), 
   new eQEDataSetDescription_t(1000200480, "Meziani:1984is",       0.240,   90.000, true), 
   new eQEDataSetDescription_t(1000200480, "Meziani:1984is",       0.240,  140.000, true), 
   new eQEDataSetDescription_t(1000200480, "Meziani:1984is",       0.280,   60.000, true), 
   new eQEDataSetDescription_t(1000200480, "Meziani:1984is",       0.280,   90.000, true), 
   new eQEDataSetDescription_t(1000200480, "Meziani:1984is",       0.280,  140.000, true), 
   new eQEDataSetDescription_t(1000200480, "Meziani:1984is",       0.320,   90.000, true), 
   new eQEDataSetDescription_t(1000200480, "Meziani:1984is",       0.320,  140.000, true), 
   new eQEDataSetDescription_t(1000200480, "Meziani:1984is",       0.340,   60.000, true), 
   new eQEDataSetDescription_t(1000200480, "Meziani:1984is",       0.360,   90.000, true), 
   new eQEDataSetDescription_t(1000200480, "Meziani:1984is",       0.360,  140.000, true), 
   new eQEDataSetDescription_t(1000200480, "Meziani:1984is",       0.400,   60.000, true), 
   new eQEDataSetDescription_t(1000200480, "Meziani:1984is",       0.400,   90.000, true), 
   new eQEDataSetDescription_t(1000200480, "Meziani:1984is",       0.400,  140.000, true), 
   new eQEDataSetDescription_t(1000200480, "Meziani:1984is",       0.440,  140.000, true), 
   new eQEDataSetDescription_t(1000200480, "Meziani:1984is",       0.480,   60.000, true), 
   new eQEDataSetDescription_t(1000200480, "Meziani:1984is",       0.480,  140.000, true), 
   new eQEDataSetDescription_t(1000200480, "Meziani:1984is",       0.520,   90.000, true), 
   new eQEDataSetDescription_t(1000200480, "Meziani:1984is",       0.560,   60.000, true), 
   new eQEDataSetDescription_t(1000200480, "Meziani:1984is",       0.560,   90.000, true), 
   new eQEDataSetDescription_t(1000200480, "Meziani:1984is",       0.620,   60.000, true),

   //56Fe
   new eQEDataSetDescription_t(1000260560, "Meziani:1984is",       0.120,   90.000, true), 
   new eQEDataSetDescription_t(1000260560, "Meziani:1984is",       0.120,  140.000, true), 
   new eQEDataSetDescription_t(1000260560, "Hotta:1994",           0.155,  180.000, true), 
   new eQEDataSetDescription_t(1000260560, "Meziani:1984is",       0.160,   60.000, true), 
   new eQEDataSetDescription_t(1000260560, "Meziani:1984is",       0.160,  140.000, true), 
   new eQEDataSetDescription_t(1000260560, "Hotta:1994",           0.174,  180.000, true), 
   new eQEDataSetDescription_t(1000260560, "Hotta:1994",           0.192,  180.000, true), 
   new eQEDataSetDescription_t(1000260560, "Meziani:1984is",       0.200,   60.000, true), 
   new eQEDataSetDescription_t(1000260560, "Meziani:1984is",       0.200,   90.000, true), 
   new eQEDataSetDescription_t(1000260560, "Meziani:1984is",       0.200,  140.000, true), 
   new eQEDataSetDescription_t(1000260560, "Hotta:1994",           0.206,  180.000, true), 
   new eQEDataSetDescription_t(1000260560, "Hotta:1994",           0.221,  180.000, true), 
   new eQEDataSetDescription_t(1000260560, "Hotta:1994",           0.236,  180.000, true), 
   new eQEDataSetDescription_t(1000260560, "Meziani:1984is",       0.240,   60.000, true), 
   new eQEDataSetDescription_t(1000260560, "Meziani:1984is",       0.240,  140.000, true), 
   new eQEDataSetDescription_t(1000260560, "Hotta:1994",           0.252,  180.000, true), 
   new eQEDataSetDescription_t(1000260560, "Hotta:1994",           0.271,  180.000, true), 
   new eQEDataSetDescription_t(1000260560, "Meziani:1984is",       0.280,   60.000, true), 
   new eQEDataSetDescription_t(1000260560, "Meziani:1984is",       0.280,  140.000, true), 
   new eQEDataSetDescription_t(1000260560, "Hotta:1994",           0.293,  180.000, true), 
   new eQEDataSetDescription_t(1000260560, "Hotta:1994",           0.312,  180.000, true), 
   new eQEDataSetDescription_t(1000260560, "Meziani:1984is",       0.320,   90.000, true), 
   new eQEDataSetDescription_t(1000260560, "Meziani:1984is",       0.320,  140.000, true), 
   new eQEDataSetDescription_t(1000260560, "Hotta:1994",           0.333,  180.000, true), 
   new eQEDataSetDescription_t(1000260560, "Meziani:1984is",       0.340,   60.000, true), 
   new eQEDataSetDescription_t(1000260560, "Meziani:1984is",       0.360,   90.000, true), 
   new eQEDataSetDescription_t(1000260560, "Meziani:1984is",       0.360,  140.000, true), 
   new eQEDataSetDescription_t(1000260560, "Meziani:1984is",       0.400,   60.000, true), 
   new eQEDataSetDescription_t(1000260560, "Meziani:1984is",       0.400,   90.000, true), 
   new eQEDataSetDescription_t(1000260560, "Meziani:1984is",       0.400,  140.000, true), 
   new eQEDataSetDescription_t(1000260560, "Meziani:1984is",       0.440,   90.000, true), 
   new eQEDataSetDescription_t(1000260560, "Meziani:1984is",       0.440,  140.000, true), 
   new eQEDataSetDescription_t(1000260560, "Meziani:1984is",       0.480,   60.000, true), 
   new eQEDataSetDescription_t(1000260560, "Meziani:1984is",       0.480,  140.000, true), 
   new eQEDataSetDescription_t(1000260560, "Meziani:1984is",       0.560,   60.000, true), 
   new eQEDataSetDescription_t(1000260560, "Meziani:1984is",       0.560,   90.000, true), 
   new eQEDataSetDescription_t(1000260560, "Meziani:1984is",       0.620,   60.000, true), 
   new eQEDataSetDescription_t(1000260560, "Chen:1990kq",          0.900,   85.000, true), 
   new eQEDataSetDescription_t(1000260560, "Sealock:1989nx",       0.961,   37.500, true), 
   new eQEDataSetDescription_t(1000260560, "Chen:1990kq",          1.100,   85.000, true), 
   new eQEDataSetDescription_t(1000260560, "Sealock:1989nx",       1.108,   37.500, true), 
   new eQEDataSetDescription_t(1000260560, "Chen:1990kq",          1.250,   85.000, true), 
   new eQEDataSetDescription_t(1000260560, "Sealock:1989nx",       1.299,   37.500, true), 
   new eQEDataSetDescription_t(1000260560, "Baran:1988tw",         1.500,   11.939, true), 
   new eQEDataSetDescription_t(1000260560, "Baran:1988tw",         1.500,   13.539, true), 
   new eQEDataSetDescription_t(1000260560, "Baran:1988tw",         1.650,   11.939, true), 
   new eQEDataSetDescription_t(1000260560, "Baran:1988tw",         1.650,   13.539, true), 
   new eQEDataSetDescription_t(1000260560, "Arrington:1995hs",     2.015,   35.510, true), 
   new eQEDataSetDescription_t(1000260560, "Day:1993md",           2.020,   15.022, true), 
   new eQEDataSetDescription_t(1000260560, "Day:1993md",           2.020,   20.016, true), 
   new eQEDataSetDescription_t(1000260560, "Chen:1990kq",          2.700,   15.000, true), 
   new eQEDataSetDescription_t(1000260560, "Arrington:1995hs",     3.188,   47.680, true), 
   new eQEDataSetDescription_t(1000260560, "Chen:1990kq",          3.300,   15.000, true), 
   new eQEDataSetDescription_t(1000260560, "Day:1993md",           3.595,   16.020, true), 
   new eQEDataSetDescription_t(1000260560, "Day:1993md",           3.595,   20.016, true), 
   new eQEDataSetDescription_t(1000260560, "Day:1993md",           3.595,   25.012, true), 
   new eQEDataSetDescription_t(1000260560, "Day:1993md",           3.595,   30.010, true), 
   new eQEDataSetDescription_t(1000260560, "Day:1993md",           3.595,   39.007, true), 
   new eQEDataSetDescription_t(1000260560, "Chen:1990kq",          3.600,   15.000, true), 
   new eQEDataSetDescription_t(1000260560, "Day:1993md",           3.605,   16.020, true), 
   new eQEDataSetDescription_t(1000260560, "Chen:1990kq",          3.900,   15.000, true), 
   new eQEDataSetDescription_t(1000260560, "Day:1993md",           3.995,   30.050, true), 
   new eQEDataSetDescription_t(1000260560, "Arrington:1998ps",     4.045,   15.000, true), 
   new eQEDataSetDescription_t(1000260560, "Arrington:1998ps",     4.045,   23.000, true), 
   new eQEDataSetDescription_t(1000260560, "Arrington:1998ps",     4.045,   30.000, true), 
   new eQEDataSetDescription_t(1000260560, "Arrington:1998ps",     4.045,   37.000, true), 
   new eQEDataSetDescription_t(1000260560, "Arrington:1998ps",     4.045,   45.000, true), 
   new eQEDataSetDescription_t(1000260560, "Arrington:1998ps",     4.045,   55.000, true), 
   new eQEDataSetDescription_t(1000260560, "Arrington:1998ps",     4.045,   74.000, true), 
   new eQEDataSetDescription_t(1000260560, "Arrington:1995hs",     4.212,   53.390, true), 
   new eQEDataSetDescription_t(1000260560, "Chen:1990kq",          4.300,   15.000, true), 
   new eQEDataSetDescription_t(1000260560, "Arrington:1995hs",     5.120,   56.640, true), 
 
   //197Au
   new eQEDataSetDescription_t(1000791970, "Arrington:1995hs",     2.015,   35.510, true), 
   new eQEDataSetDescription_t(1000791970, "Day:1993md",           2.020,   15.022, true), 
   new eQEDataSetDescription_t(1000791970, "Day:1993md",           2.020,   20.016, true), 
   new eQEDataSetDescription_t(1000791970, "Arrington:1995hs",     3.188,   47.680, true), 
   new eQEDataSetDescription_t(1000791970, "Day:1993md",           3.595,   16.020, true), 
   new eQEDataSetDescription_t(1000791970, "Day:1993md",           3.595,   20.016, true), 
   new eQEDataSetDescription_t(1000791970, "Day:1993md",           3.595,   25.012, true), 
   new eQEDataSetDescription_t(1000791970, "Day:1993md",           3.595,   30.010, true), 
   new eQEDataSetDescription_t(1000791970, "Day:1993md",           3.605,   16.020, true), 
   new eQEDataSetDescription_t(1000791970, "Arrington:1998ps",     4.045,   15.000, true), 
   new eQEDataSetDescription_t(1000791970, "Arrington:1998ps",     4.045,   23.000, true), 
   new eQEDataSetDescription_t(1000791970, "Arrington:1998ps",     4.045,   30.000, true), 
   new eQEDataSetDescription_t(1000791970, "Arrington:1998ps",     4.045,   45.000, true), 
   new eQEDataSetDescription_t(1000791970, "Arrington:1998ps",     4.045,   55.000, true), 
   new eQEDataSetDescription_t(1000791970, "Arrington:1998ps",     4.045,   74.000, true), 
   new eQEDataSetDescription_t(1000791970, "Arrington:1995hs",     4.212,   53.390, true), 
   new eQEDataSetDescription_t(1000791970, "Arrington:1995hs",     5.120,   56.640, true),

   //208Pb 
   new eQEDataSetDescription_t(1000822080, "Zghiche:1993xg",       0.140,   60.000, true), 
   new eQEDataSetDescription_t(1000822080, "Zghiche:1993xg",       0.140,   75.000, true), 
   new eQEDataSetDescription_t(1000822080, "Zghiche:1993xg",       0.140,   90.000, true), 
   new eQEDataSetDescription_t(1000822080, "Zghiche:1993xg",       0.140,  143.000, true), 
   new eQEDataSetDescription_t(1000822080, "Zghiche:1993xg",       0.206,   60.000, true), 
   new eQEDataSetDescription_t(1000822080, "Zghiche:1993xg",       0.206,   75.000, true), 
   new eQEDataSetDescription_t(1000822080, "Zghiche:1993xg",       0.206,  143.000, true), 
   new eQEDataSetDescription_t(1000822080, "Zghiche:1993xg",       0.262,   35.000, true), 
   new eQEDataSetDescription_t(1000822080, "Zghiche:1993xg",       0.262,   60.000, true), 
   new eQEDataSetDescription_t(1000822080, "Zghiche:1993xg",       0.262,   90.000, true), 
   new eQEDataSetDescription_t(1000822080, "Zghiche:1993xg",       0.262,  143.000, true), 
   new eQEDataSetDescription_t(1000822080, "Zghiche:1993xg",       0.310,   35.000, true), 
   new eQEDataSetDescription_t(1000822080, "Zghiche:1993xg",       0.310,   60.000, true), 
   new eQEDataSetDescription_t(1000822080, "Zghiche:1993xg",       0.310,   75.000, true), 
   new eQEDataSetDescription_t(1000822080, "Zghiche:1993xg",       0.310,  143.000, true), 
   new eQEDataSetDescription_t(1000822080, "Zghiche:1993xg",       0.354,   35.000, true), 
   new eQEDataSetDescription_t(1000822080, "Zghiche:1993xg",       0.354,   60.000, true), 
   new eQEDataSetDescription_t(1000822080, "Zghiche:1993xg",       0.354,   90.000, true), 
   new eQEDataSetDescription_t(1000822080, "Zghiche:1993xg",       0.354,  143.000, true), 
   new eQEDataSetDescription_t(1000822080, "Zghiche:1993xg",       0.420,   35.000, true), 
   new eQEDataSetDescription_t(1000822080, "Zghiche:1993xg",       0.420,   60.000, true), 
   new eQEDataSetDescription_t(1000822080, "Zghiche:1993xg",       0.420,   75.000, true), 
   new eQEDataSetDescription_t(1000822080, "Zghiche:1993xg",       0.420,  143.000, true), 
   new eQEDataSetDescription_t(1000822080, "Zghiche:1993xg",       0.485,   35.000, true), 
   new eQEDataSetDescription_t(1000822080, "Zghiche:1993xg",       0.485,   60.000, true), 
   new eQEDataSetDescription_t(1000822080, "Zghiche:1993xg",       0.485,   90.000, true), 
   new eQEDataSetDescription_t(1000822080, "Zghiche:1993xg",       0.485,  143.000, true), 
   new eQEDataSetDescription_t(1000822080, "Whitney:1974hr",       0.500,   60.000, true), 
   new eQEDataSetDescription_t(1000822080, "Zghiche:1993xg",       0.550,   35.000, true), 
   new eQEDataSetDescription_t(1000822080, "Zghiche:1993xg",       0.550,   60.000, true), 
   new eQEDataSetDescription_t(1000822080, "Zghiche:1993xg",       0.550,   75.000, true), 
   new eQEDataSetDescription_t(1000822080, "Zghiche:1993xg",       0.550,  143.000, true), 
   new eQEDataSetDescription_t(1000822080, "Zghiche:1993xg",       0.600,   60.000, true), 
   new eQEDataSetDescription_t(1000822080, "Zghiche:1993xg",       0.600,   75.000, true), 
   new eQEDataSetDescription_t(1000822080, "Zghiche:1993xg",       0.645,   35.000, true), 
   new eQEDataSetDescription_t(1000822080, "Zghiche:1993xg",       0.645,   60.000, true), 
   new eQEDataSetDescription_t(1000822080, "Zghiche:1993xg",       0.645,   75.000, true),

   //238U
   new eQEDataSetDescription_t(1000922380, "Blatchley:1986qd",     0.100,  160.000, true), 
   new eQEDataSetDescription_t(1000922380, "Blatchley:1986qd",     0.100,  140.000, true), 
   new eQEDataSetDescription_t(1000922380, "Blatchley:1986qd",     0.149,  160.000, true), 
   new eQEDataSetDescription_t(1000922380, "Blatchley:1986qd",     0.157,   90.000, true), 
   new eQEDataSetDescription_t(1000922380, "Blatchley:1986qd",     0.159,  140.000, true), 
   new eQEDataSetDescription_t(1000922380, "Blatchley:1986qd",     0.191,  160.000, true), 
   new eQEDataSetDescription_t(1000922380, "Blatchley:1986qd",     0.196,   90.000, true), 
   new eQEDataSetDescription_t(1000922380, "Blatchley:1986qd",     0.199,  134.500, true), 
   new eQEDataSetDescription_t(1000922380, "Blatchley:1986qd",     0.243,  134.500, true), 
   new eQEDataSetDescription_t(1000922380, "Blatchley:1986qd",     0.244,   90.000, true), 
   new eQEDataSetDescription_t(1000922380, "Blatchley:1986qd",     0.249,  140.000, true), 
   new eQEDataSetDescription_t(1000922380, "Blatchley:1986qd",     0.274,   60.000, true), 
   new eQEDataSetDescription_t(1000922380, "Blatchley:1986qd",     0.293,  160.000, true), 
   new eQEDataSetDescription_t(1000922380, "Blatchley:1986qd",     0.297,   90.000, true), 
   new eQEDataSetDescription_t(1000922380, "Blatchley:1986qd",     0.309,  134.500, true), 
   new eQEDataSetDescription_t(1000922380, "Blatchley:1986qd",     0.325,  160.000, true), 
   new eQEDataSetDescription_t(1000922380, "Blatchley:1986qd",     0.326,   60.000, true), 
   new eQEDataSetDescription_t(1000922380, "Blatchley:1986qd",     0.335,   90.000, true), 
   new eQEDataSetDescription_t(1000922380, "Blatchley:1986qd",     0.360,  140.000, true), 
   new eQEDataSetDescription_t(1000922380, "Blatchley:1986qd",     0.370,   60.000, true), 
   new eQEDataSetDescription_t(1000922380, "Blatchley:1986qd",     0.423,   60.000, true), 
   new eQEDataSetDescription_t(1000922380, "Blatchley:1986qd",     0.440,   90.000, true), 
   new eQEDataSetDescription_t(1000922380, "Blatchley:1986qd",     0.472,   60.000, true), 
   new eQEDataSetDescription_t(1000922380, "Blatchley:1986qd",     0.491,   90.000, true), 
   new eQEDataSetDescription_t(1000922380, "Blatchley:1986qd",     0.506,   60.000, true), 
   new eQEDataSetDescription_t(1000922380, "Blatchley:1986qd",     0.566,   60.000, true), 
   new eQEDataSetDescription_t(1000922380, "Blatchley:1986qd",     0.690,   60.000, true) 
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
VldTestInputs  gOptGenieInputs;
string         gOptDataFilename = "";

// dbase information
const char * kDefDataFile = "data/validation/eA/xsec/differential/qe/eQE.root";  

// globals
TFile *        gQEDataFile  = 0;
TTree *        gQEDataTree  = 0;
TPostScript *  gPS          = 0;
TCanvas *      gC           = 0;
bool           gShowModel   = false;

// consts

const int kNCx = 2; // number of columns in TCanvas::Divide()
const int kNCy = 2; // number of rows    in TCanvas::Divide()

// model line styles
const int kNMaxNumModels = 5;
const int kLStyle[kNMaxNumModels] = 
{
   1, 2,  3,  5, 6
};
string kLStyleTxt[kNMaxNumModels] = 
{
  "solid", "dashed", "dotted", "dot-dashed", "dot-dot-dashed"
};

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
  gQEDataFile = new TFile(gOptDataFilename.c_str(),"read");  
  gQEDataTree = (TTree *) gQEDataFile->Get("qent");
  if(!gQEDataTree) {
      LOG("gvldtest", pFATAL) 
         << "Can not find TTree `qent' in file: " << gOptDataFilename;
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
  gPS = new TPostScript("eQE.genie_vs_data.ps", 111);

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
  hdr.AddText("Quasielastic Electron-Nucleus Scattering: GENIE vs data");
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

  gQEDataFile->Close();
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

  return 0;
/*

  TFile * xsec_file = gOptGenieInputs.XSecFile(imodel);
  if(!xsec_file) {
     LOG("vldtest", pWARN)
        << "No corresponding cross section file";
     return 0;
  }

  TChain * event_chain = gOptGenieInputs.EvtChain(imodel);
  if(!event_chain) {
     LOG("vldtest", pWARN)
        << "No corresponding event chain.";
     return 0;
  }
  
  // electron energy and scattering angle

  double E     = kDataSet[iset]->E();
  double theta = kDataSet[iset]->Theta();

  int Z = kDataSet[iset]->TgtZ();
  int A = kDataSet[iset]->TgtA();

  double dE       = 0.01; // GeV
  double costheta = TMath::Cos(kPi*theta/180.);

  // total cross section
  TGraph * xsec_em_gr = 0; // get it from the input xsec_file
  if(!xsec_em_gr) {
     LOG("vldtest", pWARN)
        << "Null E/M cross section graph";
     return 0;
  }
  double xsec_em = xsec_em_gr->Eval(E);
  if(xsec_em <= 0) {
     LOG("vldtest", pWARN)
        << "Null E/M cross section graph";
     return 0;
  }
  
  //
  // book histograms
  //

  // E', final state primary lepton energy
  double Ep_min   = 0;
  double Ep_max   = E;
  double Ep_bin   = 0.01;
  double Ep_nbins = (Ep_max-Ep_min)/Ep_bin;

  // theta, scattering angle
  double costh_min   = -1;
  double costh_max   =  1;
  double costh_bin   =  0.01;
  double costh_nbins = (costh_max-costh_min)/costh_bin;
   
  TH2D * h2EpOmega =
     new TH2D("h2EpOmega",  "N=N(E',Omega)|{fixed:E}",
         Ep_nbins, Ep_min, Ep_max, costh_nbins, costh_min, costh_max);
     
  //
  // estimate d^2 sigma / dE' dOmega at the current incoming lepton energy E0
  //

  char cut[500];
  Form(cut,"em&&(fabs(Ev-%d)<%d)", E, dE);

  event_chain->Draw("((pxv*pxl+pyv*pyl+pzv*pzl)/(Ev*El)):El>>h2EpOmega", cut, "GOFF");

  double integral = h2EpOmega->Integral("width");
  if(integral <= 0) {
     LOG("vldtest", pWARN)
        << "Non-positive d^2N / dEp dOmega integral";
     return 0;
  }
  double normalization = 2*kPi*xsec_em/integral;
     
  h2EpOmega->Scale(normalization); // units: 1E-38 cm^2 / GeV /sterad
  
  //
  // now pick a slice at selected theta and return
  // d^2 sigma / dE' dOmega (fixed: E, theta) = f(v = E-E')
  // in the same units as the expt data (nbar/GeV/sterad)
  //

  int N = Ep_nbins;
  
  double * x = new double[N]; // v
  double * y = new double[N]; // d^2 sigma / dE' dOmega
  
  int costheta_bin = h2EpOmega->GetYaxis()->FindBin(costheta);
  
  for(int i = 0; i < h2EpOmega->GetNbinsX(); i++) {
    int Ep_bin = i+1;
    double Ep  = h2EpOmega->GetXaxis()->GetBinCenter(Ep_bin);

    double v      = E - Ep;
    double d2xsec = h2EpOmega->GetBinContent(Ep_bin, costheta_bin);
  
    x[i] = v;
    y[i] = d2xsec;
  }
     
  TGraph * gr = new TGraph(N,x,y);

  delete [] x;
  delete [] y;
    
  return gr;
*/
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

  gQEDataTree->Draw("v:xsec:xsec_err", selection, "goff");
  
  int n = gQEDataTree->GetSelectedRows();

  LOG("gvldtest", pNOTICE) 
    << "Found " << n << " data points in the xsec archive";

  if(n == 0) return 0; // return null graph

  // Data returned by TTree::Draw() are not necessarily ordered in W
  // Do the ordering here before building the graph
  int    *  idx = new int   [n];
  double *  xv  = new double[n];
  double *  yv  = new double[n];
  double *  dyv = new double[n];

  TMath::Sort(n,gQEDataTree->GetV1(),idx,false);

  for(int i=0; i<n; i++) {
     int ii = idx[i];
     xv [i] = (gQEDataTree->GetV1())[ii];
     yv [i] = (gQEDataTree->GetV2())[ii];
     dyv[i] = (gQEDataTree->GetV3())[ii];
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
  hframe->GetXaxis()->SetTitle("#nu = E-E^{'} (GeV)");
  hframe->GetYaxis()->SetTitle("d^{2}#sigma / d#Omega dE (nb/sr/GeV)");
  hframe->GetXaxis()->SetLabelFont(62);
  hframe->GetYaxis()->SetLabelFont(62);
  hframe->GetXaxis()->SetLabelSize(0.04);
  hframe->GetYaxis()->SetLabelSize(0.04);
  hframe->GetXaxis()->SetTitleSize(0.04);
  hframe->GetYaxis()->SetTitleSize(0.04);
  hframe->GetYaxis()->SetTitleOffset(1.65);

/*
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
*/

/*
  // some data show the elastic peak - mark the are to avoid confusion
  if(xmin < 1) {
    double Wm2 = 1.21; // between the QE and Delta peaks
    TBox * qe_peak = new TBox(
       scale_xmin*xmin, scale_ymin*ymin, Wm2, scale_ymax*ymax);
     qe_peak->SetFillColor(kBlue);
     qe_peak->SetFillStyle(3005);
     qe_peak->Draw();
  }
*/

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

  // get data archive
  if(parser.OptionExists('d')){
     string filename = parser.ArgAsString('d');
     gOptDataFilename = filename;
  } else {
     if(gSystem->Getenv("GENIE")) {
        string base_dir = string( gSystem->Getenv("GENIE") );
        string filename = base_dir + "/" + kDefDataFile;
        gOptDataFilename = filename;
     } else { 
        LOG("gvldtest", pFATAL) 
          << "\n Please make sure that $GENIE is defined, or use the -d option"
          << "\n You didn't specify a data file and I can not pick the default one either";
        gAbortingInErr = true;
        exit(1);
     }
  }

  // get GENIE inputs
  gShowModel = true;
  if(parser.OptionExists('g')) {
     string inputs = parser.ArgAsString('g');
     bool ok = gOptGenieInputs.LoadFromFile(inputs);
     if(!ok) {
        LOG("gvldtest", pFATAL)
           << "Could not read validation program inputs: " << inputs;
        gAbortingInErr = true;
        exit(1);
     }
  } else {
    gShowModel = false;
  }
}
//_________________________________________________________________________________
void PrintSyntax(void)
{
  LOG("gvldtest", pNOTICE)
    << "\n\n" << "Syntax:" << "\n"
    << "  gvld_e_qel_xsec -g inputs [-d data_archive_location]\n";
}
//_________________________________________________________________________________
