//____________________________________________________________________________
/*!

\program gvld_e_res_xsec

\brief   Compares GENIE with electron scattering data in the resonance region.

         The data (currently 9531 data points) from E133, E140, E140x, E49A10
         E49A6, E49B, E61, E87, E891, E8920, JLAB, NE11 and ONEN1HAF stored at 
         GENIE's MySQL dbase originate from S.Wood's online dbase.

         Syntax:
           gvld_e_res_xsec 
                [-h host] [-u user] [-p passwd] [-m model]

         Options:

           [] Denotes an optional argument.

           -h NuVld MySQL URL (eg mysql://localhost/NuScat).
           -u NuVld MySQL username.
           -p NuVld MySQL password.
           -m GENIE model

         Example:

            % gvld_e_res_xsec -u costas \
                              -p &^@7287 \
                              -h mysql://localhost/NuScat
                              -m genie::ReinSeghalRESPXSec/Default

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
#include "Conventions/GBuild.h"
#include "Conventions/Constants.h"
#include "Conventions/Units.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGUtils.h"
#include "PDG/PDGCodes.h"
#include "Utils/CmdLnArgParser.h"
#include "Utils/StringUtils.h"
#include "Utils/Style.h"
#include "ValidationTools/NuVld/DBI.h"
#include "ValidationTools/NuVld/DBStatus.h"

using std::ostringstream;
using std::string;

using namespace genie;
using namespace genie::nuvld;
using namespace genie::constants;

/* 
..............................................................................
ELECTRON SCATTERING CROSS SECTION DATA IN THE RESONANCE REGION
..............................................................................
*/
const int    kElXSecDataSets = 152;
const char * kElXSecDataSetLabel[kElXSecDataSets] = 
{
/*   0 */ "JLAB     (Hydrogen),  E =  2.445 GeV, #theta = 20.0^{o}",
/*   1 */ "JLAB     (Hydrogen),  E =  2.445 GeV, #theta = 30.0^{o}",
/*   2 */ "JLAB     (Hydrogen),  E =  2.445 GeV, #theta = 38.5^{o}",
/*   3 */ "JLAB     (Hydrogen),  E =  2.445 GeV, #theta = 70.0^{o}",
/*   4 */ "JLAB     (Hydrogen),  E =  3.245 GeV, #theta = 27.0^{o}",
/*   5 */ "JLAB     (Hydrogen),  E =  4.045 GeV, #theta = 48.0^{o}",
/*   6 */ "JLAB     (Hydrogen),  E =  4.054 GeV, #theta = 24.0^{o}",
/*   7 */ "JLAB     (Hydrogen),  E =  4.054 GeV, #theta = 30.0^{o}",
/*   8 */ "JLAB     (Hydrogen),  E =  4.054 GeV, #theta = 40.0^{o}",
/*   9 */ "JLAB     (Deuterium), E =  2.445 GeV, #theta = 20.0^{o}",
/*  10 */ "JLAB     (Deuterium), E =  2.445 GeV, #theta = 30.0^{o}",
/*  11 */ "JLAB     (Deuterium), E =  2.445 GeV, #theta = 70.0^{o}",
/*  12 */ "JLAB     (Deuterium), E =  3.245 GeV, #theta = 27.0^{o}",
/*  13 */ "JLAB     (Deuterium), E =  4.045 GeV, #theta = 30.0^{o}",
/*  14 */ "JLAB     (Deuterium), E =  4.045 GeV, #theta = 40.0^{o}",
/*  15 */ "JLAB     (Deuterium), E =  4.045 GeV, #theta = 48.0^{o}",
/*  16 */ "JLAB     (Deuterium), E =  4.054 GeV, #theta = 24.0^{o}",
/*  17 */ "NE11     (Hydrogen),  E =  5.507 GeV, #theta = 15.1^{o}",
/*  18 */ "NE11     (Hydrogen),  E =  5.507 GeV, #theta = 19.0^{o}",
/*  19 */ "NE11     (Hydrogen),  E =  5.507 GeV, #theta = 22.8^{o}",
/*  20 */ "NE11     (Hydrogen),  E =  5.507 GeV, #theta = 26.8^{o}",
/*  21 */ "NE11     (Hydrogen),  E =  9.800 GeV, #theta = 13.2^{o}",
/*  22 */ "NE11     (Hydrogen),  E =  9.800 GeV, #theta = 15.4^{o}",
/*  23 */ "NE11     (Hydrogen),  E =  9.800 GeV, #theta = 17.5^{o}",
/*  24 */ "NE11     (Hydrogen),  E =  9.800 GeV, #theta = 19.8^{o}",
/*  25 */ "NE11     (Deuterium), E =  5.507 GeV, #theta = 15.1^{o}",
/*  26 */ "NE11     (Deuterium), E =  5.507 GeV, #theta = 19.0^{o}",
/*  27 */ "NE11     (Deuterium), E =  5.507 GeV, #theta = 22.8^{o}",
/*  28 */ "NE11     (Deuterium), E =  5.507 GeV, #theta = 26.8^{o}",
/*  29 */ "E133     (Hydrogen),  E =  9.772 GeV, #theta = 10.0^{o}",
/*  30 */ "E133     (Hydrogen),  E = 12.604 GeV, #theta = 10.0^{o}",
/*  31 */ "E133     (Hydrogen),  E = 15.740 GeV, #theta = 10.0^{o}",
/*  32 */ "E133     (Hydrogen),  E = 18.520 GeV, #theta = 10.0^{o}",
/*  33 */ "E133     (Hydrogen),  E = 20.998 GeV, #theta = 10.0^{o}",
/*  34 */ "E133     (Deuterium), E =  9.766 GeV, #theta = 10.0^{o}",
/*  35 */ "E133     (Deuterium), E = 12.603 GeV, #theta = 10.0^{o}",
/*  36 */ "E133     (Deuterium), E = 15.725 GeV, #theta = 10.0^{o}",
/*  37 */ "E133     (Deuterium), E = 17.301 GeV, #theta = 10.0^{o}",
/*  38 */ "E133     (Deuterium), E = 18.498 GeV, #theta = 10.0^{o}",
/*  39 */ "E133     (Deuterium), E = 20.998 GeV, #theta = 10.0^{o}",
/*  40 */ "E140x    (Hydrogen),  E =  1.997 GeV, #theta = 46.2^{o}",
/*  41 */ "E140x    (Hydrogen),  E =  3.301 GeV, #theta = 50.0^{o}",
/*  42 */ "E140x    (Hydrogen),  E =  4.014 GeV, #theta = 50.7^{o}",
/*  43 */ "E140x    (Hydrogen),  E =  4.951 GeV, #theta = 56.0^{o}",
/*  44 */ "E140x    (Hydrogen),  E =  5.600 GeV, #theta = 58.0^{o}",
/*  45 */ "E140x    (Hydrogen),  E =  6.453 GeV, #theta = 59.5^{o}",
/*  46 */ "E140x    (Deuterium), E =  3.301 GeV, #theta = 50.0^{o}",
/*  47 */ "E140x    (Deuterium), E =  4.014 GeV, #theta = 50.7^{o}",
/*  48 */ "E140x    (Deuterium), E =  4.951 GeV, #theta = 56.0^{o}",
/*  49 */ "E140x    (Deuterium), E =  5.600 GeV, #theta = 58.0^{o}",
/*  50 */ "E140x    (Deuterium), E =  6.453 GeV, #theta = 59.5^{o}",
/*  51 */ "E49A6    (Hydrogen),  E =  4.511 GeV, #theta =  6.0^{o}",
/*  52 */ "E49A6    (Hydrogen),  E =  7.014 GeV, #theta =  6.0^{o}",
/*  53 */ "E49A6    (Hydrogen),  E = 10.027 GeV, #theta =  6.0^{o}",
/*  54 */ "E49A6    (Hydrogen),  E = 13.549 GeV, #theta =  6.0^{o}",
/*  55 */ "E49A6    (Hydrogen),  E = 16.075 GeV, #theta =  6.0^{o}",
/*  56 */ "E49A6    (Hydrogen),  E = 19.544 GeV, #theta =  6.0^{o}",
/*  57 */ "E49A6    (Deuterium), E =  4.511 GeV, #theta =  6.0^{o}",
/*  58 */ "E49A6    (Deuterium), E =  7.014 GeV, #theta =  6.0^{o}",
/*  59 */ "E49A6    (Deuterium), E = 10.027 GeV, #theta =  6.0^{o}",
/*  60 */ "E49A6    (Deuterium), E = 13.549 GeV, #theta =  6.0^{o}",
/*  61 */ "E49A6    (Deuterium), E = 16.075 GeV, #theta =  6.0^{o}",
/*  62 */ "E49A6    (Deuterium), E = 19.544 GeV, #theta =  6.0^{o}",
/*  63 */ "E49A10   (Hydrogen),  E =  4.892 GeV, #theta = 10.0^{o}",
/*  64 */ "E49A10   (Hydrogen),  E =  7.019 GeV, #theta = 10.0^{o}",
/*  65 */ "E49A10   (Hydrogen),  E =  9.022 GeV, #theta = 10.0^{o}",
/*  66 */ "E49A10   (Hydrogen),  E = 10.998 GeV, #theta = 10.0^{o}",
/*  67 */ "E49A10   (Hydrogen),  E = 13.545 GeV, #theta = 10.0^{o}",
/*  68 */ "E49A10   (Hydrogen),  E = 15.204 GeV, #theta = 10.0^{o}",
/*  69 */ "E49A10   (Hydrogen),  E = 17.706 GeV, #theta = 10.0^{o}",
/*  70 */ "E49A10   (Hydrogen),  E = 19.350 GeV, #theta = 10.0^{o}",
/*  71 */ "E49A10   (Deuterium), E =  4.892 GeV, #theta = 10.0^{o}",
/*  72 */ "E49A10   (Deuterium), E =  7.019 GeV, #theta = 10.0^{o}",
/*  73 */ "E49A10   (Deuterium), E =  9.022 GeV, #theta = 10.0^{o}",
/*  74 */ "E49A10   (Deuterium), E = 10.998 GeV, #theta = 10.0^{o}",
/*  75 */ "E49A10   (Deuterium), E = 13.545 GeV, #theta = 10.0^{o}",
/*  76 */ "E49A10   (Deuterium), E = 15.204 GeV, #theta = 10.0^{o}",
/*  77 */ "E49A10   (Deuterium), E = 17.706 GeV, #theta = 10.0^{o}",
/*  78 */ "E49A10   (Deuterium), E = 19.350 GeV, #theta = 10.0^{o}",
/*  79 */ "E49B     (Hydrogen),  E =  4.502 GeV, #theta = 34.0^{o}",
/*  80 */ "E49B     (Hydrogen),  E =  4.504 GeV, #theta = 18.0^{o}",
/*  81 */ "E49B     (Hydrogen),  E =  4.506 GeV, #theta = 26.0^{o}",
/*  82 */ "E49B     (Hydrogen),  E =  5.808 GeV, #theta = 34.0^{o}",
/*  83 */ "E49B     (Hydrogen),  E =  6.509 GeV, #theta = 18.0^{o}",
/*  84 */ "E49B     (Hydrogen),  E =  6.711 GeV, #theta = 26.0^{o}",
/*  85 */ "E49B     (Hydrogen),  E =  7.912 GeV, #theta = 34.0^{o}",
/*  86 */ "E49B     (Hydrogen),  E =  8.614 GeV, #theta = 18.0^{o}",
/*  87 */ "E49B     (Deuterium), E =  4.502 GeV, #theta = 34.0^{o}",
/*  88 */ "E49B     (Deuterium), E =  4.504 GeV, #theta = 18.0^{o}",
/*  89 */ "E49B     (Deuterium), E =  4.506 GeV, #theta = 26.0^{o}",
/*  90 */ "E49B     (Deuterium), E =  5.808 GeV, #theta = 34.0^{o}",
/*  91 */ "E49B     (Deuterium), E =  6.509 GeV, #theta = 18.0^{o}",
/*  92 */ "E49B     (Deuterium), E =  6.711 GeV, #theta = 26.0^{o}",
/*  93 */ "E49B     (Deuterium), E =  7.912 GeV, #theta = 34.0^{o}",
/*  94 */ "E49B     (Deuterium), E =  8.614 GeV, #theta = 18.0^{o}",
/*  95 */ "E61      (Hydrogen),  E =  4.499 GeV, #theta =  4.0^{o}",
/*  96 */ "E61      (Hydrogen),  E =  7.000 GeV, #theta =  4.0^{o}",
/*  97 */ "E61      (Hydrogen),  E =  9.993 GeV, #theta =  4.0^{o}",
/*  98 */ "E61      (Hydrogen),  E = 13.000 GeV, #theta =  4.0^{o}",
/*  99 */ "E61      (Hydrogen),  E = 16.000 GeV, #theta =  4.0^{o}",
/* 100 */ "E61      (Hydrogen),  E = 18.010 GeV, #theta =  4.0^{o}",
/* 101 */ "E61      (Hydrogen),  E = 20.005 GeV, #theta =  4.0^{o}",
/* 102 */ "E61      (Deuterium), E =  4.499 GeV, #theta =  4.0^{o}",
/* 103 */ "E61      (Deuterium), E =  6.999 GeV, #theta =  4.0^{o}",
/* 104 */ "E61      (Deuterium), E =  9.993 GeV, #theta =  4.0^{o}",
/* 105 */ "E61      (Deuterium), E = 12.987 GeV, #theta =  4.0^{o}",
/* 106 */ "E61      (Deuterium), E = 16.000 GeV, #theta =  4.0^{o}",
/* 107 */ "E61      (Deuterium), E = 18.010 GeV, #theta =  4.0^{o}",
/* 108 */ "E61      (Deuterium), E = 20.005 GeV, #theta =  4.0^{o}",
/* 109 */ "E891     (Hydrogen),  E =  6.500 GeV, #theta = 60.0^{o}",
/* 110 */ "E891     (Hydrogen),  E =  7.000 GeV, #theta = 50.0^{o}",
/* 111 */ "E891     (Hydrogen),  E = 10.400 GeV, #theta = 60.0^{o}",
/* 112 */ "E891     (Hydrogen),  E = 13.290 GeV, #theta = 60.0^{o}",
/* 113 */ "E891     (Hydrogen),  E = 16.002 GeV, #theta = 60.0^{o}",
/* 114 */ "E891     (Hydrogen),  E = 19.505 GeV, #theta = 60.0^{o}",
/* 115 */ "E891     (Deuterium), E =  6.500 GeV, #theta = 60.0^{o}",
/* 116 */ "E891     (Deuterium), E =  7.000 GeV, #theta = 50.0^{o}",
/* 117 */ "E891     (Deuterium), E = 10.400 GeV, #theta = 60.0^{o}",
/* 118 */ "E891     (Deuterium), E = 13.290 GeV, #theta = 60.0^{o}",
/* 119 */ "E891     (Deuterium), E = 16.002 GeV, #theta = 60.0^{o}",
/* 120 */ "E891     (Deuterium), E = 19.505 GeV, #theta = 60.0^{o}",
/* 121 */ "E8920    (Hydrogen),  E =  6.500 GeV, #theta = 18.0^{o}",
/* 122 */ "E8920    (Hydrogen),  E =  7.000 GeV, #theta =  6.0^{o}",
/* 123 */ "E8920    (Hydrogen),  E = 10.400 GeV, #theta = 18.0^{o}",
/* 124 */ "E8920    (Hydrogen),  E = 13.300 GeV, #theta = 18.0^{o}",
/* 125 */ "E8920    (Hydrogen),  E = 13.500 GeV, #theta =  6.0^{o}",
/* 126 */ "E8920    (Hydrogen),  E = 16.000 GeV, #theta = 18.0^{o}",
/* 127 */ "E8920    (Hydrogen),  E = 16.000 GeV, #theta = 15.0^{o}",
/* 128 */ "E8920    (Hydrogen),  E = 16.000 GeV, #theta =  6.0^{o}",
/* 129 */ "E8920    (Hydrogen),  E = 19.500 GeV, #theta = 20.6^{o}",
/* 130 */ "E8920    (Hydrogen),  E = 19.500 GeV, #theta = 18.0^{o}",
/* 131 */ "E8920    (Hydrogen),  E = 19.500 GeV, #theta =  6.0^{o}",
/* 132 */ "E8920    (Deuterium), E =  6.500 GeV, #theta = 18.0^{o}",
/* 133 */ "E8920    (Deuterium), E =  7.000 GeV, #theta =  6.0^{o}",
/* 134 */ "E8920    (Deuterium), E = 10.400 GeV, #theta = 18.0^{o}",
/* 135 */ "E8920    (Deuterium), E = 13.300 GeV, #theta = 18.0^{o}",
/* 136 */ "E8920    (Deuterium), E = 13.500 GeV, #theta =  6.0^{o}",
/* 137 */ "E8920    (Deuterium), E = 16.000 GeV, #theta = 18.0^{o}",
/* 138 */ "E8920    (Deuterium), E = 16.000 GeV, #theta = 15.0^{o}",
/* 139 */ "E8920    (Deuterium), E = 16.000 GeV, #theta =  6.0^{o}",
/* 140 */ "E8920    (Deuterium), E = 19.500 GeV, #theta = 20.6^{o}",
/* 141 */ "E8920    (Deuterium), E = 19.500 GeV, #theta = 18.0^{o}",
/* 142 */ "E8920    (Deuterium), E = 19.500 GeV, #theta =  6.0^{o}",
/* 143 */ "ONEN1HAF (Hydrogen),  E =  5.000 GeV, #theta =  1.5^{o}",
/* 144 */ "ONEN1HAF (Hydrogen),  E =  7.103 GeV, #theta =  1.5^{o}",
/* 145 */ "ONEN1HAF (Hydrogen),  E =  9.301 GeV, #theta =  1.5^{o}",
/* 146 */ "ONEN1HAF (Hydrogen),  E = 11.799 GeV, #theta =  1.5^{o}",
/* 147 */ "ONEN1HAF (Hydrogen),  E = 13.805 GeV, #theta =  1.5^{o}",
/* 148 */ "ONEN1HAF (Hydrogen),  E = 15.604 GeV, #theta =  1.5^{o}",
/* 149 */ "ONEN1HAF (Hydrogen),  E = 17.298 GeV, #theta =  1.5^{o}",
/* 150 */ "ONEN1HAF (Hydrogen),  E = 18.703 GeV, #theta =  1.5^{o}",
/* 151 */ "ONEN1HAF (Hydrogen),  E = 19.995 GeV, #theta =  1.5^{o}"
};
const char * kElXSecKeyList[kElXSecDataSets] = {
/*   0 */ "JLAB,0",
/*   1 */ "JLAB,0",
/*   2 */ "JLAB,0",
/*   3 */ "JLAB,0",
/*   4 */ "JLAB,0",
/*   5 */ "JLAB,0",
/*   6 */ "JLAB,0",
/*   7 */ "JLAB,0",
/*   8 */ "JLAB,0",
/*   9 */ "JLAB,1",
/*  10 */ "JLAB,1",
/*  11 */ "JLAB,1",
/*  12 */ "JLAB,1",
/*  13 */ "JLAB,1",
/*  14 */ "JLAB,1",
/*  15 */ "JLAB,1",
/*  16 */ "JLAB,1",
/*  17 */ "NE11,0",
/*  18 */ "NE11,0",
/*  19 */ "NE11,0",
/*  20 */ "NE11,0",
/*  21 */ "NE11,0",
/*  22 */ "NE11,0",
/*  23 */ "NE11,0",
/*  24 */ "NE11,0",
/*  25 */ "NE11,1",
/*  26 */ "NE11,1",
/*  27 */ "NE11,1",
/*  28 */ "NE11,1",
/*  29 */ "E133,0",
/*  30 */ "E133,0",
/*  31 */ "E133,0",
/*  32 */ "E133,0",
/*  33 */ "E133,0",
/*  34 */ "E133,1",
/*  35 */ "E133,1",
/*  36 */ "E133,1",
/*  37 */ "E133,1",
/*  38 */ "E133,1",
/*  39 */ "E133,1",
/*  40 */ "E140x,0",
/*  41 */ "E140x,0",
/*  42 */ "E140x,0",
/*  43 */ "E140x,0",
/*  44 */ "E140x,0",
/*  45 */ "E140x,0",
/*  46 */ "E140x,1",
/*  47 */ "E140x,1",
/*  48 */ "E140x,1",
/*  49 */ "E140x,1",
/*  50 */ "E140x,1",
/*  51 */ "E49A6,0",
/*  52 */ "E49A6,0",
/*  53 */ "E49A6,0",
/*  54 */ "E49A6,0",
/*  55 */ "E49A6,0",
/*  56 */ "E49A6,0",
/*  57 */ "E49A6,1",
/*  58 */ "E49A6,1",
/*  59 */ "E49A6,1",
/*  60 */ "E49A6,1",
/*  61 */ "E49A6,1",
/*  62 */ "E49A6,1",
/*  63 */ "E49A10,0",
/*  64 */ "E49A10,0",
/*  65 */ "E49A10,0",
/*  66 */ "E49A10,0",
/*  67 */ "E49A10,0",
/*  68 */ "E49A10,0",
/*  69 */ "E49A10,0",
/*  70 */ "E49A10,0",
/*  71 */ "E49A10,1",
/*  72 */ "E49A10,1",
/*  73 */ "E49A10,1",
/*  74 */ "E49A10,1",
/*  75 */ "E49A10,1",
/*  76 */ "E49A10,1",
/*  77 */ "E49A10,1",
/*  78 */ "E49A10,1",
/*  79 */ "E49B,0",
/*  80 */ "E49B,0",
/*  81 */ "E49B,0",
/*  82 */ "E49B,0",
/*  83 */ "E49B,0",
/*  84 */ "E49B,0",
/*  85 */ "E49B,0",
/*  86 */ "E49B,0",
/*  87 */ "E49B,1",
/*  88 */ "E49B,1",
/*  89 */ "E49B,1",
/*  90 */ "E49B,1",
/*  91 */ "E49B,1",
/*  92 */ "E49B,1",
/*  93 */ "E49B,1",
/*  94 */ "E49B,1",
/*  95 */ "E61,0",
/*  96 */ "E61,0",
/*  97 */ "E61,0",
/*  98 */ "E61,0",
/*  99 */ "E61,0",
/* 100 */ "E61,0",
/* 101 */ "E61,0",
/* 102 */ "E61,1",
/* 103 */ "E61,1",
/* 104 */ "E61,1",
/* 105 */ "E61,1",
/* 106 */ "E61,1",
/* 107 */ "E61,1",
/* 108 */ "E61,1",
/* 109 */ "E891,0",
/* 110 */ "E891,0",
/* 111 */ "E891,0",
/* 112 */ "E891,0",
/* 113 */ "E891,0",
/* 114 */ "E891,0",
/* 115 */ "E891,1",
/* 116 */ "E891,1",
/* 117 */ "E891,1",
/* 118 */ "E891,1",
/* 119 */ "E891,1",
/* 120 */ "E891,1",
/* 121 */ "E8920,0",
/* 122 */ "E8920,0",
/* 123 */ "E8920,0",
/* 124 */ "E8920,0",
/* 125 */ "E8920,0",
/* 126 */ "E8920,0",
/* 127 */ "E8920,0",
/* 128 */ "E8920,0",
/* 129 */ "E8920,0",
/* 130 */ "E8920,0",
/* 131 */ "E8920,0",
/* 132 */ "E8920,1",
/* 133 */ "E8920,1",
/* 134 */ "E8920,1",
/* 135 */ "E8920,1",
/* 136 */ "E8920,1",
/* 137 */ "E8920,1",
/* 138 */ "E8920,1",
/* 139 */ "E8920,1",
/* 140 */ "E8920,1",
/* 141 */ "E8920,1",
/* 142 */ "E8920,1",
/* 143 */ "ONEN1HAF,0",
/* 144 */ "ONEN1HAF,0",
/* 145 */ "ONEN1HAF,0",
/* 146 */ "ONEN1HAF,0",
/* 147 */ "ONEN1HAF,0",
/* 148 */ "ONEN1HAF,0",
/* 149 */ "ONEN1HAF,0",
/* 150 */ "ONEN1HAF,0",
/* 151 */ "ONEN1HAF,0"
};
float kElXSecEnergy[kElXSecDataSets] = {
/*   0 */   2.445,
/*   1 */   2.445,
/*   2 */   2.445,
/*   3 */   2.445,
/*   4 */   3.245,
/*   5 */   4.045,
/*   6 */   4.054,
/*   7 */   4.054,
/*   8 */   4.054,
/*   9 */   2.445,
/*  10 */   2.445,
/*  11 */   2.445,
/*  12 */   3.245,
/*  13 */   4.045,
/*  14 */   4.045,
/*  15 */   4.045,
/*  16 */   4.054,
/*  17 */   5.507,
/*  18 */   5.507,
/*  19 */   5.507,
/*  20 */   5.507,
/*  21 */   9.800,
/*  22 */   9.800,
/*  23 */   9.800,
/*  24 */   9.800,
/*  25 */   5.507,
/*  26 */   5.507,
/*  27 */   5.507,
/*  28 */   5.507,
/*  29 */   9.772,
/*  30 */  12.604,
/*  31 */  15.740,
/*  32 */  18.520,
/*  33 */  20.998,
/*  34 */   9.766,
/*  35 */  12.603,
/*  36 */  15.725,
/*  37 */  17.301,
/*  38 */  18.498,
/*  39 */  20.998,
/*  40 */   1.997,
/*  41 */   3.301,
/*  42 */   4.014,
/*  43 */   4.951,
/*  44 */   5.600,
/*  45 */   6.453,
/*  46 */   3.301,
/*  47 */   4.014,
/*  48 */   4.951,
/*  49 */   5.600,
/*  50 */   6.453,
/*  51 */   4.511,
/*  52 */   7.014,
/*  53 */  10.027,
/*  54 */  13.549,
/*  55 */  16.075,
/*  56 */  19.544,
/*  57 */   4.511,
/*  58 */   7.014,
/*  59 */  10.027,
/*  60 */  13.549,
/*  61 */  16.075,
/*  62 */  19.544,
/*  63 */   4.892,
/*  64 */   7.019,
/*  65 */   9.022,
/*  66 */  10.998,
/*  67 */  13.545,
/*  68 */  15.204,
/*  69 */  17.706,
/*  70 */  19.350,
/*  71 */   4.892,
/*  72 */   7.019,
/*  73 */   9.022,
/*  74 */  10.998,
/*  75 */  13.545,
/*  76 */  15.204,
/*  77 */  17.706,
/*  78 */  19.350,
/*  79 */   4.502,
/*  80 */   4.504,
/*  81 */   4.506,
/*  82 */   5.808,
/*  83 */   6.509,
/*  84 */   6.711,
/*  85 */   7.912,
/*  86 */   8.614,
/*  87 */   4.502,
/*  88 */   4.504,
/*  89 */   4.506,
/*  90 */   5.808,
/*  91 */   6.509,
/*  92 */   6.711,
/*  93 */   7.912,
/*  94 */   8.614,
/*  95 */   4.499,
/*  96 */   7.000,
/*  97 */   9.993,
/*  98 */  13.000,
/*  99 */  16.000,
/* 100 */  18.010,
/* 101 */  20.005,
/* 102 */   4.499,
/* 103 */   6.999,
/* 104 */   9.993,
/* 105 */  12.987,
/* 106 */  16.000,
/* 107 */  18.010,
/* 108 */  20.005,
/* 109 */   6.500,
/* 110 */   7.000,
/* 111 */  10.400,
/* 112 */  13.290,
/* 113 */  16.002,
/* 114 */  19.505,
/* 115 */   6.500,
/* 116 */   7.000,
/* 117 */  10.400,
/* 118 */  13.290,
/* 119 */  16.002,
/* 120 */  19.505,
/* 121 */   6.500,
/* 122 */   7.000,
/* 123 */  10.400,
/* 124 */  13.300,
/* 125 */  13.500,
/* 126 */  16.000,
/* 127 */  16.000,
/* 128 */  16.000,
/* 129 */  19.500,
/* 130 */  19.500,
/* 131 */  19.500,
/* 132 */   6.500,
/* 133 */   7.000,
/* 134 */  10.400,
/* 135 */  13.300,
/* 136 */  13.500,
/* 137 */  16.000,
/* 138 */  16.000,
/* 139 */  16.000,
/* 140 */  19.500,
/* 141 */  19.500,
/* 142 */  19.500,
/* 143 */   5.000,
/* 144 */   7.103,
/* 145 */   9.301,
/* 146 */  11.799,
/* 147 */  13.805,
/* 148 */  15.604,
/* 149 */  17.298,
/* 150 */  18.703,
/* 151 */  19.995
};
float kElXSecTheta[kElXSecDataSets] = {
/*   0 */  20.000,
/*   1 */  30.000,
/*   2 */  38.500,
/*   3 */  70.010,
/*   4 */  26.980,
/*   5 */  47.990,
/*   6 */  24.030,
/*   7 */  30.000,
/*   8 */  39.990,
/*   9 */  20.000,
/*  10 */  30.000,
/*  11 */  70.010,
/*  12 */  26.980,
/*  13 */  30.000,
/*  14 */  39.990,
/*  15 */  48.000,
/*  16 */  24.030,
/*  17 */  15.146,
/*  18 */  18.981,
/*  19 */  22.805,
/*  20 */  26.823,
/*  21 */  13.248,
/*  22 */  15.367,
/*  23 */  17.516,
/*  24 */  19.753,
/*  25 */  15.146,
/*  26 */  18.981,
/*  27 */  22.805,
/*  28 */  26.823,
/*  29 */  10.000, 
/*  30 */  10.000,
/*  31 */  10.000,
/*  32 */  10.000,
/*  33 */  10.000,
/*  34 */  10.000,
/*  35 */  10.000,
/*  36 */  10.000,
/*  37 */  10.000,
/*  38 */  10.000,
/*  39 */  10.000,
/*  40 */  46.160,
/*  41 */  50.000,
/*  42 */  50.700,
/*  43 */  56.010,
/*  44 */  58.000,
/*  45 */  59.510,
/*  46 */  50.000,
/*  47 */  50.700,
/*  48 */  56.010,
/*  49 */  58.000,
/*  50 */  59.510,
/*  51 */   5.988,
/*  52 */   5.988,
/*  53 */   5.988,
/*  54 */   5.988,
/*  55 */   5.988,
/*  56 */   5.988,
/*  57 */   5.988,
/*  58 */   5.988,
/*  59 */   5.988,
/*  60 */   5.988,
/*  61 */   5.988,
/*  62 */   5.988,
/*  63 */  10.000,
/*  64 */  10.000,
/*  65 */  10.000,
/*  66 */  10.000,
/*  67 */  10.000,
/*  68 */  10.000,
/*  69 */  10.000,
/*  70 */  10.000,
/*  71 */  10.000,
/*  72 */  10.000,
/*  73 */  10.000,
/*  74 */  10.000,
/*  75 */  10.000,
/*  76 */  10.000,
/*  77 */  10.000,
/*  78 */  10.000,
/*  79 */  34.009,
/*  80 */  18.020,
/*  81 */  26.015,
/*  82 */  34.009,
/*  83 */  18.020,
/*  84 */  26.015,
/*  85 */  34.009,
/*  86 */  18.020,
/*  87 */  34.009,
/*  88 */  18.020,
/*  89 */  26.015,
/*  90 */  34.009,
/*  91 */  18.020,
/*  92 */  26.015,
/*  93 */  34.009,
/*  94 */  18.020,
/*  95 */   4.000,
/*  96 */   4.000,
/*  97 */   4.000,
/*  98 */   4.000,
/*  99 */   4.000,
/* 100 */   4.000,
/* 101 */   4.000,
/* 102 */   4.000,
/* 103 */   4.000,
/* 104 */   4.000,
/* 105 */   4.000,
/* 106 */   4.000,
/* 107 */   4.000,
/* 108 */   4.000,
/* 109 */  59.999,
/* 110 */  49.998,
/* 111 */  60.000,
/* 112 */  60.000,
/* 113 */  59.999,
/* 114 */  60.000,
/* 115 */  59.999,
/* 116 */  49.998,
/* 117 */  60.000,
/* 118 */  60.000,
/* 119 */  59.999,
/* 120 */  60.000,
/* 121 */  18.000,
/* 122 */   6.000,
/* 123 */  18.000,
/* 124 */  18.000,
/* 125 */   6.000,
/* 126 */  18.000,
/* 127 */  15.000,
/* 128 */   6.000,
/* 129 */  20.600,
/* 130 */  18.000,
/* 131 */   6.000,
/* 132 */  18.000,
/* 133 */   6.000,
/* 134 */  18.000,
/* 135 */  18.000,
/* 136 */   6.000,
/* 137 */  18.000,
/* 138 */  15.000,
/* 139 */   6.000,
/* 140 */  20.600,
/* 141 */  18.000,
/* 142 */   6.000,
/* 143 */   1.482,
/* 144 */   1.482,
/* 145 */   1.482,
/* 146 */   1.482,
/* 147 */   1.482,
/* 148 */   1.482,
/* 149 */   1.482,
/* 150 */   1.482,
/* 151 */   1.482
};
int kElXSecTarget[kElXSecDataSets] = {
/*   0 */ 1000010010,           
/*   1 */ 1000010010,
/*   2 */ 1000010010,
/*   3 */ 1000010010,
/*   4 */ 1000010010,
/*   5 */ 1000010010,
/*   6 */ 1000010010,
/*   7 */ 1000010010,
/*   8 */ 1000010010,
/*   9 */ 1000020010,
/*  10 */ 1000020010,
/*  11 */ 1000020010,
/*  12 */ 1000020010,
/*  13 */ 1000020010,
/*  14 */ 1000020010,
/*  15 */ 1000020010,
/*  16 */ 1000020010,
/*  17 */ 1000010010,
/*  18 */ 1000010010,
/*  19 */ 1000010010,
/*  20 */ 1000010010,
/*  21 */ 1000010010,
/*  22 */ 1000010010,
/*  23 */ 1000010010,
/*  24 */ 1000010010,
/*  25 */ 1000020010,
/*  26 */ 1000020010,
/*  27 */ 1000020010,
/*  28 */ 1000020010,
/*  29 */ 1000010010,
/*  30 */ 1000010010,
/*  31 */ 1000010010,
/*  32 */ 1000010010,
/*  33 */ 1000010010,
/*  34 */ 1000020010,
/*  35 */ 1000020010,
/*  36 */ 1000020010,
/*  37 */ 1000020010,
/*  38 */ 1000020010,
/*  39 */ 1000020010,
/*  40 */ 1000010010, 
/*  41 */ 1000010010, 
/*  42 */ 1000010010, 
/*  43 */ 1000010010, 
/*  44 */ 1000010010, 
/*  45 */ 1000010010, 
/*  46 */ 1000020010,
/*  47 */ 1000020010,
/*  48 */ 1000020010,
/*  49 */ 1000020010,
/*  50 */ 1000020010,
/*  51 */ 1000010010,
/*  52 */ 1000010010,
/*  53 */ 1000010010,
/*  54 */ 1000010010,
/*  55 */ 1000010010,
/*  56 */ 1000010010,
/*  57 */ 1000020010,
/*  58 */ 1000020010,
/*  59 */ 1000020010,
/*  60 */ 1000020010,
/*  61 */ 1000020010,
/*  62 */ 1000020010,
/*  63 */ 1000010010,
/*  64 */ 1000010010,
/*  65 */ 1000010010,
/*  66 */ 1000010010,
/*  67 */ 1000010010,
/*  68 */ 1000010010,
/*  69 */ 1000010010,
/*  70 */ 1000010010,
/*  71 */ 1000020010,
/*  72 */ 1000020010,
/*  73 */ 1000020010,
/*  74 */ 1000020010,
/*  75 */ 1000020010,
/*  76 */ 1000020010,
/*  77 */ 1000020010,
/*  78 */ 1000020010,
/*  79 */ 1000010010,
/*  80 */ 1000010010,
/*  81 */ 1000010010,
/*  82 */ 1000010010,
/*  83 */ 1000010010,
/*  84 */ 1000010010,
/*  85 */ 1000010010,
/*  86 */ 1000010010,
/*  87 */ 1000020010,
/*  88 */ 1000020010,
/*  89 */ 1000020010,
/*  90 */ 1000020010,
/*  91 */ 1000020010,
/*  92 */ 1000020010,
/*  93 */ 1000020010,
/*  94 */ 1000020010,
/*  95 */ 1000010010, 
/*  96 */ 1000010010, 
/*  97 */ 1000010010, 
/*  98 */ 1000010010, 
/*  99 */ 1000010010, 
/* 100 */ 1000010010, 
/* 101 */ 1000010010, 
/* 102 */ 1000020010,
/* 103 */ 1000020010,
/* 104 */ 1000020010,
/* 105 */ 1000020010,
/* 106 */ 1000020010,
/* 107 */ 1000020010,
/* 108 */ 1000020010,
/* 109 */ 1000010010,
/* 110 */ 1000010010,
/* 111 */ 1000010010,
/* 112 */ 1000010010,
/* 113 */ 1000010010,
/* 114 */ 1000010010,
/* 115 */ 1000020010,
/* 116 */ 1000020010,
/* 117 */ 1000020010,
/* 118 */ 1000020010,
/* 119 */ 1000020010,
/* 120 */ 1000020010,
/* 121 */ 1000010010,
/* 122 */ 1000010010,
/* 123 */ 1000010010,
/* 124 */ 1000010010,
/* 125 */ 1000010010,
/* 126 */ 1000010010,
/* 127 */ 1000010010,
/* 128 */ 1000010010,
/* 129 */ 1000010010,
/* 130 */ 1000010010,
/* 131 */ 1000010010,
/* 132 */ 1000020010,
/* 133 */ 1000020010,
/* 134 */ 1000020010,
/* 135 */ 1000020010,
/* 136 */ 1000020010,
/* 137 */ 1000020010,
/* 138 */ 1000020010,
/* 139 */ 1000020010,
/* 140 */ 1000020010,
/* 141 */ 1000020010,
/* 142 */ 1000020010,
/* 143 */ 1000010010,
/* 144 */ 1000010010,
/* 145 */ 1000010010,
/* 146 */ 1000010010,
/* 147 */ 1000010010,
/* 148 */ 1000010010,
/* 149 */ 1000010010,
/* 150 */ 1000010010,
/* 151 */ 1000010010
};

typedef DBQueryString                 DBQ;
typedef DBTable<DBElDiffXSecTableRow> DBT;

// function prototypes
void     Init               (void);
void     Plot               (void);
void     End                (void);
void     AddCoverPage       (void);
bool     Connect            (void);
DBQ      FormQuery          (const char * key_list, float energy, float theta);
DBT *    Data               (int iset);
TGraph * Model              (int iset, int imodel);
void     Draw               (int iset);
void     GetCommandLineArgs (int argc, char ** argv);
void     PrintSyntax        (void);

// command-line arguments
string gOptDbURL;
string gOptDbUser;
string gOptDbPasswd;
string gOptGModelName;
string gOptGModelConf;

// dbase information
const char * kDefDbURL = "mysql://localhost/NuScat";  

// globals
bool            gCmpWithData  = true;
DBI *           gDBI          = 0;
TPostScript *   gPS           = 0;
TCanvas *       gC            = 0;
bool            gShowModel    = false;

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
#ifdef __GENIE_MYSQL_ENABLED__

  // connect to the NuValidator MySQL dbase
  bool ok = Connect();
  if(!ok) {
    return;
  }
 
  // loop over data sets
  for(int iset = 0; iset < kElXSecDataSets; iset++) 
  {
    Draw(iset);
  }
#endif
}
//_________________________________________________________________________________
void Init(void)
{
  LOG("vldtest", pNOTICE) << "Initializing...";

  // genie style
  utils::style::SetDefaultStyle();

  // canvas
  gC = new TCanvas("c","",20,20,500,650);
  gC->SetBorderMode(0);
  gC->SetFillColor(0);
  gC->SetGridx();
  gC->SetGridy();

  // output postscript file
  gPS = new TPostScript("genie_eres_vs_data.ps", 111);

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
  hdr.AddText("GENIE Resonance Electro-Production vs Data");
  hdr.AddText(" ");
  hdr.Draw();
  gC->Update();
}
//_________________________________________________________________________________
void End(void)
{
  LOG("vldtest", pNOTICE) << "Cleaning up...";

  gPS->Close();

  delete gC;
  delete gPS;
}
//_________________________________________________________________________________
// Corresponding GENIE prediction for the `iset' data set 
//.................................................................................
TGraph * Model(int iset, int imodel)
{
  if(!gShowModel) return 0;

  LOG("vldtest", pNOTICE) 
    << "Getting GENIE prediction (model ID = " 
    << imodel << ", data set ID = " << iset << ")";

  AlgFactory * algf = AlgFactory::Instance();
  const XSecAlgorithmI * xsec_alg = 
     dynamic_cast<const XSecAlgorithmI *> (
          algf->GetAlgorithm(gOptGModelName, gOptGModelConf));
  if(!xsec_alg) return 0;

  double M     = kNucleonMass;
  double M2    = M*M;

  double E     = (double) kElXSecEnergy[iset];
  double theta = (double) kElXSecTheta [iset];
  double costh = TMath::Cos(2*kPi*theta/360.);

  LOG("vldtest", pNOTICE) 
     << " ** E = " << E << ", theta = " << theta 
     << " (cos(theta) = " << costh << ")";


  int target_pdgc = kElXSecTarget[iset];
  int Z = pdg::IonPdgCodeToZ(target_pdgc);
  int A = pdg::IonPdgCodeToA(target_pdgc);
  int N = A-Z;
  bool   tgt_has_p = (Z>0);
  bool   tgt_has_n = (N>0);
  double frac_p    = (double) Z / (double) A;
  double frac_n    = (double) N / (double) A;

  Interaction * ep_res = 0;
  Interaction * en_res = 0;

  if(tgt_has_p) {
     ep_res = Interaction::RESEM(1000010010, kPdgProton, kPdgElectron, E);
  }
  if(tgt_has_n) {
     en_res = Interaction::RESEM(1000000010, kPdgProton, kPdgElectron, E);
  }

  const int n = 150;
  double d2sig_dEpdOmega_array[n];
  double W2_array[n];

  double Epmin = 0.;
  double Epmax = E;
  double dEp = (Epmax-Epmin)/(n-1);

  for(int i=0; i<n; i++) {

     double Ep = Epmin + i*dEp;

     LOG("vldtest", pNOTICE) << " ** Ep = " << Ep;

     double Q2 = 2*E*Ep*(1-costh);
     double W2 = M2 + 2*M*(E-Ep)-Q2;
     double W  = TMath::Sqrt( TMath::Max(0.,W2) );

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

          LOG("vldtest", pNOTICE) 
              << "d2xsec_dWdQ2(ep;" << utils::res::AsString(kResId[ires])
              << "; E=" << E << ", W=" << W << ", Q2=" << Q2 << ") = "
              << d2sig_dWdQ2_res_p << " nbarn/GeV^3";
        }
        if(tgt_has_n) {
          en_res->ExclTagPtr()->SetResonance(kResId[ires]);
          d2sig_dWdQ2_res_n = xsec_alg->XSec(en_res,kPSWQ2fE) / units::nb;

          d2sig_dWdQ2_res_n = TMath::Max(0., d2sig_dWdQ2_res_n);

          LOG("vldtest", pNOTICE) 
              << "d2xsec_dWdQ2(en;" << utils::res::AsString(kResId[ires])
              << "; E=" << E << ", W=" << W << ", Q2=" << Q2 << ") = "
              << d2sig_dWdQ2_res_p << " nbarn/GeV^3";
        }

        d2sig_dWdQ2 += (frac_p*d2sig_dWdQ2_res_p + frac_n*d2sig_dWdQ2_res_n);
     }

     // d^2 sigma / dW dQ^2 --> d^2sigma / dE' dOmega
     double jacobian = (E*Ep)*(M+2*E*(1-costh))/(kPi*W);
     double d2sig_dEpdOmega = jacobian * d2sig_dWdQ2;

     if(TMath::IsNaN(d2sig_dEpdOmega)) {
        LOG("vldtest", pWARN) << "Got a NaN!";
        d2sig_dEpdOmega = 0;
     }

     d2sig_dEpdOmega_array[i] = TMath::Max(0., d2sig_dEpdOmega);
     W2_array[i] = W2;
  }


  if(tgt_has_p) { delete ep_res; }
  if(tgt_has_n) { delete en_res; }
  
  TGraph * gr = new TGraph(n,W2_array,d2sig_dEpdOmega_array);

  return gr;
}
//_________________________________________________________________________________
// Download cross section data from NuVld MySQL dbase 
//.................................................................................
bool Connect(void)
{
  if(!gCmpWithData) return true;

  // Get a data-base interface
  TSQLServer * sql_server = TSQLServer::Connect(
      gOptDbURL.c_str(),gOptDbUser.c_str(),gOptDbPasswd.c_str());

  if(!sql_server) return false;
  if(!sql_server->IsConnected()) return false;

  gDBI = new DBI(sql_server);
  return true;
}
//_________________________________________________________________________________
DBQ FormQuery(const char * key_list, float energy, float theta)
{
// forms a DBQueryString for extracting neutrino cross section data from the input 
// key-list and for the input energy range
//  
  ostringstream query_string;
  
  query_string 
    << "KEY-LIST:" << key_list
    << "$CUTS:E_min=" << energy-0.001 << ";E_max=" << energy+0.001 
    << ";Theta_min=" << theta-0.001 << ";Theta_max=" << theta+0.001 
    << "$DRAW_OPT:none$DB-TYPE:eN-Diff-XSec";
  
  DBQ query(query_string.str());
  
  return query;
}
//_________________________________________________________________________________
DBT * Data(int iset)
{
  if(!gCmpWithData) return 0;

  DBT * dbtable = new DBT;

  const char * keylist = kElXSecKeyList[iset];
  float        energy  = kElXSecEnergy [iset];
  float        theta   = kElXSecTheta  [iset];

  DBQ query = FormQuery(keylist, energy, theta);
  assert( gDBI->FillTable(dbtable, query) == eDbu_OK );

  return dbtable;
}
//_________________________________________________________________________________
void Draw(int iset)
{
  // get all measurements for the current channel from the NuValidator MySQL dbase
  DBT * dbtable = Data(iset);

  // get the corresponding GENIE model prediction
  TGraph * model = Model(iset,0);

  if(!model && !dbtable) return;

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

  if(dbtable) {
    TGraphAsymmErrors * graph = dbtable->GetGraph("err","W2");
    xmin  = ( graph->GetX() )[TMath::LocMin(graph->GetN(),graph->GetX())];
    xmax  = ( graph->GetX() )[TMath::LocMax(graph->GetN(),graph->GetX())];
    ymin  = ( graph->GetY() )[TMath::LocMin(graph->GetN(),graph->GetY())];
    ymax  = ( graph->GetY() )[TMath::LocMax(graph->GetN(),graph->GetY())];
    xmin  = TMath::Max(xmin, 0.5); // some data go very low 
    if(model) {
       ymin  = TMath::Min(
        ymin, ( model->GetY() )[TMath::LocMin(model->GetN(),model->GetY())]);
       ymax  = TMath::Max(
        ymax, ( model->GetY() )[TMath::LocMax(model->GetN(),model->GetY())]);
    }
    hframe = (TH1F*) gC->GetPad(1)->DrawFrame(
        scale_xmin*xmin, scale_ymin*ymin, scale_xmax*xmax, scale_ymax*ymax);
    have_frame = true;

    MultiGraph * mgraph = dbtable->GetMultiGraph("err","W2");
    for(unsigned int igraph = 0; igraph < mgraph->NGraphs(); igraph++) {
       utils::style::Format(mgraph->GetGraph(igraph), 1,1,1,1,8,0.8);
       mgraph->GetGraph(igraph)->Draw("P");
    }
  }//dbtable?

  if(model) {
    if(!have_frame) {
       xmin  = ( model->GetX() )[TMath::LocMin(model->GetN(),model->GetX())];
       xmax  = ( model->GetX() )[TMath::LocMax(model->GetN(),model->GetX())];
       ymin  = ( model->GetY() )[TMath::LocMin(model->GetN(),model->GetY())];
       ymax  = ( model->GetY() )[TMath::LocMax(model->GetN(),model->GetY())];
       hframe = (TH1F*) gC->GetPad(1)->DrawFrame(
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
    1.01*scale_ymax*ymax,kElXSecDataSetLabel[iset]);
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

  gCmpWithData = true;

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


  // get DB URL
  if(parser.OptionExists('h')) {
     gOptDbURL = parser.ArgAsString('h');
  } else {
     gOptDbURL = kDefDbURL;
  }

  // get DB username
  if(parser.OptionExists('u')) {
     gOptDbUser = parser.ArgAsString('u');
  } else {
     gCmpWithData = false;
  }

  // get DB passwd
  if(parser.OptionExists('p')) {
     gOptDbPasswd = parser.ArgAsString('p');
  } else {
     gCmpWithData = false;
  }

}
//_________________________________________________________________________________
void PrintSyntax(void)
{
  LOG("gvldtest", pNOTICE)
    << "\n\n" << "Syntax:" << "\n"
    << "   gvld_nuxsec_vs_world_data [-h host] [-u user] [-p passwd] [-m model]\n";
}
//_________________________________________________________________________________

