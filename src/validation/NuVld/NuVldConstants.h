//_____________________________________________________________________________
/*!

\namespace  genie::nuvld::constants

\brief      NuValidator constants

\author     Costas Andreopoulos (Rutherford Lab.)  <costas.andreopoulos \at stfc.ac.uk>

\created    January 12, 2004
*/
//_____________________________________________________________________________

#ifndef _NUVLD_CONSTANTS_H_
#define _NUVLD_CONSTANTS_H_

namespace genie {
namespace nuvld {
namespace constants {

static const char * kSavedPlotExtensions[] =
{
  "GIF", "*.gif", "ROOT macro", "*.C", "Encapsulated postscript", "*.eps", 0, 0
};

static const char * kSavedSplineExtensions[] =
{
  "XML", "*.xml", "ROOT", "*root", "text", "*.txt", 0, 0
};

static const char * kFitters[] =
{
   "NONE", "SIMPLE", "XSEC-NORM", "NORM-FLOAT", 0,
};

static const char * kMajorLabel =
     "GENIE Object-Oriented Neutrino MC Generator";
static const char * kMinorLabel =
     "http://www.genie-mc.org";


} // namespace constants
} // namespace nuvld
} // namespace genie

#endif


