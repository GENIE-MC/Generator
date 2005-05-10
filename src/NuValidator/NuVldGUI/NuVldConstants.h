//_____________________________________________________________________________
/*!

\namespace  genie::nuvld::constants

\brief      NuValidator constants

\author     Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>

\created    January 12, 2004
*/
//_____________________________________________________________________________

#ifndef _NUVLD_CONSTANTS_H_
#define _NUVLD_CONSTANTS_H_

namespace genie {
namespace nuvld {
namespace constants {

//static const char * kSqlFileExtensions[] =
//{
// "All files", "*", "SQL files", "*.sql", 0, 0
//};

//static const char * kRootFileExtensions[] =
//{
//  "All files", "*", "ROOT files", "*.root", 0, 0
//};

static const char * kSavedPlotExtensions[] =
{
  "GIF", "*.gif", "ROOT macro", "*.C", "Encapsulated postscript", "*.eps", 0, 0
};

static const char * kFitters[] =
{
   "NONE", "SIMPLE", "NORM-FLOAT", 0, 
};

static const char * kMajorLabel =
     "GENIE Universal Object-Oriented Neutrino Generator Collaboration";
static const char * kMinorLabel =
                   "http://hepunx.rl.ac.uk/~candreop/generators/GENIE/";


} // namespace constants
} // namespace nuvld
} // namespace genie

#endif


