//____________________________________________________________________________
/*!

\class   genie::geometry::PointGeomAnalyzer

\brief   The PointGeomAnalyzer class is the simplest implementation of the
         GeomAnalyserI interface and defines a simple 'point-like' geometry.

         Use this geometry analyzer to generate events when you do not want
         to use a detailed GEANT/ROOT geometry description but you only need
         to generate events for a 'single' nuclear target while you still want
         to use the GENIE MC job driver 'loaded' with a GENIE flux driver.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created July 14, 2005

*/
//____________________________________________________________________________

#include <TLorentzVector.h>
#include <TVector3.h>

#include "Geo/PointGeomAnalyzer.h"
#include "EVGDrivers/PathLengthList.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodeList.h"
#include "PDG/PDGLibrary.h"
#include "Utils/PrintUtils.h"

using namespace genie;
using namespace genie::geometry;

//___________________________________________________________________________
PointGeomAnalyzer::PointGeomAnalyzer(int tgtpdgc) :
GeomAnalyzerI()
{
  fCurrVertex         = new TVector3(0,0,0);
  fCurrPathLengthList = new PathLengthList(1);
  fCurrPDGCodeList    = new PDGCodeList(1);

  (*fCurrPDGCodeList)[0]    = tgtpdgc;
  (*fCurrPathLengthList)[0] = 1.;
}
//___________________________________________________________________________
PointGeomAnalyzer::~PointGeomAnalyzer()
{
  if( fCurrPathLengthList ) delete fCurrPathLengthList;
  if( fCurrPDGCodeList    ) delete fCurrPDGCodeList;
}
//___________________________________________________________________________
const PDGCodeList & PointGeomAnalyzer::ListOfTargetNuclei(void)
{
// pdg code list contains a single code corresponding to the material passed
// at the geom analyser ctor

  return *fCurrPDGCodeList;
}
//___________________________________________________________________________
const PathLengthList & PointGeomAnalyzer::ComputeMaxPathLengths(void)
{
// this is irrelevant for the 'point' geometry - return a path length of 1.
// for the only defined material

  return *fCurrPathLengthList;
}
//___________________________________________________________________________
const PathLengthList & PointGeomAnalyzer::ComputePathLengths(
                  const TLorentzVector & /*x*/, const TLorentzVector & /*p*/)
{
// this is irrelevant for the 'point' geometry - return a path length of 1.
// for the only defined material

  return *fCurrPathLengthList;
}
//___________________________________________________________________________
const TVector3 & PointGeomAnalyzer::GenerateVertex(
  const TLorentzVector & /*x*/, const TLorentzVector & /*p*/, int /*tgtpdg*/)
{
// this is irrelevant for the 'point' geometry - return a vtx at (0,0,0)

  return *fCurrVertex;
}
//___________________________________________________________________________
