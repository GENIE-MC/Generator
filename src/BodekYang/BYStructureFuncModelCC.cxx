//____________________________________________________________________________
/*!

\class    genie::BYStructureFuncModelCC

\brief    Computes CC vN DIS Form Factors according to the Bodek-Yang model.
          Inherits part of its implemenation from the BYStructureFuncModel
          abstract class.

          Check out BYStructureFuncModel for comments and references.

          Is a concrete implementation of the DISFormFactorsModelI interface.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  September 28, 2004

*/
//____________________________________________________________________________

#include "BodekYang/BYStructureFuncModelCC.h"
#include "Conventions/Constants.h"
#include "Messenger/Messenger.h"
#include "PDF/PDF.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
BYStructureFuncModelCC::BYStructureFuncModelCC() :
BYStructureFuncModel("genie::BYStructureFuncModelCC")
{

}
//____________________________________________________________________________
BYStructureFuncModelCC::BYStructureFuncModelCC(string config):
BYStructureFuncModel("genie::BYStructureFuncModelCC", config)
{

}
//____________________________________________________________________________
BYStructureFuncModelCC::~BYStructureFuncModelCC()
{

}
//____________________________________________________________________________
