//____________________________________________________________________________
/*!

\class    genie::QELFormFactors

\brief    A class holding Quasi Elastic (QEL) Form Factors.

          This class is using the \b Strategy Pattern. \n
          It can accept requests to calculate itself, for a given interaction,
          that it then delegates to the algorithmic object, implementing the
          QELFormFactorsModelI interface, that it finds attached to itself.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  April 20, 2004

*/
//____________________________________________________________________________

#include "Base/QELFormFactors.h"
#include "Messenger/Messenger.h"

using namespace genie;

using std::cout;
using std::endl;

//____________________________________________________________________________
namespace genie
{
  ostream & operator << (ostream & stream, const QELFormFactors & ff)
  {
     ff.Print(stream);     
     return stream;
  }  
}
//____________________________________________________________________________
QELFormFactors::QELFormFactors()
{
  fModel = 0;
  InitFormFactors();
}
//____________________________________________________________________________
QELFormFactors::QELFormFactors(const QELFormFactors & form_factors)
{
  fModel = form_factors.fModel;

  fF1V   = form_factors.fF1V;
  fxiF2V = form_factors.fxiF2V;
  fFA    = form_factors.fFA;
  fFp    = form_factors.fFp;
}
//____________________________________________________________________________
void QELFormFactors::SetModel(const QELFormFactorsModelI * model)
{
  InitFormFactors();

  fModel = model;
}
//____________________________________________________________________________
void QELFormFactors::Calculate(const Interaction * interaction)
{
  fF1V   = fModel -> F1V   (interaction);
  fxiF2V = fModel -> xiF2V (interaction);
  fFA    = fModel -> FA    (interaction);
  fFp    = fModel -> Fp    (interaction);
}
//____________________________________________________________________________
void QELFormFactors::InitFormFactors()
{
  fF1V   = 0;
  fxiF2V = 0;
  fFA    = 0;
  fFp    = 0;
}
//____________________________________________________________________________
void QELFormFactors::Print(ostream & stream) const
{
  stream << endl;
  stream << "F1V    = " << fF1V   << endl;
  stream << "xi*F2V = " << fxiF2V << endl;
  stream << "FA     = " << fFA    << endl;
  stream << "Fp     = " << fFp    << endl;
}
//____________________________________________________________________________

