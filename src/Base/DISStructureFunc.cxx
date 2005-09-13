//____________________________________________________________________________
/*!

\class    genie::DISStructureFunc

\brief    A class holding Deep Inelastic Scattering (DIS) Form Factors
          (invariant structure funstions)

          This class is using the \b Strategy Pattern. \n
          It can accept requests to calculate itself, for a given interaction,
          that it then delegates to the algorithmic object, implementing the
          DISStructureFuncModelI interface, that it finds attached to itself.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 05, 2004

*/
//____________________________________________________________________________

#include "Base/DISStructureFunc.h"

using namespace genie;

using std::endl;

//____________________________________________________________________________
namespace genie
{
  ostream & operator << (ostream & stream, const DISStructureFunc & ff)
  {
     ff.Print(stream);
     return stream;
  }
}
//____________________________________________________________________________
DISStructureFunc::DISStructureFunc()
{
  this->fModel = 0;
  this->InitFormFactors();
}
//____________________________________________________________________________
DISStructureFunc::DISStructureFunc(const DISStructureFunc & form_factors)
{
  this->fModel = form_factors.fModel;
  this->fxF1   = form_factors.fxF1;
  this->fF2    = form_factors.fF2;
  this->fxF3   = form_factors.fxF3;
  this->fF4    = form_factors.fF4;
  this->fxF5   = form_factors.fxF5;
  this->fF6    = form_factors.fF6;
}
//____________________________________________________________________________
void DISStructureFunc::SetModel(const DISStructureFuncModelI * model)
{
  this->InitFormFactors();

  this->fModel = model;
}
//____________________________________________________________________________
void DISStructureFunc::Calculate(const Interaction * interaction)
{
  fModel->CalculatePDFs(interaction);

  this->fxF1   = fModel->xF1 (interaction);
  this->fF2    = fModel->F2  (interaction);
  this->fxF3   = fModel->xF3 (interaction);
  this->fF4    = fModel->F4  (interaction);
  this->fxF5   = fModel->xF5 (interaction);
  this->fF6    = fModel->F6  (interaction);
}
//____________________________________________________________________________
void DISStructureFunc::InitFormFactors(void)
{
  this->fxF1 = 0.0;
  this->fF2  = 0.0;
  this->fxF3 = 0.0;
  this->fF4  = 0.0;
  this->fxF5 = 0.0;
  this->fF6  = 0.0;
}
//____________________________________________________________________________
void DISStructureFunc::Print(ostream & stream) const
{
  stream << "xF1  = " << this->fxF1   << endl;
  stream << "F2   = " << this->fF2    << endl;
  stream << "xF3  = " << this->fxF3   << endl;
  stream << "F4   = " << this->fF4    << endl;
  stream << "xF5  = " << this->fxF5   << endl;
  stream << "F6   = " << this->fF6    << endl;
}
//____________________________________________________________________________

