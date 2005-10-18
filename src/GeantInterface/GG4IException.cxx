//____________________________________________________________________________
/*!

\class   genie::exceptions::GG4IException

\brief   GENIE / GEANT4 Interface Exception

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created October 03, 2005

*/
//____________________________________________________________________________

#include "GeantInterface/GG4IException.h"

using std::endl;
using namespace genie::exceptions;

//___________________________________________________________________________
namespace genie {
 namespace exceptions {
   ostream & operator<< (ostream& stream, const GG4IException & e)
   {
      e.Print(stream);
      return stream;
   }
 }
}
//___________________________________________________________________________
GG4IException::GG4IException()
{
  this->Init();
}
//___________________________________________________________________________
GG4IException::GG4IException(const GG4IException & exc)
{
  this->Copy(exc);
}
//___________________________________________________________________________
GG4IException::~GG4IException()
{

}
//___________________________________________________________________________
void GG4IException::Init(void)
{
  this->fReason       = "";
  this->fFileNotFound = false;
  this->fEmptyTree    = false;
  this->fEndOfFile    = false;
  this->fFormatError  = false;
}
//___________________________________________________________________________
void GG4IException::Copy(const GG4IException & e)
{
  this->fReason       = e.fReason;
  this->fFileNotFound = e.fFileNotFound;
  this->fEmptyTree    = e.fEmptyTree;
  this->fEndOfFile    = e.fEndOfFile;
  this->fFormatError  = e.fFormatError;
}
//___________________________________________________________________________
void GG4IException::Print(ostream & stream) const
{
  stream << "**EXCEPTION Reason: " << this->ShowReason() << endl;
}
//___________________________________________________________________________
