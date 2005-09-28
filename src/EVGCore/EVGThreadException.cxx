//____________________________________________________________________________
/*!

\class   genie::exceptions::EVGThreadException

\brief   An exception thrown by EventRecordVisitorI when the normal processing
         sequence has to be disrupted (fast-fwd at the end or step-back)

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created September 27, 2005

*/
//____________________________________________________________________________

#include "EVGCore/EVGThreadException.h"
#include "Messenger/Messenger.h"

using std::endl;
using namespace genie::exceptions;

//___________________________________________________________________________
namespace genie {
 namespace exceptions {
  ostream & operator<< (ostream& stream, const EVGThreadException & exc)
  {
   exc.Print(stream);
   return stream;
  }
 }
}
//___________________________________________________________________________
EVGThreadException::EVGThreadException()
{
  this->Init();
}
//___________________________________________________________________________
EVGThreadException::EVGThreadException(const EVGThreadException & exc)
{
  this->Copy(exc);
}
//___________________________________________________________________________
EVGThreadException::~EVGThreadException()
{

}
//___________________________________________________________________________
void EVGThreadException::Init(void)
{
  fReason     = "";
  fFastFwd    = false;
  fStepBack   = false;
  fReturnStep = 999999;
}
//___________________________________________________________________________
void EVGThreadException::Copy(const EVGThreadException & exc)
{
  fReason     = exc.fReason;
  fFastFwd    = exc.fFastFwd;
  fStepBack   = exc.fStepBack;
  fReturnStep = exc.fReturnStep;
}
//___________________________________________________________________________
void EVGThreadException::Print(ostream & stream) const
{
  stream << "**EXCEPTION Reason: " << this->ShowReason() << endl;
}
//___________________________________________________________________________
