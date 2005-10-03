//____________________________________________________________________________
/*!

\class   genie::utils::clap::CmdLineArgParserException

\brief   An exception thrown by the command line argument parser when command
         line arguments can not be found.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created October 03, 2005

*/
//____________________________________________________________________________

#include "Utils/CmdLineArgParserException.h"

using std::endl;
using namespace genie::utils::clap;

//___________________________________________________________________________
namespace genie {
 namespace utils {
  namespace clap {
   ostream & operator<< (ostream& stream, const CmdLineArgParserException & e)
   {
      e.Print(stream);
      return stream;
   }
  }
 }
}
//___________________________________________________________________________
CmdLineArgParserException::CmdLineArgParserException()
{
  this->Init();
}
//___________________________________________________________________________
CmdLineArgParserException::CmdLineArgParserException(
                                       const CmdLineArgParserException & exc)
{
  this->Copy(exc);
}
//___________________________________________________________________________
CmdLineArgParserException::~CmdLineArgParserException()
{

}
//___________________________________________________________________________
void CmdLineArgParserException::Init(void)
{
  this->fReason    = "";
  this->fArgFound  = true;

}
//___________________________________________________________________________
void CmdLineArgParserException::Copy(const CmdLineArgParserException & e)
{
  this->fReason   = e.fReason;
  this->fArgFound = e.fArgFound;
}
//___________________________________________________________________________
void CmdLineArgParserException::Print(ostream & stream) const
{
  stream << "**EXCEPTION Reason: " << this->ShowReason() << endl;
}
//___________________________________________________________________________
