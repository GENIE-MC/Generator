//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

         Steve Dytman <dytman \at pitt.edu>
	 Univ. of Pittsburgh         

\created October 10, 2011

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include "Physics/HadronTransport/INukeException.h"
#include "Framework/Messenger/Messenger.h"

using std::endl;
using namespace genie::exceptions;

//___________________________________________________________________________
namespace genie {
 namespace exceptions {
  ostream & operator<< (ostream& stream, const INukeException & exc)
  {
   exc.Print(stream);
   return stream;
  }
 }
}
//___________________________________________________________________________
INukeException::INukeException()
{
  this->Init();
}
//___________________________________________________________________________
INukeException::INukeException(const INukeException & exc)
{
  this->Copy(exc);
}
//___________________________________________________________________________
INukeException::~INukeException()
{

}
//___________________________________________________________________________
void INukeException::Init(void)
{
  fReason = "";
}
//___________________________________________________________________________
void INukeException::Copy(const INukeException & exc)
{
  fReason = exc.fReason;
}
//___________________________________________________________________________
void INukeException::Print(ostream & stream) const
{
  stream << "**EXCEPTION Reason: " << this->ShowReason() << endl;
}
//___________________________________________________________________________
