//____________________________________________________________________________
/*
 Copyright (c) 2003-2011, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory

         Steve Dytman <dytman \at pitt.edu>
	 Univ. of Pittsburgh         

\created October 10, 2011

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include "HadronTransport/INukeException_inel.h"
#include "Messenger/Messenger.h"

using std::endl;
using namespace genie::exceptions;

//___________________________________________________________________________
namespace genie {
 namespace exceptions {
  ostream & operator<< (ostream& stream, const INukeException_inel & exc)
  {
   exc.Print(stream);
   return stream;
  }
 }
}
//___________________________________________________________________________
INukeException_inel::INukeException_inel()
{
  this->Init();
}
//___________________________________________________________________________
INukeException_inel::INukeException_inel(const INukeException_inel & exc)
{
  this->Copy(exc);
}
//___________________________________________________________________________
INukeException_inel::~INukeException_inel()
{

}
//___________________________________________________________________________
void INukeException_inel::Init(void)
{
  fReason = "";
}
//___________________________________________________________________________
void INukeException_inel::Copy(const INukeException_inel & exc)
{
  fReason = exc.fReason;
}
//___________________________________________________________________________
void INukeException_inel::Print(ostream & stream) const
{
  stream << "**EXCEPTION Reason: " << this->ShowReason() << endl;
}
//___________________________________________________________________________
