//____________________________________________________________________________
/*!

\class    genie::exceptions::InteractionException

\brief    Exception used inside Interaction classes.

\author   Jeremy Wolcott <jwolcott \at fnal.gov>
          Tufts University

\created  July 15, 2016

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________


#include <iostream>

#include "InteractionException.h"

namespace genie
{
  namespace exceptions
  {

    InteractionException::InteractionException ()
      : fReason("")
    {}

    InteractionException::InteractionException (const std::string& reason)
      : fReason(reason)
    {}

    void InteractionException::Print (std::ostream& stream) const
    {
      stream << "**EXCEPTION Reason: " << this->ShowReason() << std::endl;
    }

  } /* namespace exceptions */
} /* namespace genie */

std::ostream & operator<< (std::ostream& stream, const genie::exceptions::InteractionException & exc)
{
 exc.Print(stream);
 return stream;
}
