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

#ifndef INTERACTIONEXCEPTION_H_
#define INTERACTIONEXCEPTION_H_

#include <exception>
#include <iostream>
#include <string>

namespace genie
{
  namespace exceptions
  {
    class InteractionException : public std::exception
    {
      public:
        InteractionException ();
        InteractionException (const std::string & reason);
        ~InteractionException() throw() {};

        void                Print       (std::ostream & stream) const;
        const std::string & ShowReason  ()                      const { return fReason;     }

        // from std::exception
        const char * what  () const throw()  { return this->fReason.c_str(); };

        friend std::ostream & operator << (std::ostream & stream, const InteractionException & exception);

      private:
        std::string fReason;
    }; /* class InteractionException */

  } /* namespace exceptions */
} /* namespace genie */

std::ostream & operator<< (std::ostream& stream, const genie::exceptions::InteractionException & exc);

#endif /* INTERACTIONEXCEPTION_H_ */
