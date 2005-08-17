//____________________________________________________________________________
/*!

\class   genie::EventGeneratorI

\brief   Defines the EventGeneratorI interface.

         The concrete implementations of this interface are Event Record
         Visitors (subclasses of the EventRecordVisitorI pABC) that,
         additionally, declare a 'Validity Context'. \n

         The declared validity context is used for selecting the appropriate
         concrete EventGeneratorI to generate the interacion at hand using
         the 'chain-of-responsibility' design pattern.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created November 22, 2004

*/
//____________________________________________________________________________


#ifndef _EVENT_GENERATOR_I_H_
#define _EVENT_GENERATOR_I_H_

#include "EVGCore/GVldContext.h"
#include "EVGCore/EventRecordVisitorI.h"

namespace genie {

class InteractionListGeneratorI;
class XSecAlgorithmI;

class EventGeneratorI: public EventRecordVisitorI {

//-- define an extension to the public EventRecordVisitorI interface

public :

  virtual const GVldContext &               ValidityContext  (void) const = 0;
  virtual const InteractionListGeneratorI * IntListGenerator (void) const = 0;
  virtual const XSecAlgorithmI *            CrossSectionAlg  (void) const = 0;

//-- define an extension to the private EventRecordVisitorI interface

protected:

  virtual void  InstantiateValidityContext (void) = 0;

  GVldContext * fVldContext;

  //-- dummy ctors & dtor
  EventGeneratorI();
  EventGeneratorI(const char * param_set);
  ~EventGeneratorI();
};

}      // genie namespace

#endif // _EVENT_GENERATOR_I_H_
