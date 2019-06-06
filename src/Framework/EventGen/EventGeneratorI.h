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

\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

\created November 22, 2004

\cpright Copyright (c) 2003-2019, The GENIE Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________


#ifndef _EVENT_GENERATOR_I_H_
#define _EVENT_GENERATOR_I_H_

#include "Framework/EventGen/EventRecordVisitorI.h"
#include "Framework/EventGen/GVldContext.h"

namespace genie {

class InteractionListGeneratorI;
class XSecAlgorithmI;

class EventGeneratorI: public EventRecordVisitorI {

public :

  virtual ~EventGeneratorI();

  //-- define an extension to the public EventRecordVisitorI interface
  virtual const GVldContext &               ValidityContext  (void) const = 0;
  virtual const InteractionListGeneratorI * IntListGenerator (void) const = 0;
  virtual const XSecAlgorithmI *            CrossSectionAlg  (void) const = 0;

protected:

  //-- dummy ctors & dtor
  EventGeneratorI();
  EventGeneratorI(string name);
  EventGeneratorI(string name, string config);
};

}      // genie namespace

#endif // _EVENT_GENERATOR_I_H_
