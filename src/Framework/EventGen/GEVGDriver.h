//____________________________________________________________________________
/*!

\class   genie::GEVGDriver

\brief   GENIE Event Generation Driver.
         A minimalist user interface object for generating neutrino interactions.
         Each such object is configured for a given initial state and it drives all
         relevant GENIE neutrino interaction physics simulation code for that state.
         To set-up MC jobs involving a multitude of possible initial states,
         including arbitrarily complex neutrino flux and detector geometry 
         descriptions, see the GMCJDriver object.

\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

\created August 06, 2004

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _GEVG_DRIVER_H_
#define _GEVG_DRIVER_H_

#include <ostream>
#include <string>

#include <TLorentzVector.h>
#include <TBits.h>

#include "Framework/Utils/Range1.h"

using std::ostream;
using std::string;

namespace genie {

class GEVGDriver;
class EventRecord;
class EventGeneratorList;
class EventGeneratorI;
class InteractionSelectorI;
class InteractionGeneratorMap;
class InteractionList;
class Interaction;
class InitialState;
class Target;
class Spline;

ostream & operator << (ostream & stream, const GEVGDriver & driver);

class GEVGDriver {

public :
  GEVGDriver();
 ~GEVGDriver();

  // Driver options:
  // - Set before calling Configure()
  void UseSplines (void);
  void SetEventGeneratorList(string listname);
  // - Set before GenerateEvent()
  void SetUnphysEventMask(const TBits & mask);

  // Configure the driver
  void Configure (int nu_pdgc, int Z, int A);
  void Configure (const InitialState & init_state);

  // Generate single event
  EventRecord * GenerateEvent (const TLorentzVector & nu4p);

  // Get the list of all interactions that can be simulated for the specified 
  // initial state (depends on which event generation threads were loaded into
  // the event generation driver driver)
  const InteractionList * Interactions(void) const;

  // Get event generator thread list
  const EventGeneratorList * EventGenerators (void) const { return fEvGenList; }

  // Get the event generator that is responsible for generating the input event
  const EventGeneratorI * FindGenerator(const Interaction * interaction) const;

  // Cross section splines for input interaction and for the sum of all
  // simulated interactions for the specified initial state
  const Spline * XSecSumSpline       (void) const { return fXSecSumSpl; }
  const Spline * XSecSpline          (const Interaction * interaction) const;

  // Instruct the driver to create all the splines it needs
  void CreateSplines (int nknots=-1, double emax=-1, bool inLogE=true);

  // Methods used for building the 'total' cross section spline
  double XSecSum             (const TLorentzVector & nup4);
  void   CreateXSecSumSpline (int nk, double Emin, double Emax, bool inlogE=true);

  // Get validity range (combined validity range of loaded evg threads)
  Range1D_t ValidEnergyRange (void) const;

  // Reset, Print etc
  void Reset (void);
  void Print (ostream & stream) const;

  friend ostream & operator << (ostream & stream, const GEVGDriver & driver);

private:

  // Private initialization, configuration & input validation methods
  void Init                         (void);
  void CleanUp                      (void);
  void BuildInitialState            (const InitialState & init_state);
  void BuildGeneratorList           (void);
  void BuildInteractionGeneratorMap (void);
  void BuildInteractionSelector     (void);
  void AssertIsValidInitState       (void) const;

  // Private data members
  InitialState *            fInitState;       ///< initial state information for driver instance
  EventRecord *             fCurrentRecord;   ///< ptr to the event record being processed
  EventGeneratorList *      fEvGenList;       ///< all Event Generators available at this job
  InteractionSelectorI *    fIntSelector;     ///< interaction selector
  InteractionGeneratorMap * fIntGenMap;       ///< interaction -> generator assosiative container
  TBits *                   fUnphysEventMask; ///< controls whether unphysical events are returned
  bool                      fUseSplines;      ///< controls whether xsecs are computed or interpolated
  Spline *                  fXSecSumSpl;      ///< sum{xsec(all interactions | this init state)}
  unsigned int              fNRecLevel;       ///< recursive mode depth counter
  string                    fEventGenList;    ///< list of event generators loaded by this driver (what used to be the $GEVGL setting)
};

}      // genie namespace

#endif // _GENIE_H_
