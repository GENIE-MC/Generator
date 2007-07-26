//____________________________________________________________________________
/*!

\class   genie::GEVGDriver

\brief   Minimal interface object for generating neutrino interactions for
         a given initial state.

         If you need to generate events for a given neutrino flux and detector 
         geometry (and therefore for a multitude of possible initial states) 
         then use the GMCJDriver.
         It is worth noting that GEVGDriver is the piece of code that puts the
         actual event generation framework into motion, and GMCDriver itself
         assembles a list of GEVGDrivers (1 / possible initial state).

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         STFC, Rutherford Appleton Laboratory

\created August 06, 2004

\cpright  Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
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

#include "Utils/Range1.h"

using std::ostream;
using std::string;

namespace genie {

class EventRecord;
class EventGeneratorList;
class InteractionSelectorI;
class InteractionGeneratorMap;
class InteractionList;
class Interaction;
class InitialState;
class Target;
class Spline;

class GEVGDriver {

public :
  GEVGDriver();
  ~GEVGDriver();

  //-- Set driver options before calling Configure()
  void FilterUnphysical (const TBits & unphysmask);
  void UseSplines       (void);

  //-- Configure the driver
  void Configure (int nu_pdgc, int Z, int A);
  void Configure (const InitialState & init_state);

  //-- Generate single event
  EventRecord * GenerateEvent (const TLorentzVector & nu4p);

  //-- Get loaded event generator list
  const EventGeneratorList * EventGenerators (void) const { return fEvGenList; }

  //-- Get the list of all interactions that can be simulated for the specified 
  //   initial state (depends on which event generation threads were loaded into
  //   the event generation driver driver)
  const InteractionList * Interactions(void) const;

  //-- Cross section splines for input interaction and for the sum of all
  //   simulated interactions for the specified initial state
  const Spline * XSecSumSpline       (void) const { return fXSecSumSpl; }
  const Spline * XSecSpline          (const Interaction * interaction) const;

  //-- Instruct the driver to create all the splines it needs
  void CreateSplines (int nknots=-1, double emax=-1, bool inLogE=true);

  //-- Methods used for building the 'total' cross section spline
  double XSecSum             (const TLorentzVector & nup4);
  void   CreateXSecSumSpline (int nk, double Emin, double Emax, bool inlogE=true);

  //-- Get validity range (combined validity range of loaded evg threads)
  Range1D_t ValidEnergyRange (void) const;

  //-- Reset, Print etc
  void Reset (void);
  void Print (ostream & stream) const;

  friend ostream & operator << (ostream & stream, const GEVGDriver & driver);

private:

  //-- Private initialization, configuration & input validation methods
  void Init                         (void);
  void CleanUp                      (void);
  void BuildInitialState            (const InitialState & init_state);
  void BuildGeneratorList           (void);
  void BuildInteractionGeneratorMap (void);
  void BuildInteractionSelector     (void);
  void AssertIsValidInitState       (void) const;

  //-- Private data members
  InitialState *            fInitState;       ///< initial state information for driver instance
  EventRecord *             fCurrentRecord;   ///< ptr to the event record being processed
  EventGeneratorList *      fEvGenList;       ///< all Event Generators available at this job
  InteractionSelectorI *    fIntSelector;     ///< interaction selector
  InteractionGeneratorMap * fIntGenMap;       ///< interaction -> generator assosiative container
  TBits *                   fFiltUnphysMask;  ///< controls whether unphysical events are returned
  bool                      fUseSplines;      ///< controls whether xsecs are computed or interpolated
  Spline *                  fXSecSumSpl;      ///< sum{xsec(all interactions | this init state)}
  unsigned int              fNRecLevel;       ///< recursive mode depth counter
};

}      // genie namespace

#endif // _GENIE_H_
