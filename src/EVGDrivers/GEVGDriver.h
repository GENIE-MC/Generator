//____________________________________________________________________________
/*!

\class   genie::GEVGDriver

\brief   Minimal interface object for generating neutrino interactions for
         a given initial state.

         When GMC is used, a GEVGDriver object list is assembled for all possible
         initial states (corresponding to combinations of all neutrino types
         -declared by the input GFluxI- and all target nuclei types -found
         in the input ROOT geometry-.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created August 06, 2004

*/
//____________________________________________________________________________

#ifndef _GEVG_DRIVER_H_
#define _GEVG_DRIVER_H_

#include <ostream>
#include <string>

#include <TLorentzVector.h>
#include "Utils/Range1.h"

using std::ostream;
using std::string;

namespace genie {

class EventRecord;
class EventGeneratorList;
class InteractionSelectorI;
class EGResponsibilityChain;
class InitialState;
class Target;
class Spline;
class XSecAlgorithmMap;

class GEVGDriver {

public :

  GEVGDriver();
  ~GEVGDriver();

  //! Set driver options (before calling Configure())
  void FilterUnphysical       (bool on_off);
  void UseSplines             (void);
  //void UseInteractionSelector (string name, string config);

  //! Configure the driver
  void Configure (int nu_pdgc, int Z, int A);
  void Configure (const InitialState & init_state);

  //! Generate single event
  EventRecord * GenerateEvent (const TLorentzVector & nu4p);

  //! Instruct the driver to create all the splines it needs
  void CreateSplines (bool useLogE = true);

  //! Cross section sum for all interactions that can be generated for
  //! the current init-state.
  double         XSecSum             (const TLorentzVector & nup4);
  void           CreateXSecSumSpline (int nk, double Emin, double Emax, bool inlogE=true);
  const Spline * XSecSumSpline       (void) const { return fXSecSumSpl; }

  //! Get loaded event generator list
  const EventGeneratorList * EventGenerators (void) const { return fEvGenList; }

  //! Get validity range (combined validity range of loaded evg threads)
  Range1D_t ValidEnergyRange (void) const;

  //! Reset, Print etc
  void Reset (void);
  void Print (ostream & stream) const;

  friend ostream & operator << (ostream & stream, const GEVGDriver & driver);

private:

  //! Private initialization, configuration & input validation methods
  void Init                     (void);
  void CleanUp                  (void);
  void BuildInitialState        (const InitialState & init_state);
  void BuildGeneratorList       (void);
  void BuildXSecAlgorithmMap    (void);
  void BuildResponsibilityChain (void);
  void BuildInteractionSelector (void);
  void AssertIsValidInitState   (void) const;

  //! Private data members
  InitialState *          fInitState;       ///< initial state information for driver instance
  EventRecord *           fCurrentRecord;   ///< ptr to the event record being processed
  EventGeneratorList *    fEvGenList;       ///< all Event Generators available at this job
  EGResponsibilityChain * fChain;           ///< an Event Generator chain of responsibility
  InteractionSelectorI *  fIntSelector;     ///< interaction selector
  XSecAlgorithmMap *      fXSecAlgorithmMap;///< interaction -> xsec alg. assosiative container
  bool                    fFilterUnphysical;///< controls whether unphysical events are returned
  bool                    fUseSplines;      ///< controls whether xsecs are computed or interpolated
  Spline *                fXSecSumSpl;      ///< sum{xsec(all interactions | this init state)}
  unsigned int            fNRecLevel;       ///< recursive mode depth counter
};

}      // genie namespace

#endif // _GENIE_H_
