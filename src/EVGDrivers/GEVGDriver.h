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

using std::ostream;
using std::string;

namespace genie {

class EventRecord;
class EventGeneratorList;
class InteractionFilter;
class InteractionSelectorI;
class EGResponsibilityChain;
class InitialState;
class Target;
class Spline;

class GEVGDriver {

public :

  GEVGDriver();
  ~GEVGDriver();

  //-- Set initial state, interaction filter and option to create splines
  void SetInitialState   (int nu_pdgc, int Z, int A);
  void SetInitialState   (const InitialState & init_state);
  void SetFilter         (const InteractionFilter & filter);
  void UseSplines        (void);
  void CreateSplines     (bool useLogE = true);
  void FilterUnphysical  (bool on_off);

  //-- Generate single event
  EventRecord * GenerateEvent (const TLorentzVector & nu4p);

  //-- Cross section sum for all interactions that can be generated for
  //   the current init-state, and the cross section for the maximum
  //   cross section interaction that can be generated for the current
  //   init-state.
  //   The cross sections for specific interactions and the interaction
  //   list for the current init-state is accessed from the EventGeneratorList
  //   object).
  //   These methods are used by the GENIE MC job driver to select a target
  //   nucleus from a ROOT geometry and estimate the maximum interaction
  //   probability that scales all interaction probabilities.
  double SumCrossSection(const TLorentzVector & nup4);
  double MaxCrossSection(const TLorentzVector & nup4);
  void   CreateXSecSumSpline(int nk, double Emin, double Emax, bool inlogE=true);
  const Spline * XSecSumSpline(void) const { return fXSecSumSpl; }

  //-- Get state
  const EventGeneratorList * EventGenerators (void) const { return fEvGenList; }
  const InteractionFilter *  Filter          (void) const { return fFilter;    }

  //-- Print state
  void Print (ostream & stream) const;
  friend ostream & operator << (ostream & stream, const GEVGDriver & driver);

private:

  //-- Private initialization, configuration & input validation methods
  void Initialize       (void);
  void Configure        (void);
  bool IsValidInitState (void) const;

  //-- Minimal Initial State Information
  bool     fUseSplines;
  int      fNuPDG;
  Target * fNuclTarget;

  //-- Private data members
  EventRecord *           fCurrentRecord;   ///< ptr to the event record being processed
  EventGeneratorList *    fEvGenList;       ///< all Event Generators available at this job
  EGResponsibilityChain * fChain;           ///< an Event Generator chain of responsibility
  InteractionSelectorI *  fIntSelector;     ///< interaction selector
  InteractionFilter *     fFilter;          ///< interaction filter
  bool                    fFilterUnphysical;///< controls whether unphysical events are returned
  Spline *                fXSecSumSpl;      ///< sum{xsec(all interactions | this init state)}
  unsigned int            fNRecLevel;       ///< recursive mode depth counter
};

}      // genie namespace

#endif // _GENIE_H_
