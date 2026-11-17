//____________________________________________________________________________
/*!

\file     Pythia8StatusCodes.h

\brief    PYTHIA 8 particle status codes used by the GENIE PYTHIA 8
          hadronization interfaces (Pythia8Hadro2019, AGCharmPythia8Hadro2023)
          when they copy the PYTHIA 8 event record into the GHEP record.

          PYTHIA 8 status codes are signed: a positive code marks a particle
          that is still present at the end of the event, and a negative code
          marks a particle that has been fragmented or decayed. The absolute
          value says how the particle was produced (see "Particle Properties"
          in the PYTHIA 8 manual):
            21-29 : particles of the hardest subprocess. GENIE appends the
                    string ends (leading quark and remnant diquark) with
                    status 23 ("outgoing"); after the string fragments PYTHIA
                    flips it to -23.
            81-89 : primary hadrons produced by the hadronization.
            91-99 : particles produced in a decay (or by Bose-Einstein
                    effects).
          So -81..-89 are primary hadrons that PYTHIA decayed during the
          hadronization (eg rho, omega, eta, K*, Delta) and -91..-99 are decay
          products that decayed in turn (eg an eta from an eta' decay).

\author   Adam Aurisano <adam.aurisano \at gmail.com>
          University of Cincinnati

\created  November 17, 2026

\cpright  Copyright (c) 2003-2026, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _PYTHIA8_STATUS_CODES_H_
#define _PYTHIA8_STATUS_CODES_H_

namespace genie {

// String ends (q, qq) after the string has been fragmented
static const int kPythia8StFragmentedStringEnd   = -23;

// Primary hadrons that PYTHIA decayed during the hadronization
static const int kPythia8StDecayedPrimaryFirst   = -81;
static const int kPythia8StDecayedPrimaryLast    = -89;

// Decay products that PYTHIA decayed in turn
static const int kPythia8StDecayedSecondaryFirst = -91;
static const int kPythia8StDecayedSecondaryLast  = -99;

// True for a hadron that PYTHIA decayed during the hadronization
inline bool IsPythia8DecayedHadron(int status)
{
  return
    (status <= kPythia8StDecayedPrimaryFirst   &&
     status >= kPythia8StDecayedPrimaryLast  ) ||
    (status <= kPythia8StDecayedSecondaryFirst &&
     status >= kPythia8StDecayedSecondaryLast);
}

} // genie namespace

#endif // _PYTHIA8_STATUS_CODES_H_
