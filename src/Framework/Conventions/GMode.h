//____________________________________________________________________________
/*!

\class    genie::GEvGenMode_t

\brief    Enumeration of GENIE event generation modes

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  Nov 10, 2011

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _GENIE_MODE_H_
#define _GENIE_MODE_H_

namespace genie {

typedef enum EGEvGenMode {

  kGMdUnknown = 0,
  kGMdLeptonNucleus,     // chg.lepton/neutrino + nucleon/nucleus scattering
  kGMdHadronNucleus,     // hadron + nucleon/nucleus scattering
  kGMdPhotonNucleus,     // photon + nucleon/nucleus scattering
  kGMdDarkMatterNucleus, // dark matter + nucleon/nucleus scattering
  kGMdNucleonDecay,      // nucleon decay
  kGMdNeutronOsc         // neutron oscillation

} GEvGenMode_t;

} // genie namespace

#endif 
