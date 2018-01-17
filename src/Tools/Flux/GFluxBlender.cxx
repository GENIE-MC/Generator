////////////////////////////////////////////////////////////////////////
/// \file  GFluxBlender.h
/// \brief GENIE GFluxI adapter to allow flavor modification
///
/// \version $Id: GFluxBlender.cxx,v 1.1.1.1 2010/12/22 16:18:52 p-nusoftart Exp $
/// \author  Robert Hatcher <rhatcher \at fnal.gov>
///          Fermi National Accelerator Laboratory
///
/// \update  2010-10-31 initial version
///
/// \update  2011-02-22 - JD
///   Implemented dummy versions of the new GFluxI::Clear, GFluxI::Index 
///   and GFluxI::GenerateWeighted methods needed for pre-generation of 
///   flux interaction probabilities in GMCJDriver.
///
////////////////////////////////////////////////////////////////////////
#include <math.h>
#include <iostream>
#include <iomanip>

//GENIE includes
#include "PDG/PDGCodes.h"
#include "FluxDrivers/GNuMIFlux.h"
#include "FluxDrivers/GSimpleNtpFlux.h"
#include "Numerical/RandomGen.h"

#include "FluxDrivers/GFluxBlender.h"
#include "FluxDrivers/GFlavorMixerI.h"
#include "Messenger/Messenger.h"
#define  LOG_BEGIN(a,b)   LOG(a,b)
#define  LOG_END ""

namespace genie {
namespace flux {

//____________________________________________________________________________

GFluxBlender::GFluxBlender() : 
  GFluxI(),
  fRealGFluxI(0),
  fGNuMIFlux(0),
  fGSimpleFlux(0),
  fFlavorMixer(0),
  fPDGListGenerator(),
  fPDGListMixed(),
  fNPDGOut(0),
  fBaselineDist(0),
  fEnergy(0),
  fDistance(0),
  fPdgCGenerated(0),
  fPdgCMixed(0)
{ ; }

GFluxBlender::~GFluxBlender()
{
  if ( fRealGFluxI  ) { delete fRealGFluxI;  fRealGFluxI  = 0; }
  if ( fFlavorMixer ) { delete fFlavorMixer; fFlavorMixer = 0; }
}

//____________________________________________________________________________
const PDGCodeList& GFluxBlender::FluxParticles(void)
{
  /// Return (reference to) list of possible neutrinos *after* mixing.
  /// These are PDG code that might be interacted.

  // this is only ever called by GMCJDriver::GetParticleLists()
  // during GMCJDriver::Configure() which should happen just once
  // so don't try to be too clever.
  fPDGListGenerator = fRealGFluxI->FluxParticles();

  // okay, really stupid at this time ...
  fPDGListMixed.clear();
  fPDGListMixed.push_back(kPdgNuE);
  fPDGListMixed.push_back(kPdgNuMu);
  fPDGListMixed.push_back(kPdgNuTau);
  fPDGListMixed.push_back(kPdgAntiNuE);
  fPDGListMixed.push_back(kPdgAntiNuMu);
  fPDGListMixed.push_back(kPdgAntiNuTau);

  // size the probability arrays to the same number of entries
  fNPDGOut = fPDGListMixed.size();
  fProb.resize(fNPDGOut);
  fSumProb.resize(fNPDGOut);

  if ( ! fFlavorMixer ) return fRealGFluxI->FluxParticles();
  else                  return fPDGListMixed;
}

//____________________________________________________________________________
bool GFluxBlender::GenerateNext(void)
{
  
  bool gen1 = false;
  while ( ! gen1 ) {
    if ( ! fRealGFluxI->GenerateNext() ) return false;
    // have a new entry
    fPdgCGenerated = fRealGFluxI->PdgCode();
    if ( ! fFlavorMixer ) {
      // simple case when not configured with a mixing model
      fPdgCMixed = fPdgCGenerated;
      gen1 = true;
    } else {
      // now pick a new flavor 
      fDistance = fBaselineDist; 
      if ( fGNuMIFlux   ) fDistance = fGNuMIFlux->GetDecayDist();
      if ( fGSimpleFlux ) fDistance = fGSimpleFlux->GetDecayDist();
      fEnergy = fRealGFluxI->Momentum().Energy();
      fPdgCMixed = ChooseFlavor(fPdgCGenerated,fEnergy,fDistance);
      // we may have to generate a new neutrino if it oscillates away
      // don't pass non-particles to GENIE ... it won't like it
      gen1 = ( fPdgCMixed != 0 );
    }
  }
  return true;
}
//____________________________________________________________________________
void GFluxBlender::Clear(Option_t * opt)
{
// Clear method needed to conform to GFluxI interface 
//
  fRealGFluxI->Clear(opt);
}
//____________________________________________________________________________
long int GFluxBlender::Index(void)
{
// Index method needed to conform to GFluxI interface 
//
  return fRealGFluxI->Index();
}
//____________________________________________________________________________
void GFluxBlender::GenerateWeighted(bool gen_weighted)
{
// Dummy implementation needed to conform to GFluxI interface
//
  LOG("FluxBlender", pERROR) <<
      "No GenerateWeighted(bool gen_weighted) method implemented for " <<
      "gen_weighted: " << gen_weighted;
}
//____________________________________________________________________________
GFluxI* GFluxBlender::AdoptFluxGenerator(GFluxI* generator)
{
  GFluxI* oldgen = fRealGFluxI;
  fRealGFluxI  = generator;
  // avoid re-casting
  fGNuMIFlux   = dynamic_cast<GNuMIFlux*>(fRealGFluxI);
  fGSimpleFlux = dynamic_cast<GSimpleNtpFlux*>(fRealGFluxI);
  // force evaluation of particle lists
  this->FluxParticles();

  return oldgen;
}

//____________________________________________________________________________
GFlavorMixerI* GFluxBlender::AdoptFlavorMixer(GFlavorMixerI* mixer)
{
  GFlavorMixerI* oldmix = fFlavorMixer;
  fFlavorMixer = mixer;
  return oldmix;
}

//____________________________________________________________________________
int GFluxBlender::ChooseFlavor(int pdg_init, double energy, double dist)
{
  // choose a new flavor
  bool   isset = false;
  int    pdg_out = 0;
  double sumprob = 0;
    
  fRndm = RandomGen::Instance()->RndFlux().Rndm();
  for (size_t indx = 0; indx < fNPDGOut; ++indx ) {
    int pdg_test = fPDGListMixed[indx];
    double prob = fFlavorMixer->Probability(pdg_init,pdg_test,energy,dist);
    fProb[indx] = prob;
    sumprob += fProb[indx];
    fSumProb[indx] = sumprob;
    if ( ! isset && fRndm < sumprob ) {
      isset   = true;
      pdg_out = pdg_test;
    }
  }

  return pdg_out;
}

//____________________________________________________________________________
void GFluxBlender::PrintConfig(void)
{
  LOG_BEGIN("FluxBlender", pINFO) << "GFluxBlender::PrintConfig()" << LOG_END;
  if ( fRealGFluxI ) {
    LOG_BEGIN("FluxBlender", pINFO) 
      << "   fRealGFluxI is a \"" 
      << typeid(fRealGFluxI).name() << "\"" 
      << LOG_END;
  } else {
    LOG_BEGIN("FluxBlender", pINFO) 
      << "   fRealGFluxI is not initialized" << LOG_END;
  }
  if ( fFlavorMixer ) {
    LOG_BEGIN("FluxBlender", pINFO) 
      << "   fFlavorMixer is a \""
      << typeid(fFlavorMixer).name() << "\"" 
      << LOG_END;
  } else {
    LOG_BEGIN("FluxBlender", pINFO)
      << "   fFlavorMixer is not initialized" << LOG_END;
  }
  LOG_BEGIN("FluxBlender", pINFO) 
    << "   BaselineDist " << fBaselineDist << LOG_END;
  LOG_BEGIN("FluxBlender", pINFO) 
    << "PDG List from Generator" << fPDGListGenerator << LOG_END;
  LOG_BEGIN("FluxBlender", pINFO)
    << "PDG List after mixing (n=" << fNPDGOut << ")"
    << fPDGListMixed << LOG_END;

}

//____________________________________________________________________________
void GFluxBlender::PrintState(bool verbose)
{
  LOG_BEGIN("FluxBlender", pINFO) << "GFluxBlender::PrintState()" << LOG_END;
  LOG_BEGIN("FluxBlender", pINFO) 
    << "  Flavor " << fPdgCGenerated 
    << " ==> " << fPdgCMixed 
    << " (E=" << fEnergy << ", dist=" << fDistance << ")" << LOG_END;
  if ( verbose ) {
    LOG_BEGIN("FluxBlender", pINFO) << "  Rndm = " << fRndm << LOG_END;
    for (size_t indx = 0; indx < fNPDGOut; ++indx ) {
      LOG_BEGIN("FluxBlender", pINFO)
        << "   [" << indx << "] "
        << std::setw(3) << fPdgCGenerated << " => "
        << std::setw(3) << fPDGListMixed[indx]
        << "  p = " << std::setw(10) << fProb[indx]
        << "  sum_p = " << std::setw(10) << fSumProb[indx] << LOG_END;
    }
  }
}

//____________________________________________________________________________
} // namespace flux
} // namespace genie
