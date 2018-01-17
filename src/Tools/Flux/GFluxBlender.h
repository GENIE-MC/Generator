////////////////////////////////////////////////////////////////////////
/// \file  GFluxBlender.h
/// \brief GENIE GFluxI adapter to allow flavor modification
///
///        This adapter intervenes between the GENIE GMCJDriver class
///        (MC job driver) and a concrete GFluxI flux generator, to
///        allow user modification of the neutrino flavor.  This
///        modification could be a fixed "swap" or an energy and/or
///        distance dependent (standard oscillations) one.
///
///        Because the GMCJDriver only queries the flavor of a 
///        generated neutrino once, prior to propagation through
///        the geometry, this approach is _not_ appropriate with
///        use of an oscillatory model in situations where the flavor 
///        might change significantly over the scale of the geometry.
///        In such cases one would have to generate with a fixed flavor
///        (energy/distance independent) swap and reweight after the fact.
///
///        Do not use this as a means of selecting only certain flavor
///        from flux generators that support other means (e.g. GNuMIFlux,
///        GSimpleNtpFlux which have SetFluxParticles(PDGCodeList)) as
///        those will be more efficient.
///
/// \version $Id: GFluxBlender.h,v 1.1.1.1 2010/12/22 16:18:52 p-nusoftart Exp $
/// \author  Robert Hatcher <rhatcher \at fnal.gov>
///          Fermi National Accelerator Laboratory
///
////////////////////////////////////////////////////////////////////////
#ifndef GENIE_FLUX_GFLUXBLENDER_H
#define GENIE_FLUX_GFLUXBLENDER_H

#include <vector>
#include "EVGDrivers/GFluxI.h"
#include "PDG/PDGCodeList.h"

namespace genie {
namespace flux {

  // forward declarations within namespace
  class GFlavorMixerI;
  class GNuMIFlux;
  class GSimpleNtpFlux;

  class GFluxBlender : public GFluxI {
    
  public:
  
    GFluxBlender();
    ~GFluxBlender();

    //
    // implement the GFluxI interface:
    //   overriding real GFluxI methods:
    //      FluxParticles()  [final list of flavor might be more than source]
    //      GenerateNext()   [choose flavor during generation, avoid sterile]
    //
    const PDGCodeList &    FluxParticles (void); ///< declare list of flux neutrinos that can be generated (for init. purposes)
    double                 MaxEnergy     (void) { return fRealGFluxI->MaxEnergy(); } ///< declare the max flux neutrino energy that can be generated (for init. purposes)
    bool                   GenerateNext  (void); ///< generate the next flux neutrino (return false in err)
    int                    PdgCode       (void) { return fPdgCMixed; } ///< returns the flux neutrino pdg code
    double                 Weight        (void) { return fRealGFluxI->Weight(); } ///< returns the flux neutrino weight (if any)
    const TLorentzVector & Momentum      (void) { return fRealGFluxI->Momentum(); } ///< returns the flux neutrino 4-momentum 
    const TLorentzVector & Position      (void) { return fRealGFluxI->Position(); } ///< returns the flux neutrino 4-position (note: expect SI rather than physical units)
    bool                   End           (void) { return fRealGFluxI->End(); }  ///< true if no more flux nu's can be thrown (eg reaching end of beam sim ntuples)
    long int               Index            (void);
    void                   Clear            (Option_t * opt);
    void                   GenerateWeighted (bool gen_weighted);

    //
    // Additions to the GFluxI interface:
    //
    int             PdgCodeGenerated (void) { return fPdgCGenerated; } ///< returns the flux neutrino original pdg code
    double          Energy           (void) { return fEnergy; } //< returns the current neutrino's energy
    double          TravelDist       (void) { return fDistance; } ///< returns the distance used in the flavor mixing
    //
    // For flux methods that report a distance from the point of origin
    // to the chosen starting point for the ray use that distance
    // (supported in GNuMIFlux and GSimpleNtpFlux).
    // For other cases use a fixed travel distance (meters) set via:
    //
    void            SetBaselineDist  (double dist) { fBaselineDist = dist; }
    double          GetBaselineDist  (void) { return fBaselineDist; }

    //
    // Configuration:
    //
    GFluxI*         AdoptFluxGenerator(GFluxI* generator);    ///< return previous
    GFlavorMixerI*  AdoptFlavorMixer(GFlavorMixerI* mixer);   ///< return previous
    GFluxI*         GetFluxGenerator() { return fRealGFluxI; }  ///< access, not ownership
    GFlavorMixerI*  GetFlavorMixer()   { return fFlavorMixer; } ///< access, not ownership

    void            PrintConfig(void);
    void            PrintState(bool verbose=true);

  private:
    int             ChooseFlavor(int pdg_init, double energy, double dist);

    GFluxI*         fRealGFluxI;        ///< actual flux generator
    GNuMIFlux*      fGNuMIFlux;         ///< ref to avoid repeat dynamic_cast
    GSimpleNtpFlux* fGSimpleFlux;       ///< ref to avoid repeat dynamic_cast

    GFlavorMixerI*  fFlavorMixer;       ///< flavor modification schema

    PDGCodeList     fPDGListGenerator;  ///< possible flavors from generator
    PDGCodeList     fPDGListMixed;      ///< possible flavors after mixing
    size_t          fNPDGOut;           ///< # of possible output flavors

    double          fBaselineDist;      ///< travel dist for mixing (if flux doesn't support GetDecayDist())

    double          fEnergy;            ///< current neutrino's energy
    double          fDistance;          ///< current neutrino's travel distance
    int             fPdgCGenerated;     ///< current neutrino's original flavor
    int             fPdgCMixed;         ///< current neutrino's new flavor

    std::vector<double> fProb;          ///< individual transition probs
    std::vector<double> fSumProb;       ///< cummulative probability
    double              fRndm;          ///< random # used to make choice

  };

} // namespace flux
} // namespace genie
#endif //GENIE_FLUX_GFLUXBLENDER_H
