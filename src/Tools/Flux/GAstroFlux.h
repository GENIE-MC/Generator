//____________________________________________________________________________
/*!

\class    genie::flux::GAstroFlux

\brief    A base class for the concrete astrophysical neutrino flux drivers.


          <<<< NOTE: CODE UNDER DEVELOPMENT >>>>


          COORDINATE SYSTEMS / NEUTRINO GENERATION & PROPAGATION :

          Neutrinos are generated on a sphere with radius R_{earth}. 
          Especially,
          - For diffuse fluxes: 
            Neutrinos can be generated anywhere on that surface. 
          - For point sources: 
            Neutrinos are generated at fixed right ascension and declination. 
            Then time is randomized, to account for Earth's rotation, and the 
            Equatorial Coordinates are converted to GEF. So, neutrinos from each
            point source are generated on circles parallel to the Earth's Equator.

          Initially, neutrino coordinates are generated in the Geocentric Earth-
          Fixed (GEF) Coordinate System (later to be converted to the appropriate
          detector coordinate system - See further below).

          * Definition: 
            Geocentric Earth-Fixed (GEF) Coordinate System
              +z: Points to North Pole.
              xy: Equatorial plane.
              +x: Points to the Prime Meridian.
              +y: As needed to make a right-handed coordinate system.

          Neutrinos are then propagated towards the detector.
          The Earth opaqueness to ultra high energy neutrinos is taken into
          account. The Earth density profile is modelled using the PREM 
          (Preliminary Earth Model, The Encyclopedia of Solid Earth Geophysics,
          David E. James, ed., Van Nostrand Reinhold, New York, 1989, p.331).

          The detector position is determined in the Spherical/Geographic System 
          by its geographic latitude (angle relative to Equator), its geographic 
          longitude (angle relative to Prime Meridian) and its depth from the 
          surface.

          The generated flux neutrinos, propagated through the Earth towards the
          detector) are then positioned on the surface of a sphere with radius Rd 
          which should fully enclose the neutrino detector. The centre of that
          sphere is taken to be the origin of the detector coordinate system.
          The transverse coords are appropriately randomized so that neutrinos 
          from any given direction bath the entire sphere enclosing the detector. 

          The final flux neutrino coordinates are given in the detector coordinate 
          system. The default detector coordinate system is the Topocentric Horizontal
          (THZ) Coordinate System. Alternative user-defined topocentric systems can
          be defined by specifying the appropriate rotation from THZ.

          * Definition: 
            Topocentric Horizontal (THZ) Coordinate System 
            (default detector coordinate system)
                +z: Points towards the local zenith.
                +x: On same plane as local meridian, pointing south.
                +y: As needed to make a right-handed coordinate system.
		origin: detector centre

          WEIGHTING SCHEMES:

          For a detector with geometrical cross section ~ 1km^2, the solid
          angle acceptance changes by ~10 orders of magnitude across the 
          surface of the Earth.

          The driver supports both weighted and un-weighted flux generation
          schemes. However, because of the enormous changes in solid angle
          acceptance and energy, only the weighted scheme is practical.

          PHYSICS:

          The relative neutrino population needs to be set by the user using
          the SetRelNuPopulations() method.
          If run without arguments, the following relative populations are set:
          nue:numu:nutau:nuebar:numubar:nutaubar = 1:2:0:1:2:0

          The energy spectrum is follows a power law. The user needs to 
          specify the power-law index by calling SetEnergyPowLawIdx().

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  March 27, 2010

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _GASTRO_FLUX_H_
#define _GASTRO_FLUX_H_

#include <string>
#include <map>

#include <TLorentzVector.h>
#include <TVector3.h>
#include <TRotation.h>

#include "Framework/Conventions/Units.h"
#include "Framework/EventGen/GFluxI.h"

class TH1D;
class TH2D;

using std::string;
using std::map;

namespace genie {
namespace flux  {

const double kAstroDefMaxEv      = 1E+20 * units::GeV; ///< 
const double kAstroDefMinEv      = 1E-3  * units::GeV; ///< 
const int    kAstroNlog10EvBins  = 1000;               ///<
const int    kAstroNCosThetaBins = 500;                ///<
const int    kAstroNPhiBins      = 500;                ///<

//
//
//

class GAstroFlux: public GFluxI {

public :
  virtual ~GAstroFlux();

  //
  // methods implementing the GENIE GFluxI interface
  //
  virtual const PDGCodeList &    FluxParticles (void) { return *fPdgCList; }
  virtual double                 MaxEnergy     (void);
  virtual bool                   GenerateNext  (void);
  virtual int                    PdgCode       (void) { return fgPdgC;     }
  virtual double                 Weight        (void) { return fgWeight;   }
  virtual const TLorentzVector & Momentum      (void) { return fgP4;       }
  virtual const TLorentzVector & Position      (void) { return fgX4;       }
  virtual bool                   End           (void) { return false;      }
  virtual long int               Index         (void) { return -1;         }
  virtual void                   Clear            (Option_t * opt);
  virtual void                   GenerateWeighted (bool gen_weighted);

  //
  // configuration methods specific to all astrophysical neutrino flux drivers
  //
  void ForceMinEnergy      (double emin);
  void ForceMaxEnergy      (double emax);
  void SetDetectorPosition (double latitude, double longitude, double depth, double size);
  void SetRelNuPopulations (double nnue=1, double nnumu=2, double nnutau=0, double nnuebar=1, double nnumubar=2, double nnutaubar=0);
  void SetEnergyPowLawIdx  (double n);
  void SetUserCoordSystem  (TRotation & rotation); ///< rotation Topocentric Horizontal -> User-defined Topocentric Coord System

protected:

  class NuGenerator;
  class NuPropagator;

  // abstract class, ctor hidden
  GAstroFlux();

  // protected methods
  void Initialize               (void);
  void CleanUp                  (void);
  void ResetSelection           (void);

  //
  // protected data members
  //

  PDGCodeList *    fPdgCList;             ///< declared list of neutrino pdg-codes that can be thrown by current instance
  int              fgPdgC;                ///< (current) generated nu pdg-code
  TLorentzVector   fgP4;                  ///< (current) generated nu 4-momentum
  TLorentzVector   fgX4;                  ///< (current) generated nu 4-position
  double           fgWeight;              ///< (current) generated nu weight
  // configuration properties set by the user
  double           fMaxEvCut;             ///< (config) user-defined maximum energy cut
  double           fMinEvCut;             ///< (config) user-defined minimum energy cut
  bool             fGenWeighted;          ///< (config) generate a weighted or unweighted flux?
  double           fDetGeoLatitude;       ///< (config) detector: geographic latitude
  double           fDetGeoLongitude;      ///< (config) detector: geographic longitude
  double           fDetGeoDepth;          ///< (config) detector: depth from surface
  double           fDetSize;              ///< (config) detector: size (detector should be enclosed in sphere of this radius)
  map<int,double>  fRelNuPopulations;     ///< (config) relative neutrino populations
  TRotation        fRotGEF2THz;         ///< (config) coord. system rotation: GEF translated to detector centre -> THZ
  TRotation        fRotTHz2User;         ///< (config) coord. system rotation: THZ -> Topocentric user-defined 
  // internal flags and utility objects
  TVector3         fDetCenter;            ///<
  TH1D *           fEnergySpectrum;       ///<
  TH2D *           fSolidAngleAcceptance; ///<
  NuGenerator *    fNuGen;                ///<
  NuPropagator *   fNuPropg;              ///<

  //
  // utility classes
  //
  class NuGenerator {
  public:
    NuGenerator() {}
   ~NuGenerator() {}
    bool SelectNuPdg (bool weighted, const map<int,double> & nupdgpdf, int & nupdg, double & wght);
    bool SelectEnergy(bool weighted, TH1D & log10epdf, double log10emin, double log10emax, double & log10e, double & wght);
    bool SelectOrigin(bool weighted, TH2D & opdf, double & phi, double & costheta, double & wght);
  };
  class NuPropagator {
  public:
    NuPropagator(double stepsz) : fStepSize(stepsz/units::km) { }
   ~NuPropagator() { }
    bool Go(double phi_start, double costheta_start, const TVector3 & detector_centre, double detector_sz, int nu_pdg, double Ev);
    int        NuPdgAtDetVolBoundary (void) { return fNuPdg; }
    TVector3 & X3AtDetVolBoundary    (void) { return fX3;    }
    TVector3 & P3AtDetVolBoundary    (void) { return fP3;    }
  private:
    double   fStepSize;
    int      fNuPdg;
    TVector3 fX3;
    TVector3 fP3;
  };

};


//
// Concrete astrophysical flux drivers
//

//............................................................................
// GENIE diffuse astrophysical neutrino flux driver 
//
class GDiffuseAstroFlux: public GAstroFlux {
public :
  GDiffuseAstroFlux();
 ~GDiffuseAstroFlux();

  //
  //
}; 

//............................................................................
// GENIE concrete flux driver for astrophysical point neutrino sources
//
class GPointSourceAstroFlux: public GAstroFlux {
public :
  GPointSourceAstroFlux();
 ~GPointSourceAstroFlux();

  //
  bool GenerateNext (void);

  //
  void AddPointSource(string name, double ra, double dec, double rel_intensity);

private:

  bool SelectSource(void);

  map<int, string> fPntSrcName;  ///< point source name
  map<int, double> fPntSrcRA;    ///< right ascension
  map<int, double> fPntSrcDec;   ///< declination
  map<int, double> fPntSrcRelI;  ///< relative intensity
  double           fPntSrcTotI;  ///< sum of all relative intensities

  unsigned int fSelSourceId;
}; 
//............................................................................

} // flux namespace
} // genie namespace

#endif // _GASTRO_FLUX_H_

