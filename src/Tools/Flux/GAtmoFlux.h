//____________________________________________________________________________
/*!

\class    genie::flux::GAtmoFlux

\brief    A base class for the FLUKA, BGLRS and ATMNC atmo. nu. flux drivers.
          The driver depends on data files provided by the atmospheric neutrino
          flux simulation authors in order to determine the angular and energy
          dependence for each neutrino species.
          The position of each flux neutrino [going towards a detector centered 
          at (0,0,0)] is generated uniformly on a plane that is perpendicular 
          to a sphere of radius Rl at the point that is determined by the 
          generated neutrino direction (theta,phi). The size of the area of 
          that plane, where flux neutrinos are generated, is determined by the 
          transverse radius Rt. You can tweak Rl, Rt to match the size of your 
          detector. 
          Initially, neutrino coordinates are generated in a default detector 
          coordinate system (Topocentric Horizontal Coordinate -THZ-):
             +z: Points towards the local zenith.
             +x: On same plane as local meridian, pointing south.
             +y: As needed to make a right-handed coordinate system.
             origin: detector centre
          Alternative user-defined topocentric systems can
          be defined by specifying the appropriate rotation from THZ.
          The driver allows minimum and maximum energy cuts.
          Also it provides the options to generate wither unweighted or weighted 
          flux neutrinos (the latter giving smoother distributions at the tails).

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  January 26, 2008

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _GATMO_FLUX_H_
#define _GATMO_FLUX_H_

#include <string>
#include <map>
#include <vector>

#include <TLorentzVector.h>
#include <TRotation.h>

#include "Framework/EventGen/GFluxI.h"

class TH3D;

using std::string;
using std::map;
using std::vector;

namespace genie {
namespace flux  {

class GAtmoFlux: public GFluxI {

public :
  virtual ~GAtmoFlux();

  // methods implementing the GENIE GFluxI interface
  virtual const PDGCodeList &    FluxParticles (void) { return *fPdgCList; }
  virtual double                 MaxEnergy     (void);
  virtual bool                   GenerateNext  (void);
  virtual int                    PdgCode       (void) { return fgPdgC;     }
  virtual double                 Weight        (void) { return fWeight;    }
  virtual const TLorentzVector & Momentum      (void) { return fgP4;       }
  virtual const TLorentzVector & Position      (void) { return fgX4;       }
  virtual bool                   End           (void) { return false;      }
  virtual long int               Index         (void) { return -1;         }
  virtual void                   Clear            (Option_t * opt);
  virtual void                   GenerateWeighted (bool gen_weighted);

  // get neutrino energy/direction of generated events
  double Enu        (void) { return fgP4.Energy(); }
  double Energy     (void) { return fgP4.Energy(); }
  double CosTheta   (void) { return -fgP4.Pz()/fgP4.Energy(); }
  double CosZenith  (void) { return -fgP4.Pz()/fgP4.Energy(); }

  // methods specific to the atmospheric flux drivers
  long int NFluxNeutrinos     (void) const; ///< Number of flux nu's generated. Not the same as the number of nu's thrown towards the geometry (if there are cuts).
  void     ForceMinEnergy     (double emin);
  void     ForceMaxEnergy     (double emax);
  void     SetSpectralIndex   (double index); 
  void     SetRadii           (double Rlongitudinal, double Rtransverse);
  void     SetUserCoordSystem (TRotation & rotation); ///< Rotation: Topocentric Horizontal -> User-defined Topocentric Coord System.
  void     AddFluxFile        (int neutrino_pdg, string filename);
  void     AddFluxFile        (string filename);
  bool     LoadFluxData       (void);

  TH3D*    GetFluxHistogram   (int flavour);
  double   GetFlux            (int flavour);
  double   GetFlux            (int flavour, double energy);
  double   GetFlux            (int flavour, double energy, double costh);
  double   GetFlux            (int flavour, double energy, double costh, double phi);

protected:

  // abstract class, ctor hidden
  GAtmoFlux();

  // protected methods
  bool    GenerateNext_1try (void);
  void    Initialize        (void);
  void    CleanUp           (void);
  void    ResetSelection    (void);
  double  MinEnergy         (void) { return fMinEvCut; }
  TH3D *  CreateFluxHisto   (string name, string title);
  void    ZeroFluxHisto     (TH3D * hist);
  void    AddAllFluxes      (void);
  int     SelectNeutrino    (double Ev, double costheta, double phi); 
  TH3D*   CreateNormalisedFluxHisto ( TH3D* hist);  // normalise flux files

  // pure virtual methods; to be implemented by concrete flux drivers
  virtual bool FillFluxHisto (int nu_pdg, string filename) = 0;

  // protected data members
  double           fMaxEv;              ///< maximum energy (in input flux files)
  PDGCodeList *    fPdgCList;           ///< input list of neutrino pdg-codes
  int              fgPdgC;              ///< current generated nu pdg-code
  TLorentzVector   fgP4;                ///< current generated nu 4-momentum
  TLorentzVector   fgX4;                ///< current generated nu 4-position
  double           fWeight;             ///< current generated nu weight
  long int         fNNeutrinos;         ///< number of flux neutrinos thrown so far
  double           fMaxEvCut;           ///< user-defined cut: maximum energy 
  double           fMinEvCut;           ///< user-defined cut: minimum energy  
  double           fRl;                 ///< defining flux neutrino generation surface: longitudinal radius
  double           fRt;                 ///< defining flux neutrino generation surface: transverse radius
  TRotation        fRotTHz2User;        ///< coord. system rotation: THZ -> Topocentric user-defined
  unsigned int     fNumPhiBins;         ///< number of phi bins in input flux data files
  unsigned int     fNumCosThetaBins;    ///< number of cos(theta) bins in input flux data files
  unsigned int     fNumEnergyBins;      ///< number of energy bins in input flux data files
  double *         fPhiBins;            ///< phi bins in input flux data files
  double *         fCosThetaBins;       ///< cos(theta) bins in input flux data files
  double *         fEnergyBins;         ///< energy bins in input flux data files
  bool             fGenWeighted;        ///< generate a weighted or unweighted flux?
  double           fSpectralIndex;      ///< power law function used for weighted flux
  bool             fInitialized;        ///< flag to check that initialization is run
  TH3D *           fTotalFluxHisto;     ///< flux = f(Ev,cos8,phi) summed over neutrino species
  double           fTotalFluxHistoIntg; ///< fFluxSum2D integral 
  map<int, TH3D*>  fFluxHistoMap;       ///< flux = f(Ev,cos8,phi) for each neutrino species
  map<int, TH3D*>  fRawFluxHistoMap;    ///< flux = f(Ev,cos8,phi) for each neutrino species
  vector<int>      fFluxFlavour;        ///< input flux file for each neutrino species
  vector<string>   fFluxFile;           ///< input flux file for each neutrino species
};

} // flux namespace
} // genie namespace

#endif // _GATMO_FLUX_H_

