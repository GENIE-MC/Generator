//____________________________________________________________________________
/*!

\class    genie::flux::GAtmoFlux

\brief    A base class for the concrete FLUKA and BGLRS atmospheric neutrino 
          flux drivers.
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
          STFC, Rutherford Appleton Laboratory

\created  January 26, 2008

\cpright  Copyright (c) 2003-2010, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _GATMO_FLUX_H_
#define _GATMO_FLUX_H_

#include <string>
#include <map>

#include <TLorentzVector.h>
#include <TRotation.h>

#include "EVGDrivers/GFluxI.h"

class TH2D;

using std::string;
using std::map;

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

  // methods specific to the atmospheric flux drivers
  void ForceMinEnergy     (double emin);
  void ForceMaxEnergy     (double emax);
  void GenerateWeighted   (bool gen_weighted);
  void SetRadii           (double Rlongitudinal, double Rtransverse);
  void SetUserCoordSystem (TRotation & rotation); ///< rotation: Topocentric Horizontal -> User-defined Topocentric Coord System
  void SetFluxFile        (int neutrino_pdg, string filename);
  bool LoadFluxData       (void);

protected:

  // abstract class, ctor hidden
  GAtmoFlux();

  // protected methods
  bool    GenerateNext_1try (void);
  void    Initialize        (void);
  void    CleanUp           (void);
  void    ResetSelection    (void);
  double  MinEnergy         (void) { return fMinEvCut; }
  TH2D *  CreateFluxHisto2D (string name, string title);
  void    ZeroFluxHisto2D   (TH2D * h2);
  void    AddAllFluxes      (void);
  int     SelectNeutrino    (double Ev, double costheta);
  
  // pure virtual protected methods; to be implemented by concrete flux drivers
  virtual bool FillFluxHisto2D   (TH2D * h2, string filename) = 0;

  // protected data members
  double           fMaxEv;            ///< (input) maximum energy (in input flux files)
  PDGCodeList *    fPdgCList;         ///< (input) list of neutrino pdg-codes
  int              fgPdgC;            ///< (current) generated nu pdg-code
  TLorentzVector   fgP4;              ///< (current) generated nu 4-momentum
  TLorentzVector   fgX4;              ///< (current) generated nu 4-position
  double           fWeight;           ///< (current) generated nu weight
  double           fMaxEvCut;         ///< (config) user-defined maximum energy cut
  double           fMinEvCut;         ///< (config) user-defined minimum energy cut
  map<int, string> fFluxFile;         ///< (config) input flux file for each neutrino species
  double           fRl;               ///< (config) flux neutrino generation surface: longitudinal radius
  double           fRt;               ///< (config) flux neutrino generation surface: transverse radius
  TRotation        fRotTHz2User;      ///< (config) coord. system rotation: THZ -> Topocentric user-defined
  unsigned int     fNumCosThetaBins;  ///< (config) number of cos(theta) bins in input flux data files
  unsigned int     fNumEnergyBins;    ///< (config) number of energy bins in input flux data files
  double *         fCosThetaBins;     ///< (config) cos(theta) bins in input flux data files
  double *         fEnergyBins;       ///< (config) energy bins in input flux data files
  bool             fGenWeighted;      ///< (config) generate a weighted or unweighted flux?
  map<int, TH2D*>  fFlux2D;           ///< flux = f(Ev,cos8) for each neutrino species
  TH2D *           fFluxSum2D;        ///< flux = f(Ev,cos8) summed over neutrino species
  double           fFluxSum2DIntg;    ///< fFluxSum2D integral 
};

} // flux namespace
} // genie namespace

#endif // _GATMO_FLUX_H_

