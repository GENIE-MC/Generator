//____________________________________________________________________________
/*!

\class    genie::flux::GNuMIFlux

\brief    A GENIE flux driver encapsulating the NuMI neutrino flux.
          It reads-in the official GNUMI neutrino flux ntuples.

\ref      http://www.hep.utexas.edu/~zarko/wwwgnumi/v19/

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

          Robert Hatcher <rhatcher@fnal.gov>
          Fermi National Accelerator Laboratory

\created  Jun 27, 2008

\cpright  Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _GNUMI_NEUTRINO_FLUX_H_
#define _GNUMI_NEUTRINO_FLUX_H_

#include <string>
#include <iostream>

#include <TVector3.h>
#include <TLorentzVector.h>
#include <TLorentzRotation.h>

#include "EVGDrivers/GFluxI.h"
#include "PDG/PDGUtils.h"

class TFile;
class TChain;
class TBranch;

// MakeClass created classes for handling flux files
class g3numi;
class g4numi;

using std::string;
using std::ostream;

namespace genie {
namespace flux  {

class GNuMIFluxPassThroughInfo;

class GNuMIFlux: public GFluxI {

public :
  GNuMIFlux();
 ~GNuMIFlux();

  // Methods implementing the GENIE GFluxI interface, required for integrating
  // the NuMI neutrino flux simulations with the GENIE event generation drivers

  const PDGCodeList &    FluxParticles (void) { return *fPdgCList;            }
  double                 MaxEnergy     (void) { return  fMaxEv;               }
  bool                   GenerateNext  (void);
  int                    PdgCode       (void) { return  fgPdgC;               }
  double                 Weight        (void) { return  fWeight;              }
  const TLorentzVector & Momentum      (void) { return  fgP4;                 }
  const TLorentzVector & Position      (void) { return  fgX4;                 }
  bool                   End           (void) { return  fIEntry >= fNEntries
                                                     && fICycle == fNCycles;  }

  // Methods specific to the NuMI flux driver,
  // for configuration/initialization of the flux & event generation drivers and
  // and for passing-through flux information (eg neutrino parent decay kinematics)
  // not used by the generator but required by analyses/processing further upstream

  void LoadBeamSimData  (string filename);                     ///< load a gnumi root flux ntuple
  void SetFluxParticles (const PDGCodeList & particles);       ///< specify list of flux neutrino species
  void SetMaxEnergy     (double Ev);                           ///< specify maximum flx neutrino energy
  void SetFilePOT       (double pot);                          ///< POTs per input flux file
  void SetGenWeighted   (bool genwgt=false) { fGenWeighted = genwgt; } ///< toggle whether GenerateNext() returns weight=1 flux (initial default false)

  void SetNumOfCycles   (int n);                               ///< set how many times to cycle through the ntuple (default: 1 / n=0 means 'infinite')
  void SetTreeName      (string name);                         ///< set input tree name (default: "h10")
  void ScanForMaxWeight (void);                                ///< scan for max flux weight (before generating unweighted flux neutrinos)

  double   POT_curr       (void);                              ///< current average POT
  long int NFluxNeutrinos (void) const { return fNNeutrinos; } ///< number of flux neutrinos looped so far
  double   SumWeight      (void) const { return fSumWeight;  } ///< intergated weight for flux neutrinos looped so far

  const GNuMIFluxPassThroughInfo &
     PassThroughInfo(void) { return *fCurrentEntry; } ///< GNuMIFluxPassThroughInfo

  void PrintCurrent (void);  ///< print current entry from leaves

  // configure a flux window (or point) where E_nu and weight are evaluated

  void UseFluxAtNearDetCenter (void);  ///< MINOS NearDet "center" (as found in ntuple)
  void UseFluxAtFarDetCenter  (void);  ///< MINOS FarDet "center" (as found in ntuple)

  // GNuMIFlux uses "cm" as the length unit consistently internally (this the is
  // units used by both the g3 and g4 ntuples).  User interactions setting the
  // beam-to-detector coordinate transform, flux window, and the returned position
  // might need to be in other units.  
  typedef enum ELengthUnits { kmeter, kcm, kmm, kfm } LengthUnits_t;
  LengthUnits_t GetLengthUnits() const { return fLengthUnits; }
  void          SetLengthUnits(LengthUnits_t units) { fLengthUnits = units; SetLengthScale(); }

  typedef enum EStdFluxWindow {
    kMinosNearDet,      // appropriate for Near Detector
    kMinosFarDet,       // appropriate for Far Detector
    kMinosNearRock,     // appropriate for Near rock generation
    kMinosFarRock,      // appropriate for Far rock generation
    kMinosNearCenter,   // point location mimic near value in ntuple
    kMinosFarCenter     // point location mimic far value in ntuple
  } StdFluxWindow_t;

  // set both flux window and coordinate transform for some standard conditions
  bool SetBeamFluxWindow(StdFluxWindow_t stdwindow, double padding=0);  ///< return false if unhandled
  
  // rwh: upgrade allow flux window set/get in beam coords as optional flag to *etFluxWindow
  void SetFluxWindow(TLorentzVector  p1, TLorentzVector  p2, TLorentzVector  p3); ///< 3 points define a plane (by default in beam coordinates)
  void GetFluxWindow(TLorentzVector& p1, TLorentzVector& p2, TLorentzVector& p3) const; ///< 3 points define a plane in beam coordinate 
  void SetDetectorCoord(TLorentzVector det0, TLorentzRotation detrot);  ///< set detector position / orientation in beam system
  void GetDetectorCoord(TLorentzVector& det0, TLorentzRotation& detrot) const;  ///< set detector position / orientation in beam system

  void SetUpstreamZ     (double z0);                           ///< set flux neutrino initial z position (upstream of the detector)


private:

  // Private methods
  //
  bool GenerateNext_weighted (void);
  void Initialize            (void);
  void SetDefaults           (void);
  void CleanUp               (void);
  void ResetCurrent          (void);

  void SetLengthScale        (void);

  // Private data members
  //
  double         fMaxEv;          ///< maximum energy
  PDGCodeList *  fPdgCList;       ///< list of neutrino pdg-codes

  int            fgPdgC;          ///< running generated nu pdg-code
  TLorentzVector fgP4;            ///< running generated nu 4-momentum
  TLorentzVector fgX4;            ///< running generated nu 4-position

  string    fNuFluxFilePattern;   ///< wildcarded path
  string    fNuFluxTreeName;      ///< Tree name "h10" (g3) or "nudata" (g4)
  TChain*   fNuFluxTree;          ///< TTree in g3numi or g4numi // REF ONLY!
  g3numi*   fG3NuMI;              ///< g3numi ntuple
  g4numi*   fG4NuMI;              ///< g4numi ntuple
  int       fNFiles;              ///< number of files in chain

  long int  fNEntries;            ///< number of flux ntuple entries
  long int  fIEntry;              ///< current flux ntuple entry
  double    fMaxWeight;           ///< max flux neutrino weight in input file
  double    fFilePOT;             ///< file POT normalization
  double    fZ0;                  ///< configurable starting z position for each flux neutrino (in detector coord system)
  int       fNCycles;             ///< how many times to cycle through the flux ntuple
  int       fICycle;              ///< current cycle
  double    fWeight;              ///< current neutrino weight
  double    fSumWeight;           ///< sum of weights for neutrinos thrown so far
  long int  fNNeutrinos;          ///< number of flux neutrinos thrown so far
  bool      fGenWeighted;         ///< does GenerateNext() give weights?
  bool      fDetLocIsSet;         ///< is a flux location (near/far) set?
  // 
  LengthUnits_t fLengthUnits;       ///< units for coordinates in user exchanges
  double        fLengthScaleB2U;    ///< scale factor beam (cm) --> user
  TLorentzVector   fFluxWindowBase; ///< base point for flux window
  TLorentzVector   fFluxWindowDir1; ///< extent for flux window (direction 1)
  TLorentzVector   fFluxWindowDir2; ///< extent for flux window (direction 2)
  TLorentzVector   fDetectorZero;   ///< detector (0,0,0) in beam coordinates
  TLorentzRotation fDetectorRot;    ///< detector axis' unit vectors in beam coordinate

  int       fUseFluxAtDetCenter;  ///< use flux at near (-1) or far (+1) det center from ntuple?

  GNuMIFluxPassThroughInfo* fCurrentEntry;  ///< copy of current ntuple entry info (owned structure)

};

#define GNUMI_TEST_XY_WGT
#ifdef  GNUMI_TEST_XY_WGT
class xypartials { 
  // intermediate partial info from xy reweighting for comparison w/ f77 version
public:
  xypartials() { ; }
  void ReadStream(ifstream& myfile);
  int  Compare(const xypartials& other) const;
  void Print() const;
  // actual data
  double parent_mass, parentp, parent_energy;
  double gamma, beta_mag, enuzr, rad;
  double costh_pardet, theta_pardet, emrat, eneu;
  double sangdet, wgt;
  double betanu[3], p_nu[3], partial1, p_dcm_nu[4];
  double muparent_px, muparent_py, muparent_pz;
  double gammamp, betamp[3], partial2, p_pcm_mp[3], p_pcm;
  double costhmu, wgt_ratio;
  int ptype, ntype;
};
#endif

// A small persistable C-struct -like class that mirrors (some of) the structure of 
// the gnumi ntuples.  This can then be stored as an extra branch of the output event 
// tree -alongside with the generated event branch- for use further upstream in the 
// analysis chain - e.g. beam reweighting etc.
//
class GNuMIFluxPassThroughInfo: public TObject {
public:
   GNuMIFluxPassThroughInfo();
   /* allow default copy constructor ... for now nothing special
   GNuMIFluxPassThroughInfo(const GNuMIFluxPassThroughInfo & info);
   */
   virtual ~GNuMIFluxPassThroughInfo() { };

   void Copy(const g3numi*);
   void Copy(const g4numi*);

   void Reset();  //
   void ConvertPartCodes();

   int CalcEnuWgt(double xpos, double ypos, double zpos, double& enu, double& wgt_xy
#ifdef  GNUMI_TEST_XY_WGT
                  , xypartials& partials
#endif
                  ) const;

   friend ostream & operator << (ostream & stream, const GNuMIFluxPassThroughInfo & info);

   int   pcodes;  // 0=original GEANT particle codes, 1=converted to PDG
   int   units;   // 0=original GEANT cm, 1=meters

   // maintained variable names from gnumi ntuples
   // see http://www.hep.utexas.edu/~zarko/wwwgnumi/v19/[/v19/output_gnumi.html]

   Int_t    run;
   Int_t    evtno;
   Double_t ndxdz;
   Double_t ndydz;
   Double_t npz;
   Double_t nenergy;
   Double_t ndxdznea;
   Double_t ndydznea;
   Double_t nenergyn;
   Double_t nwtnear;
   Double_t ndxdzfar;
   Double_t ndydzfar;
   Double_t nenergyf;
   Double_t nwtfar;
   Int_t    norig;
   Int_t    ndecay;
   Int_t    ntype;
   Double_t vx;
   Double_t vy;
   Double_t vz;
   Double_t pdpx;
   Double_t pdpy;
   Double_t pdpz;
   Double_t ppdxdz;
   Double_t ppdydz;
   Double_t pppz;
   Double_t ppenergy;
   Int_t    ppmedium;
   Int_t    ptype;     // converted to PDG
   Double_t ppvx;
   Double_t ppvy;
   Double_t ppvz;
   Double_t muparpx;
   Double_t muparpy;
   Double_t muparpz;
   Double_t mupare;
   Double_t necm;
   Double_t nimpwt;
   Double_t xpoint;
   Double_t ypoint;
   Double_t zpoint;
   Double_t tvx;
   Double_t tvy;
   Double_t tvz;
   Double_t tpx;
   Double_t tpy;
   Double_t tpz;
   Int_t    tptype;   // converted to PDG
   Int_t    tgen;
   Int_t    tgptype;  // converted to PDG
   Double_t tgppx;
   Double_t tgppy;
   Double_t tgppz;
   Double_t tprivx;
   Double_t tprivy;
   Double_t tprivz;
   Double_t beamx;
   Double_t beamy;
   Double_t beamz;
   Double_t beampx;
   Double_t beampy;
   Double_t beampz;    

ClassDef(GNuMIFluxPassThroughInfo,3)
};

} // flux namespace
} // genie namespace

#endif // _GNUMI_NEUTRINO_FLUX_H_
