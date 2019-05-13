//____________________________________________________________________________
/*!

\class    genie::flux::GNuMIFlux

\brief    A GENIE flux driver encapsulating the NuMI neutrino flux.
          It reads-in the official GNUMI neutrino flux ntuples.
          Supports both geant3 and geant4 formats.

\ref      http://www.hep.utexas.edu/~zarko/wwwgnumi/v19/

\author   Robert Hatcher <rhatcher \at fnal.gov>
          Fermi National Accelerator Laboratory

\created  Jun 27, 2008

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _GNUMI_NEUTRINO_FLUX_H_
#define _GNUMI_NEUTRINO_FLUX_H_

#include <string>
#include <iostream>
#include <vector>
#include <set>

#include <TVector3.h>
#include <TLorentzVector.h>
#include <TLorentzRotation.h>

#include "Framework/EventGen/GFluxI.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Tools/Flux/GFluxExposureI.h"
#include "Tools/Flux/GFluxFileConfigI.h"

class TFile;
class TChain;
class TTree;
class TBranch;

// MakeClass created classes for handling NuMI flux files
class g3numi;
class g4numi;
class flugg;

using std::string;
using std::ostream;

namespace genie {
namespace flux  {

class GNuMIFluxPassThroughInfo;
ostream & operator << (ostream & stream, const GNuMIFluxPassThroughInfo & info);

/// GNuMIFluxPassThroughInfo:
/// =========================
/// A small persistable C-struct -like class that mirrors (some of) the 
/// structure of the gnumi ntuples.  This can then be stored as an extra 
/// branch of the output event tree -alongside with the generated event 
/// branch- for use further upstream in the analysis chain - e.g. beam 
/// reweighting etc.
/// To do future x-y reweighting users must retain the info found in:
//     Ntype   Vx      Vy      Vz      
//     pdPx    pdPy    pdPz    
//     ppdxdz  ppdydz  pppz    ppenergy ptype
//     muparpx muparpy muparpz mupare   Necm
//     Nimpwt  
///
class GNuMIFluxPassThroughInfo: public TObject {
public:
   GNuMIFluxPassThroughInfo();
   /* allow default copy constructor ... for now nothing special
   GNuMIFluxPassThroughInfo(const GNuMIFluxPassThroughInfo & info);
   */
   virtual ~GNuMIFluxPassThroughInfo() { };

   void MakeCopy(const g3numi*);  ///< pull in from g3 ntuple
   void MakeCopy(const g4numi*);  ///< pull in from g4 ntuple
   void MakeCopy(const flugg*);   ///< pull in from flugg ntuple

   void ResetCopy();     // reset portion copied from ntuple
   void ResetCurrent();  // reset generated xy positioned info
   void ConvertPartCodes();
   void Print(const Option_t* opt = "") const;

   int CalcEnuWgt(const TLorentzVector& xyz, double& enu, double& wgt_xy) const;

   friend ostream & operator << (ostream & stream, const GNuMIFluxPassThroughInfo & info);

   int   pcodes;  // 0=original GEANT particle codes, 1=converted to PDG
   int   units;   // 0=original GEANT cm, 1=meters

   // Values for GNuMIFlux chosen x-y-z position, not from flux ntuple
   int            fgPdgC;   ///< generated nu pdg-code
   double         fgXYWgt;  ///< generated nu x-y weight
                            ///   not the same as GNuMIFlux::Weight()
                            ///   which include importance wgt and deweighting
   TLorentzVector fgP4;     ///< generated nu 4-momentum beam coord
   TLorentzVector fgX4;     ///< generated nu 4-position beam coord
   TLorentzVector fgP4User; ///< generated nu 4-momentum user coord
   TLorentzVector fgX4User; ///< generated nu 4-position user coord

   // Values copied from gnumi ntuples (generally maintained variable names)
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

#ifndef SKIP_MINERVA_MODS
   //=========================================
   // The following was inserted by MINERvA
   //=========================================
   int getProcessID(TString sval);
   int getVolID(TString sval);

   static const unsigned int MAX_N_TRAJ = 10; ///< Maximum number of trajectories to store

   Int_t    ntrajectory;
   Bool_t   overflow;
   int pdgcode[MAX_N_TRAJ];
   int trackId[MAX_N_TRAJ];
   int parentId[MAX_N_TRAJ];

   double startx[MAX_N_TRAJ];
   double starty[MAX_N_TRAJ];
   double startz[MAX_N_TRAJ];
   double startpx[MAX_N_TRAJ];
   double startpy[MAX_N_TRAJ];
   double startpz[MAX_N_TRAJ];
   double stopx[MAX_N_TRAJ];
   double stopy[MAX_N_TRAJ];
   double stopz[MAX_N_TRAJ];
   double stoppx[MAX_N_TRAJ];
   double stoppy[MAX_N_TRAJ];
   double stoppz[MAX_N_TRAJ];
   double pprodpx[MAX_N_TRAJ];
   double pprodpy[MAX_N_TRAJ];
   double pprodpz[MAX_N_TRAJ];

   int proc[MAX_N_TRAJ];
   int ivol[MAX_N_TRAJ];
   int fvol[MAX_N_TRAJ];
   //END of minerva additions
#endif

ClassDef(GNuMIFluxPassThroughInfo,5)
};

/// GNuMIFlux:
/// ==========
/// An implementation of the GFluxI interface that provides NuMI flux
///
class GNuMIFlux 
  : public genie::GFluxI
  , public genie::flux::GFluxExposureI
  , public genie::flux::GFluxFileConfigI 
{

public :
  GNuMIFlux();
 ~GNuMIFlux();

  // Methods implementing the GENIE GFluxI interface, required for integrating
  // the NuMI neutrino flux simulations with the GENIE event generation drivers

  const PDGCodeList &    FluxParticles (void) { return *fPdgCList;            }
  double                 MaxEnergy     (void) { return  fMaxEv;               }
  bool                   GenerateNext  (void);
  int                    PdgCode       (void) { return  fCurEntry->fgPdgC;    }
  double                 Weight        (void) { return  fWeight;              }
  const TLorentzVector & Momentum      (void) { return  fCurEntry->fgP4User;  }
  const TLorentzVector & Position      (void) { return  fCurEntry->fgX4User;  }
  bool                   End           (void) { return  fEnd;                 }
  long int               Index         (void) { return  fIEntry;              }
  void                   Clear            (Option_t * opt);
  void                   GenerateWeighted (bool gen_weighted);

  // Methods specific to the NuMI flux driver,
  // for configuration/initialization of the flux & event generation drivers 
  // and and for passing-through flux information (e.g. neutrino parent decay
  // kinematics) not used by the generator but required by analyses/processing 
  // further downstream

  //
  // information about or actions on current entry
  //
  const GNuMIFluxPassThroughInfo &
     PassThroughInfo(void) { return *fCurEntry; } ///< GNuMIFluxPassThroughInfo
  Long64_t GetEntryNumber() { return fIEntry; }   ///< index in chain

  double    GetDecayDist() const; ///< dist (user units) from dk to current pos
  void      MoveToZ0(double z0);  ///< move ray origin to user coord Z0

  //
  // information about the current state
  //
  virtual double    GetTotalExposure() const;  // GFluxExposureI interface
  virtual long int  NFluxNeutrinos() const;    ///< # of rays generated

  double    POT_curr(void);             ///< current average POT (RWH?)
  double    UsedPOTs(void) const;       ///< # of protons-on-target used

  double    SumWeight(void) const { return fSumWeight;  } ///< integrated weight for flux neutrinos looped so far

  void      PrintCurrent(void);         ///< print current entry from leaves
  void      PrintConfig();              ///< print the current configuration

  std::vector<std::string> GetFileList();  ///< list of files currently part of chain

  // 
  // GFluxFileConfigI interface
  //
  virtual void  LoadBeamSimData(const std::vector<std::string>& filenames,
                                const std::string&              det_loc);
  using GFluxFileConfigI::LoadBeamSimData; // inherit the rest
  virtual void GetBranchInfo(std::vector<std::string>& branchNames,
                             std::vector<std::string>& branchClassNames,
                             std::vector<void**>&      branchObjPointers);
  virtual TTree* GetMetaDataTree();

  //
  // configuration of GNuMIFlux
  //

  bool      LoadConfig(string cfg);                               ///< load a named configuration
  void      SetMaxEnergy(double Ev);                              ///< specify maximum flx neutrino energy

  void      SetGenWeighted(bool genwgt=false) { fGenWeighted = genwgt; } ///< toggle whether GenerateNext() returns weight=1 flux (initial default false)

  void      SetEntryReuse(long int nuse=1);                       ///<  # of times to use entry before moving to next

  void      SetTreeName(string name);                             ///< set input tree name (default: "h10")
  void      ScanForMaxWeight(void);                               ///< scan for max flux weight (before generating unweighted flux neutrinos)
  void      SetMaxWgtScan(double fudge = 1.05, long int nentries = 2500000)      ///< configuration when estimating max weight
            { fMaxWgtFudge = fudge; fMaxWgtEntries = nentries; }
  void      SetMaxEFudge(double fudge = 1.05)                     ///< extra fudge factor in estimating maximum energy
            { fMaxEFudge = fudge; }
  void      SetApplyWindowTiltWeight(bool apply = true)           ///< apply wgt due to tilt of flux window relative to beam
            { fApplyTiltWeight = apply; }


  // GNuMIFlux uses "cm" as the length unit consistently internally (this is 
  // the length units used by both the g3 and g4 ntuples).  User interactions 
  // setting the beam-to-detector coordinate transform, flux window, and the 
  // returned position might need to be in other units.  Use:
  //     double scale = genie::utils::units::UnitFromString("cm");
  // ( #include "Utils/UnitUtils.h for declaration )
  // to get the correct scale factor to pass in.  This should get set
  // FIRST before setting detector position/rotation

  void   SetLengthUnits(double user_units);  ///< Set units assumed by user
  double    LengthUnits(void) const;         ///< Return user units
  
  // set relative orientation of user coords vs. beam system, i.e.
  //  x3_user = ( beamrot * x3_beam ) + x0beam_user
  //  p3_user =   beamrot * p3_beam 

  ///< beam (0,0,0) relative to user frame, beam direction in user frame
  void      SetBeamRotation(TRotation beamrot);
  void      SetBeamCenter(TVector3 beam0);
  TRotation GetBeamRotation() const; ///< rotation to apply from beam->user
  TVector3  GetBeamCenter() const;   ///< beam origin in user frame

  // configure a flux window (or point) where E_nu and weight are evaluated

  typedef enum EStdFluxWindow {
    kMinosNearDet,      // appropriate for Near Detector
    kMinosFarDet,       // appropriate for Far Detector
    kMinosNearRock,     // appropriate for Near rock generation
    kMinosFarRock,      // appropriate for Far rock generation
    kMinosNearCenter,   // point location mimic near value in ntuple
    kMinosFarCenter     // point location mimic far value in ntuple
  } StdFluxWindow_t;

  // set both flux window in user coord and coordinate transform 
  // for some standard conditions
  bool      SetFluxWindow(StdFluxWindow_t stdwindow, double padding=0);  ///< return false if unhandled
  
  // rwh: potential upgrade: allow flux window set/get in beam coords 
  // as optional flag to *etFluxWindow
  void      SetFluxWindow(TVector3  p1, TVector3  p2, TVector3  p3); ///< 3 points define a plane (by default in user coordinates)
  void      GetFluxWindow(TVector3& p1, TVector3& p2, TVector3& p3) const; ///< 3 points define a plane in beam coordinate 

  /// force weights at MINOS detector "center" as found in ntuple
  void      UseFluxAtNearDetCenter(void);
  void      UseFluxAtFarDetCenter(void);

  //
  // Actual coordinate transformations  b=beam, u=user (e.g. detector)
  //
  void      Beam2UserPos(const TLorentzVector& beamxyz,
                               TLorentzVector& usrxyz  ) const;
  void      Beam2UserDir(const TLorentzVector& beamdir,
                               TLorentzVector& usrdir  ) const;
  void      Beam2UserP4 (const TLorentzVector& beamp4,
                               TLorentzVector& usrp4   ) const;
  void      User2BeamPos(const TLorentzVector& usrxyz,
                               TLorentzVector& beamxyz ) const;
  void      User2BeamDir(const TLorentzVector& usrdir,
                               TLorentzVector& beamdir ) const;
  void      User2BeamP4 (const TLorentzVector& usrp4,
                               TLorentzVector& beamp4  ) const;

  TVector3  FluxWindowNormal() { return fWindowNormal; }

private:

  // Private methods
  //
  bool GenerateNext_weighted (void);
  void Initialize            (void);
  void SetDefaults           (void);
  void CleanUp               (void);
  void ResetCurrent          (void);
  void AddFile               (TTree* tree, string fname);
  void CalcEffPOTsPerNu      (void);
  
  // Private data members
  //
  double         fMaxEv;          ///< maximum energy
  bool           fEnd;            ///< end condition reached

  std::vector<string> fNuFluxFilePatterns;   ///< (potentially wildcarded) path(s)
  string    fNuFluxTreeName;      ///< Tree name "h10" (g3) or "nudata" (g4)
  TChain*   fNuFluxTree;          ///< TTree in g3numi or g4numi // REF ONLY!
  string    fNuFluxGen;           ///< "g3numi" "g4numi" or "flugg"
  g3numi*   fG3NuMI;              ///< g3numi ntuple
  g4numi*   fG4NuMI;              ///< g4numi ntuple
  flugg*    fFlugg;               ///< flugg ntuple
  int       fNFiles;              ///< number of files in chain
  Long64_t  fNEntries;            ///< number of flux ntuple entries
  Long64_t  fIEntry;              ///< current flux ntuple entry
  Long64_t  fNuTot;               ///< cummulative # of entries (=fNEntries)
  Long64_t  fFilePOTs;            ///< # of protons-on-target represented by all files

  double    fWeight;              ///< current neutrino weight, =1 if generating unweighted entries
  double    fMaxWeight;           ///< max flux neutrino weight in input file
  double    fMaxWgtFudge;         ///< fudge factor for estimating max wgt
  long int  fMaxWgtEntries;       ///< # of entries in estimating max wgt
  double    fMaxEFudge;           ///< fudge factor for estmating max enu (0=> use fixed 120GeV)

  long int  fNUse;                ///< how often to use same entry in a row
  long int  fIUse;                ///< current # of times an entry has been used
  double    fSumWeight;           ///< sum of weights for nus thrown so far
  long int  fNNeutrinos;          ///< number of flux neutrinos thrown so far
  double    fEffPOTsPerNu;        ///< what a entry is worth ...
  double    fAccumPOTs;           ///< POTs used so far

  bool      fGenWeighted;         ///< does GenerateNext() give weights?
  bool      fApplyTiltWeight;     ///< wgt due to window normal not || beam 
  bool      fDetLocIsSet;         ///< is a flux location (near/far) set?
  int       fUseFluxAtDetCenter;  ///< use flux at near (-1) or far (+1) det center from ntuple?
  
  double           fLengthUnits;    ///< units for coord in user exchanges
  double           fLengthScaleB2U; ///< scale factor beam (cm) --> user
  double           fLengthScaleU2B; ///< scale factor beam user --> (cm)

  TLorentzVector   fBeamZero;       ///< beam origin in user coords
  TLorentzRotation fBeamRot;        ///< rotation applied beam --> user coord
  TLorentzRotation fBeamRotInv;

  TVector3         fFluxWindowPtUser[3]; ///<  user points of flux window
  TLorentzVector   fFluxWindowBase; ///< base point for flux window - beam coord
  TLorentzVector   fFluxWindowDir1; ///< extent for flux window (direction 1)
  TLorentzVector   fFluxWindowDir2; ///< extent for flux window (direction 2)
  double           fFluxWindowLen1;
  double           fFluxWindowLen2;
  TVector3         fWindowNormal;   ///< normal direction for flux window

  TLorentzVector   fgX4dkvtx;       ///< decay 4-position beam coord

  GNuMIFluxPassThroughInfo* fCurEntry;  ///< copy of current ntuple entry info (owned structure)

};

//#define GNUMI_TEST_XY_WGT
#ifdef  GNUMI_TEST_XY_WGT
class xypartials { 
  // intermediate partial info from xy reweighting for comparison w/ f77 version
   friend ostream & operator << (ostream & stream, const xypartials & info);
public:
  xypartials() { ; }
  void ReadStream(std::ifstream& myfile);
  int  Compare(const xypartials& other) const;
  void Print(const Option_t* opt = "") const;
  static xypartials& GetStaticInstance(); // copy used by CalcEnuWgt()
  // actual data
  double xdet, ydet, zdet;
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

} // flux namespace
} // genie namespace

#endif // _GNUMI_NEUTRINO_FLUX_H_
