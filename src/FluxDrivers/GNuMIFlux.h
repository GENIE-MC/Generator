//____________________________________________________________________________
/*!

\class    genie::flux::GNuMIFlux

\brief    A GENIE flux driver encapsulating the NuMI neutrino flux.
          It reads-in the official GNUMI neutrino flux ntuples.

\ref      http://www.hep.utexas.edu/~zarko/wwwgnumi/v19/

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          STFC, Rutherford Appleton Laboratory

          Robert Hatcher <rhatcher@fnal.gov>
          Fermi National Accelerator Laboratory

\created  Jun 27, 2008

\cpright  Copyright (c) 2003-2008, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _GNUMI_NEUTRINO_FLUX_H_
#define _GNUMI_NEUTRINO_FLUX_H_

#include <string>
#include <iostream>

#include <TLorentzVector.h>
#include <TVector3.h>

#include "EVGDrivers/GFluxI.h"
#include "PDG/PDGUtils.h"

class TFile;
class TTree;
class TBranch;

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
  void SetUpstreamZ     (double z0);                           ///< set flux neutrino initial z position (upstream of the detector)
  void SetNumOfCycles   (int n);                               ///< set how many times to cycle through the ntuple (default: 1 / n=0 means 'infinite')
  void SetTreeName      (string name);                         ///< set input tree name (default: "h10")
  void ScanForMaxWeight (void);                                ///< scan for max flux weight (before generating unweighted flux neutrinos)
  void UseFluxAtFarDetCenter  (void);
  void UseFluxAtNearDetCenter (void);

  void PrintCurrent (void);  ///< print current entry from leaves

  double   POT_curr       (void);                              ///< current average POT
  long int NFluxNeutrinos (void) const { return fNNeutrinos; } ///< number of flux neutrinos looped so far
  double   SumWeight      (void) const { return fSumWeight;  } ///< intergated weight for flux neutrinos looped so far

  const GNuMIFluxPassThroughInfo &
     PassThroughInfo(void) { return *fCurrentEntry; } ///< GNuMIFluxPassThroughInfo

private:

  // Private methods
  //
  bool GenerateNext_weighted (void);
  void Initialize            (void);
  void SetDefaults           (void);
  void CleanUp               (void);
  void ResetCurrent          (void);

  // Private data members
  //
  double         fMaxEv;       ///< maximum energy
  PDGCodeList *  fPdgCList;    ///< list of neutrino pdg-codes

  int            fgPdgC;       ///< running generated nu pdg-code
  TLorentzVector fgP4;         ///< running generated nu 4-momentum
  TLorentzVector fgX4;         ///< running generated nu 4-position

  TFile *   fNuFluxFile;       ///< input flux file
  TTree *   fNuFluxTree;       ///< input flux ntuple
  string    fNuFluxTreeName;   ///< input flux ntuple name
  long int  fNEntries;         ///< number offlux ntuple entries
  long int  fIEntry;           ///< current flux ntuple entry
  double    fMaxWeight;        ///< max flux neutrino weight in input file
  double    fFilePOT;          ///< file POT normalization
  double    fZ0;               ///< configurable starting z position for each flux neutrino (in detector coord system)
  int       fNCycles;          ///< how many times to cycle through the flux ntuple
  int       fICycle;           ///< current cycle
  double    fWeight;           ///< current neutrino weight
  double    fSumWeight;        ///< sum of weights for neutrinos thrown so far
  long int  fNNeutrinos;       ///< number of flux neutrinos thrown so far
  int       fUseFluxAtDetCenter;  ///< use flux at near (-1) or far (+1) det center?
  bool      fDetLocIsSet;         ///< is a flux location (near/far) set?

  //-- gnumi ntuple branches
  //   maintained names/comments from gnumi documentation
  //   http://www.hep.utexas.edu/~zarko/wwwgnumi/v19/
  //
  TBranch * fBr_run;           ///< Run number
  TBranch * fBr_evtno;         ///< Event number (proton on target)
  TBranch * fBr_Ndxdz;         ///< Neutrino direction slopes for a random decay
  TBranch * fBr_Ndydz;         ///< see above
  TBranch * fBr_Npz;           ///< Neutrino momentum (GeV/c) along z direction (beam axis)
  TBranch * fBr_Nenergy;       ///< Neutrino energy (GeV/c) for a random decay
  TBranch * fBr_NdxdzNea;      ///< Neutrino direction slopes for a decay forced at center of near detector
  TBranch * fBr_NdydzNea;      ///< see above
  TBranch * fBr_NenergyN;      ///< Neutrino energy for a decay forced at center of near detector
  TBranch * fBr_NWtNear;       ///< Neutrino weight for a decay forced at center of near detector
  TBranch * fBr_NdxdzFar;      ///< Repeats of above, for a decay forces to the center of the far detector.
  TBranch * fBr_NdydzFar;      ///< see above
  TBranch * fBr_NenergyF;      ///< see above
  TBranch * fBr_NWtFar;        ///< see above
  TBranch * fBr_Norig;         ///< Obsolete...
  TBranch * fBr_Ndecay;        ///< Decay mode that produced neutrino
                               ///   1  K0L -> nue pi- e+
                               ///   2  K0L -> nuebar pi+ e-
                               ///   3  K0L -> numu pi- mu+
                               ///   4  K0L -> numubar pi+ mu-
                               ///   5  K+  -> numu mu+
                               ///   6  K+  -> nue pi0 e+
                               ///   7  K+  -> numu pi0 mu+
                               ///   8  K-  -> numubar mu-
                               ///   9  K-  -> nuebar pi0 e-
                               ///  10  K-  -> numubar pi0 mu-
                               ///  11  mu+ -> numubar nue e+
                               ///  12  mu- -> numu nuebar e-
                               ///  13  pi+ -> numu mu+
                               ///  14  pi- -> numubar mu-
  TBranch * fBr_Ntype;         ///< Neutrino flavor (56=numu, 55=numu-bar, 53=nue, 52=nue-bar)
  TBranch * fBr_Vx;            ///< X position of hadron/muon decay
  TBranch * fBr_Vy;            ///< Y position of hadron/muon decay
  TBranch * fBr_Vz;            ///< Z position of hadron/muon decay
  TBranch * fBr_pdPx;          ///< Parent X momentum at decay point
  TBranch * fBr_pdPy;          ///< Parent Y momentum at decay point
  TBranch * fBr_pdPz;          ///< Parent Z momentum at decay point
  TBranch * fBr_ppdxdz;        ///< Parent dxdz direction at production
  TBranch * fBr_ppdydz;        ///< Parent dydz direction at production
  TBranch * fBr_pppz;          ///< Parent Z momentum at production
  TBranch * fBr_ppenergy;      ///< Parent energy at production
  TBranch * fBr_ppmedium;      ///< Tracking medium number where parent was produced
  TBranch * fBr_ptype;         ///< Parent GEANT code particle ID
  TBranch * fBr_ppvx;          ///< Parent production vertex X,Y,Z (cm)
  TBranch * fBr_ppvy;          ///< see above
  TBranch * fBr_ppvz;          ///< see above
  TBranch * fBr_muparpx;       ///< Repeat of information above, but for muon neutrino parents
  TBranch * fBr_muparpy;       ///< see above
  TBranch * fBr_muparpz;       ///< see above
  TBranch * fBr_mupare;        ///< see above
  TBranch * fBr_Necm;          ///< Neutrino energy in COM frame
  TBranch * fBr_Nimpwt;        ///< Weight of neutrino parent
  TBranch * fBr_xpoint;        ///< Debugging hook -- unused
  TBranch * fBr_ypoint;        ///< Debugging hook -- unused
  TBranch * fBr_zpoint;        ///< Debugging hook -- unused
  TBranch * fBr_tvx;           ///< X exit point of parent particle at the target
  TBranch * fBr_tvy;           ///< Y exit point of parent particle at the target
  TBranch * fBr_tvz;           ///< Z exit point of parent particle at the target
  TBranch * fBr_tpx;           ///< Parent momentum exiting the target (X)
  TBranch * fBr_tpy;           ///< Parent momentum exiting the target (Y)
  TBranch * fBr_tpz;           ///< Parent momentum exiting the target (Z)
  TBranch * fBr_tptype;        ///< Parent particle ID exiting the target (GEANT code)
  TBranch * fBr_tgen;          ///< Parent generation in cascade
                               ///  1 = primary proton
                               ///  2 = particles produced by proton interaction
                               ///  3 = particles produced by interactions of the 2's, ...
  TBranch * fBr_tgptype;       ///< Type of particle that created a particle flying of the target
  TBranch * fBr_tgppx;         ///< Momentum of a particle, that created a particle that flies off the target, at the interaction point.
  TBranch * fBr_tgppy;         ///< see above
  TBranch * fBr_tgppz;         ///< see above
  TBranch * fBr_tprivx;        ///< Primary particle interaction vertex
  TBranch * fBr_tprivy;        ///< see above
  TBranch * fBr_tprivz;        ///< see above
  TBranch * fBr_beamx;         ///< Primary proton origin
  TBranch * fBr_beamy;         ///< see above
  TBranch * fBr_beamz;         ///< see above
  TBranch * fBr_beampx;        ///< Primary proton momentum
  TBranch * fBr_beampy;        ///< see above
  TBranch * fBr_beampz;        ///< see above
  //
  // corresponding leaves are part of owned structure
  //
  GNuMIFluxPassThroughInfo* fCurrentEntry;

};


// A small persistable C-struct -like class that mirrors the structure of 
// the gnumi ntuples.  This can then be stored as an extra branch of
// the output event tree -alongside with the generated event branch- 
// for use further upstream in the analysis chain -eg beam reweighting etc-)
//
class GNuMIFluxPassThroughInfo: public TObject {
public:
   GNuMIFluxPassThroughInfo();
   /* allow default copy constructor ... for now nothing special
   GNuMIFluxPassThroughInfo(const GNuMIFluxPassThroughInfo & info);
   */
   virtual ~GNuMIFluxPassThroughInfo() { };

   void Reset();  //
   void ConvertPartCodes();

   friend ostream & operator << (ostream & stream, const GNuMIFluxPassThroughInfo & info);

   int   pcodes;  // 0=original GEANT particle codes, 1=converted to PDG
   int   units;   // 0=original GEANT cm, 1=meters

   // maintained variable names from gnumi ntuples
   // see http://www.hep.utexas.edu/~zarko/wwwgnumi/v19/[/v19/output_gnumi.html]
   int   run;
   int   evtno;
   float ndxdz;
   float ndydz;
   float npz;
   float nenergy;
   float ndxdznea;
   float ndydznea;
   float nenergyn;
   float nwtnear;
   float ndxdzfar;
   float ndydzfar;
   float nenergyf;
   float nwtfar;
   int   norig;
   int   ndecay;
   int   ntype;
   float vx;
   float vy;
   float vz;
   float pdpx;
   float pdpy;
   float pdpz;
   float ppdxdz;
   float ppdydz;
   float pppz;
   float ppenergy;
   int   ppmedium;
   int   ptype;     // converted to PDG
   float ppvx;
   float ppvy;
   float ppvz;
   float muparpx;
   float muparpy;
   float muparpz;
   float mupare;
   float necm;
   float nimpwt;
   float xpoint;
   float ypoint;
   float zpoint;
   float tvx;
   float tvy;
   float tvz;
   float tpx;
   float tpy;
   float tpz;
   int   tptype;   // converted to PDG
   int   tgen;
   int   tgptype;  // converted to PDG
   float tgppx;
   float tgppy;
   float tgppz;
   float tprivx;
   float tprivy;
   float tprivz;
   float beamx;
   float beamy;
   float beamz;
   float beampx;
   float beampy;
   float beampz;    

ClassDef(GNuMIFluxPassThroughInfo,2)
};

} // flux namespace
} // genie namespace

#endif // _GNUMI_NEUTRINO_FLUX_H_
