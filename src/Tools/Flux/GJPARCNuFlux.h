//____________________________________________________________________________
/*!

\class    genie::flux::GJPARCNuFlux

\brief    A GENIE flux driver encapsulating the JPARC neutrino flux.
          It reads-in the official JPARC neutrino flux ntuples.

\ref      See http://jnusrv01.kek.jp/internal/t2k/nubeam/flux/ 
          (Note: T2K internal)

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  Feb 04, 2008

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _GJPARC_NEUTRINO_FLUX_H_
#define _GJPARC_NEUTRINO_FLUX_H_

#include <string>
#include <iostream>

#include <TLorentzVector.h>
#include <TVector3.h>
#include <TH1D.h>

#include "Framework/EventGen/GFluxI.h"
#include "Framework/ParticleData/PDGUtils.h"

class TFile;
class TTree;
class TChain;
class TBranch;

using std::string;
using std::ostream;

namespace genie {
namespace flux  {

class GJPARCNuFluxPassThroughInfo;

ostream & operator << (ostream & stream, const GJPARCNuFluxPassThroughInfo & info);

class GJPARCNuFlux: public GFluxI {

public :
  GJPARCNuFlux();
 ~GJPARCNuFlux();

  // Methods implementing the GENIE GFluxI interface, required for integrating
  // the JPARC neutrino flux simulations with the GENIE event generation drivers

  const PDGCodeList &    FluxParticles (void) { return *fPdgCList;             }
  double                 MaxEnergy     (void) { return  fMaxEv;                }
  bool                   GenerateNext  (void);
  int                    PdgCode       (void) { return  fgPdgC;                }
  double                 Weight        (void) { return  fNorm / fMaxWeight;    }
  const TLorentzVector & Momentum      (void) { return  fgP4;                  }
  const TLorentzVector & Position      (void) { return  fgX4;                  }
  bool                   End           (void) { return  fEntriesThisCycle >= fNEntries 
                                                     && fICycle == fNCycles && fNCycles > 0;   }
  long int               Index         (void);                              
  void                   Clear            (Option_t * opt); 
  void                   GenerateWeighted (bool gen_weighted = true);

  // Methods specific to the JPARC flux driver, 
  // for configuration/initialization of the flux & event generation drivers and
  // and for passing-through flux information (eg neutrino parent decay kinematics)
  // not used by the generator but required by analyses/processing further upstream

  bool LoadBeamSimData  (string filename, string det_loc);     ///< load a jnubeam root flux ntuple
  void SetFluxParticles (const PDGCodeList & particles);       ///< specify list of flux neutrino species
  void SetMaxEnergy     (double Ev);                           ///< specify maximum flx neutrino energy
  void SetFilePOT       (double pot);                          ///< flux file norm is in /N POT/det [ND] or /N POT/cm^2 [FD]. Specify N (typically 1E+21)
  void SetUpstreamZ     (double z0);                           ///< set flux neutrino initial z position (upstream of the detector)
  void SetNumOfCycles   (int n);                               ///< set how many times to cycle through the ntuple (default: 1 / n=0 means 'infinite')
  void DisableOffset    (void){fUseRandomOffset = false;}      ///< switch off random offset, must be called before LoadBeamSimData to have any effect 
  void RandomOffset     (void);                                ///< choose a random offset as starting entry in flux ntuple 

  double   POT_1cycle     (void);                              ///< flux POT per cycle
  double   POT_curravg    (void);                              ///< current average POT
  long int NFluxNeutrinos (void) const { return fNNeutrinos; } ///< number of flux neutrinos looped so far
  double   SumWeight      (void) const { return fSumWeight;  } ///< intergated weight for flux neutrinos looped so far

  const GJPARCNuFluxPassThroughInfo & 
     PassThroughInfo(void) { return *fPassThroughInfo; } ///< GJPARCNuFluxPassThroughInfo

private:

  // Private methods
  //
  bool GenerateNext_weighted (void);
  void Initialize            (void);
  void SetDefaults           (void);  
  void CleanUp               (void);
  void ResetCurrent          (void);
  int  DLocName2Id           (string name);

  // Private data members
  //
  double         fMaxEv;       ///< maximum energy
  PDGCodeList *  fPdgCList;    ///< list of neutrino pdg-codes

  int            fgPdgC;       ///< running generated nu pdg-code
  TLorentzVector fgP4;         ///< running generated nu 4-momentum
  TLorentzVector fgX4;         ///< running generated nu 4-position

  TFile *   fNuFluxFile;       ///< input flux file
  TTree *   fNuFluxTree;       ///< input flux ntuple
  TChain *  fNuFluxChain;      ///< input flux ntuple
  TTree *   fNuFluxSumTree;    ///< input summary ntuple for flux file. This tree is only present for later flux versions > 10a
  TChain *  fNuFluxSumChain;   ///< input summary ntuple for flux file. This tree is only present for later flux versions > 10a
  bool      fNuFluxUsingTree;  ///< are we using a TTree or a TChain to view the input flux file?
  string    fDetLoc;           ///< input detector location ('sk','nd1','nd2',...)
  int       fDetLocId;         ///< input detector location id (fDetLoc -> jnubeam idfd)
  int       fNDetLocIdFound;   ///< per cycle keep track of the number of fDetLocId are found - if this is zero will exit job 
  bool      fIsFDLoc;          ///< input location is a 'far'  detector location?
  bool      fIsNDLoc;          ///< input location is a 'near' detector location?
  long int  fNEntries;         ///< number of flux ntuple entries
  long int  fIEntry;           ///< current flux ntuple entry
  long int  fEntriesThisCycle; ///< keep track of number of entries used so far for this cycle   
  long int  fOffset;           ///< start looping at entry fOffset
  double    fNorm;             ///< current flux ntuple normalisation
  double    fMaxWeight;        ///< max flux  neutrino weight in input file for the specified detector location
  double    fFilePOT;          ///< file POT normalization, typically 1E+21
  double    fZ0;               ///< configurable starting z position for each flux neutrino (in detector coord system)
  int       fNCycles;          ///< how many times to cycle through the flux ntuple
  int       fICycle;           ///< current cycle
  double    fSumWeight;        ///< sum of weights for neutrinos thrown so far
  long int  fNNeutrinos;       ///< number of flux neutrinos thrown so far
  double    fSumWeightTot1c;   ///< total sum of weights for neutrinos to be thrown / cycle
  long int  fNNeutrinosTot1c;  ///< total number of flux neutrinos to be thrown / cycle
  bool      fGenerateWeighted; ///< generate weighted/deweighted flux neutrinos (default is false)
  bool      fUseRandomOffset;  ///< whether set random starting point when looping over flux ntuples
  bool      fLoadedNeutrino;   ///< set to true when GenerateNext_weighted has been called successfully

  GJPARCNuFluxPassThroughInfo * fPassThroughInfo;
};


// A small persistable C-struct - like class that may be stored at an extra branch of 
// the output event tree -alongside with the generated event branch- for use further 
// upstream in the t2k analysis chain -eg beam reweighting etc-)
//
const int fNgmax = 12;
class GJPARCNuFluxPassThroughInfo: public TObject { 
public:
   GJPARCNuFluxPassThroughInfo();
   GJPARCNuFluxPassThroughInfo(const GJPARCNuFluxPassThroughInfo & info);
   virtual ~GJPARCNuFluxPassThroughInfo() { };

   void Reset();
   friend ostream & operator << (ostream & stream, const GJPARCNuFluxPassThroughInfo & info);

   long   fluxentry;
   string fluxfilename;
   // Using an instance of this class the following datamembers are set 
   // directly as the branch addresses of jnubeam flux ntuples tree:
   float  Enu;            // set to "Enu/F": Nu energy (GeV)  
   int    ppid;           // set to "ppid/I": Nu parent GEANT particle id 
   int    mode;           // set to "mode/I": Nu parent decay mode (see http://jnusrv01.kek.jp/internal/t2k/nubeam/flux/nemode.h)
   float  ppi;            // set to "ppi/F": Nu parent momentum at its decay point (GeV)
   float  xpi[3];         // set to "xpi[3]/F": Nu parent position vector at decay (cm, in t2k global coord system)
   float  npi[3];         // set to "npi[3]/F": Nu parent direction vector at decay (in t2k global coord system) 
   float  norm;           // set to "norm/F": Weight to give flux in /N POT/det. [ND] or /N POT/cm^2 [FD], where is N is typically 1E+21
   int    nvtx0;          // set to "nvtx0/I": Number of vtx where the nu. parent was produced (made obsolete by nd variable inroduced in 10d flux version)
   float  ppi0;           // set to "ppi0/F": Nu parent momentum at its production point (GeV)
   float  xpi0[3];        // set to "xpi0[3]/F": Nu parent position vector at production (cm, in t2k global coord system)
   float  npi0[3];        // set to "npi0[3]/F": Nu parent direction vector at production (in t2k global coord system)
   float  rnu;            // set to "rnu/F": Nu radial position (cm, in detector coord system)
   float  xnu;            // set to "xnu/F": Nu x position (cm, in detector coord system)
   float  ynu;            // set to "ynu/F": Nu y position (cm, in detector coord system)
   float  nnu[3];         // set to "nnu[3]/F": Nu direction (in t2k global coord system)
   // New since JNuBeam 10a flux version.
   float  cospibm;        // set to "cospibm/F": Nu parent direction cosine at decay (with respect to the beam direction) 
   float  cospi0bm;       // set to "cospi0bm/F": Nu parent direction cosine at production (with respect to the beam direction)
   int    idfd;           // set to "idfd/I": Detector ID
   unsigned char gipart;  // set to "gipart/B": Primary particle ID
   float  gpos0[3];       // set to "gpos0[3]/F": Primary particle starting point
   float  gvec0[3];       // set to "gvec0[3]/F": Primary particle direction at the starting point
   float  gamom0;         // set to "gamom0/F": Momentum of the primary particle at the starting point
   // New since JNuBeam 10d and 11a flux version updates
   int    ng;             // set to "ng/I": Number of parents (number of generations)
   float  gpx[fNgmax];    // set to "gpx[20]/F":  Momentum X of each ancestor particle
   float  gpy[fNgmax];    // set to "gpy[20]/F":  Momentum Y of each ancestor particle
   float  gpz[fNgmax];    // set to "gpz[20]/F":  Momentum Z of each ancestor particle
   float  gcosbm[fNgmax]; // set to "gcosbm[20]/F": Cosine of the angle between the ancestor particle direction and the beam direction
   float  gvx[fNgmax];    // set to "gvx[20]/F": Vertex X position of each ancestor particle 
   float  gvy[fNgmax];    // set to "gvy[20]/F": Vertex Y position of each ancestor particle 
   float  gvz[fNgmax];    // set to "gvz[20]/F": Vertex Z position of each ancestor particle 
   int    gpid[fNgmax];   // set to "gpid[20]/I": Particle ID of each ancestor particles
   int    gmec[fNgmax];   // set to "gmec[20]/I": Particle production mechanism of each ancestor particle 
   // Next five only present since 11a flux
   int    gmat[fNgmax];   // set to "gmat[fNgmax]/I": Material in which the particle originates 
   float  gdistc[fNgmax]; // set to "gdistc[fNgmax]/F": Distance traveled through carbon 
   float  gdistal[fNgmax]; // set to "gdista[fNgmax]/F": Distance traveled through aluminum
   float  gdistti[fNgmax];// set to "gdistti[fNgmax]/F": Distance traveled through titanium
   float  gdistfe[fNgmax];// set to "gdistte[fNgmax]/F": Distance traveled through iron
   float  Enusk;          // set to "Enusk/F": "Enu" for SK
   float  normsk;         // set to "normsk/F": "norm" for SK 
   float  anorm;          // set to "anorm/F": Norm component from ND acceptance calculation
   // The following do not change per flux entry as is summary info for the flux
   // file. For simplicity we just store per flux entry and accept the duplication.
   float  version;        // set to "version/F": Jnubeam version
   int    tuneid;         // set to "tuneid/I": Parameter set identifier
   int    ntrig;          // set to "ntrig/I": Number of Triggers in simulation
   int    pint;           // set to "pint/I": Interaction model ID
   float  bpos[2];        // set to "bpos[2]/F": Beam center position
   float  btilt[2];       // set to "btilt[2]/F": Beam Direction
   float  brms[2];        // set to "brms[2]/F": Beam RMS Width
   float  emit[2];        // set to "emit[2]/F": Beam Emittance 
   float  alpha[2];       // set to "alpha[2]/F": Beam alpha parameter
   float  hcur[3];        // set to "hcur[3]/F": Horns 1, 2 and 3 Currents
   int    rand;           // set to "rand/I": Random seed
   int    rseed[2];          // set to "rseed/I": Random seed

ClassDef(GJPARCNuFluxPassThroughInfo,3)
};

} // flux namespace
} // genie namespace

#endif // _GJPARC_NEUTRINO_FLUX_H_
