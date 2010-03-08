//____________________________________________________________________________
/*!

\class    genie::flux::GJPARCNuFlux

\brief    A GENIE flux driver encapsulating the JPARC neutrino flux.
          It reads-in the official JPARC neutrino flux ntuples.

\ref      See http://jnusrv01.kek.jp/internal/t2k/nubeam/flux/ 
          (Note: T2K internal)

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  Feb 04, 2008

\cpright  Copyright (c) 2003-2010, GENIE Neutrino MC Generator Collaboration
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

#include "EVGDrivers/GFluxI.h"
#include "PDG/PDGUtils.h"

class TFile;
class TTree;
class TBranch;

using std::string;
using std::ostream;

namespace genie {
namespace flux  {

class GJPARCNuFluxPassThroughInfo;

class GJPARCNuFlux: public GFluxI {

public :
  GJPARCNuFlux();
 ~GJPARCNuFlux();

  // Methods implementing the GENIE GFluxI interface, required for integrating
  // the JPARC neutrino flux simulations with the GENIE event generation drivers

  const PDGCodeList &    FluxParticles (void) { return *fPdgCList;            }
  double                 MaxEnergy     (void) { return  fMaxEv;               }
  bool                   GenerateNext  (void);
  int                    PdgCode       (void) { return  fgPdgC;               }
  double                 Weight        (void) { return  fLfNorm;              }
  const TLorentzVector & Momentum      (void) { return  fgP4;                 }
  const TLorentzVector & Position      (void) { return  fgX4;                 }
  bool                   End           (void) { return  fIEntry >= fNEntries 
                                                     && fICycle == fNCycles;  }

  // Methods specific to the JPARC flux driver, 
  // for configuration/initialization of the flux & event generation drivers and
  // and for passing-through flux information (eg neutrino parent decay kinematics)
  // not used by the generator but required by analyses/processing further upstream

  void LoadBeamSimData  (string filename, string det_loc);     ///< load a jnubeam root flux ntuple
  void SetFluxParticles (const PDGCodeList & particles);       ///< specify list of flux neutrino species
  void SetMaxEnergy     (double Ev);                           ///< specify maximum flx neutrino energy
  void SetFilePOT       (double pot);                          ///< flux file norm is in /N POT/det [ND] or /N POT/cm^2 [FD]. Specify N (typically 1E+21)
  void SetUpstreamZ     (double z0);                           ///< set flux neutrino initial z position (upstream of the detector)
  void SetNumOfCycles   (int n);                               ///< set how many times to cycle through the ntuple (default: 1 / n=0 means 'infinite')

  double   POT_1cycle     (void);                              ///< flux POT per cycle
  double   POT_curravg    (void);                              ///< current average POT
  long int NFluxNeutrinos (void) const { return fNNeutrinos; } ///< number of flux neutrinos looped so far
  double   SumWeight      (void) const { return fSumWeight;  } ///< intergated weight for flux neutrinos looped so far
  string   FluxVersion    (void) const { return fFluxVersion;} ///< return the string describing the current flux version


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
  string    fDetLoc;           ///< input detector location ('sk','nd1','nd2',...)
  int       fDetLocId;         ///< input detector location id (fDetLoc -> jnubeam idfd)
  int       fNDetLocIdFound;   ///< per cycle keep track of the number of fDetLocId are found - if this is zero will exit job 
  bool      fIsFDLoc;          ///< input location is a 'far'  detector location?
  bool      fIsNDLoc;          ///< input location is a 'near' detector location?
  long int  fNEntries;         ///< number of flux ntuple entries
  long int  fIEntry;           ///< current flux ntuple entry
  double    fMaxWeight;        ///< max flux  neutrino weight in input file for the specified detector location
  double    fFilePOT;          ///< file POT normalization, typically 1E+21
  double    fZ0;               ///< configurable starting z position for each flux neutrino (in detector coord system)
  int       fNCycles;          ///< how many times to cycle through the flux ntuple
  int       fICycle;           ///< current cycle
  double    fSumWeight;        ///< sum of weights for neutrinos thrown so far
  long int  fNNeutrinos;       ///< number of flux neutrinos thrown so far
  double    fSumWeightTot1c;   ///< total sum of weights for neutrinos to be thrown / cycle
  long int  fNNeutrinosTot1c;  ///< total number of flux neutrinos to be thrown / cycle

  string    fFluxVersion;      ///< string representing current flux version, e.g., "07a", "10a", etc...

  //-- jnubeam ntuple branches
  //   branches marked with [f] can be found in SK flux ntuples only
  //   branches marked with [n] can be found in near detector flux ntuples only
  //   branches marked with [a] can be found in both ntuples
  //   branches marked with [10a] first appeared in JNUBEAM version 10a flux ntuples
  TBranch * fBrNorm;           ///< 'norm'     branch [a]: Weight to give flux in /N POT/det. [ND] or /N POT/cm^2 [FD], where is N is typically 1E+21
  TBranch * fBrIdfd;           ///< 'idfd'     branch [n]: Detector ID
  TBranch * fBrEnu;            ///< 'Enu'      branch [a]: Nu energy (GeV)
  TBranch * fBrRnu;            ///< 'rnu'      branch [n]: Nu radial position (cm, in detectro coord system)
  TBranch * fBrXnu;            ///< 'xnu'      branch [n]: Nu x position (cm, in detector coord system)
  TBranch * fBrYnu;            ///< 'ynu'      branch [n]: Nu y position (cm, in detector coord system)
  TBranch * fBrNnu;            ///< 'nnu'      branch [n]: Nu direction (in t2k global coord system)
  TBranch * fBrPpid;           ///< 'ppid'     branch [a]: Nu parent GEANT particle id 
  TBranch * fBrMode;           ///< 'mode'     branch [a]: Nu parent decay mode (see http://jnusrv01.kek.jp/internal/t2k/nubeam/flux/nemode.h)
  TBranch * fBrPpi;            ///< 'ppi'      branch [a]: Nu parent momentum at its decay point (GeV)
  TBranch * fBrXpi;            ///< 'xpi'      branch [a]: Nu parent position vector at decay (cm, in t2k global coord system)
  TBranch * fBrNpi;            ///< 'npi'      branch [a]: Nu parent direction vector at decay (in t2k global coord system) 
  TBranch * fBrCospibm;        ///< 'cospibm'  branch [a]: Nu parent direction cosine at decay (with respect to the beam direction)
  TBranch * fBrPpi0;           ///< 'ppi0'     branch [a]: Nu parent momentum at its production point (GeV)
  TBranch * fBrXpi0;           ///< 'xpi0'     branch [a]: Nu parent position vector at production (cm, in t2k global coord system)
  TBranch * fBrNpi0;           ///< 'npi0'     branch [a]: Nu parent direction vector at production (in t2k global coord system)
  TBranch * fBrCospi0bm;       ///< 'cospi0bm' branch [n]: Nu parent direction cosine at production (with respect to the beam direction)
  TBranch * fBrNVtx0;          ///< 'nvtx0'    branch [f]: Number of vtx where the nu. parent was produced
  TBranch * fBrGipart;         ///< 'gipart'   branch [a][10a]: Primary particle ID
  TBranch * fBrGpos0;          ///< 'gpos0'    branch [a][10a]: Primary particle starting point
  TBranch * fBrGvec0;          ///< 'gvec0'    branch [a][10a]: Primary particle direction at the starting point
  TBranch * fBrGamom0;         ///< 'gamom0'   branch [a][10a]: Momentum of the primary particle at the starting point

  float     fLfNorm;           ///< leaf on branch 'norm'
  int       fLfIdfd;           ///< leaf on branch 'idfd'
  float     fLfEnu;            ///< leaf on branch 'Enu'
  float     fLfRnu;            ///< leaf on branch 'rnu'
  float     fLfXnu;            ///< leaf on branch 'xnu'
  float     fLfYnu;            ///< leaf on branch 'ynu'
  float     fLfNnu[3];         ///< leaf on branch 'nnu'
  int       fLfPpid;           ///< leaf on branch 'ppid'
  int       fLfMode;           ///< leaf on branch 'mode'
  float     fLfPpi;            ///< leaf on branch 'ppi'
  float     fLfXpi[3];         ///< leaf on branch 'xpi'
  float     fLfNpi[3];         ///< leaf on branch 'npi'
  float     fLfCospibm;        ///< leaf on branch 'cospibm'
  float     fLfPpi0;           ///< leaf on branch 'ppi0'
  float     fLfXpi0[3];        ///< leaf on branch 'xpi0'
  float     fLfNpi0[3];        ///< leaf on branch 'npi0'
  float     fLfCospi0bm;       ///< leaf on branch 'cospi0bm'
  int       fLfNVtx0;          ///< leaf on branch 'nvtx0'
  float     fLfGpos0[3];       ///< leaf on branch 'gpos0'
  float     fLfGvec0[3];       ///< leaf on branch 'gvec0'
  float     fLfGamom0;         ///< leaf on branch 'gamom0'
  unsigned char fLfGipart;     ///< leaf on branch 'gipart'
  
  GJPARCNuFluxPassThroughInfo * fPassThroughInfo;
};


// A small persistable C-struct - like class that may be stored at an extra branch of 
// the output event tree -alongside with the generated event branch- for use further 
// upstream in the t2k analysis chain -eg beam reweighting etc-)
//
class GJPARCNuFluxPassThroughInfo: public TObject { 
public:
   GJPARCNuFluxPassThroughInfo();
   GJPARCNuFluxPassThroughInfo(const GJPARCNuFluxPassThroughInfo & info);
   virtual ~GJPARCNuFluxPassThroughInfo() { };

   friend ostream & operator << (ostream & stream, const GJPARCNuFluxPassThroughInfo & info);

   long   fluxentry;
   int    ppid;
   int    mode;
   float  ppi;
   float  xpi[3];
   float  npi[3];
   float  ppi0;
   float  xpi0[3];
   float  npi0[3];
   int    nvtx0;
   float  cospibm;
   float  cospi0bm;
   int    idfd;
   float  gpos0[3];
   float  gvec0[3];
   float  gamom0;
   int    gipart; 

ClassDef(GJPARCNuFluxPassThroughInfo,2)
};

} // flux namespace
} // genie namespace

#endif // _GJPARC_NEUTRINO_FLUX_H_
