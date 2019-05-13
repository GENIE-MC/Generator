//____________________________________________________________________________
/*!

\class    genie::GMCJDriver

\brief    A GENIE `MC Job Driver'. Can be used for setting up complicated event 
          generation cases involving detailed flux descriptions and detector 
          geometry descriptions.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  May 25, 2005

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _GENIE_MC_JOB_DRIVER_H_
#define _GENIE_MC_JOB_DRIVER_H_

#include <string>
#include <map>

#include <TH1D.h>
#include <TLorentzVector.h>
#include <TFile.h>
#include <TTree.h>
#include <TBits.h>

#include "Framework/EventGen/PathLengthList.h"
#include "Framework/ParticleData/PDGCodeList.h"

using std::string;
using std::map;

namespace genie {

class EventRecord;
class GFluxI;
class GeomAnalyzerI;
class GENIE;
class GEVGPool;

class GMCJDriver {

public :
  GMCJDriver();
 ~GMCJDriver();

  // configure MC job
  void SetEventGeneratorList       (string listname);
  void SetUnphysEventMask          (const TBits & mask);
  void UseFluxDriver               (GFluxI * flux);
  void UseGeomAnalyzer             (GeomAnalyzerI * geom);
  void UseSplines                  (bool useLogE = true);
  bool UseMaxPathLengths           (string xml_filename);
  void KeepOnThrowingFluxNeutrinos (bool keep_on);
  void ForceSingleProbScale        (void);
  void PreSelectEvents             (bool preselect = true);
  bool PreCalcFluxProbabilities    (void);
  bool LoadFluxProbabilities       (string filename);
  void SaveFluxProbabilities       (string outfilename);
  void Configure                   (bool calc_prob_scales = true);

  // generate single neutrino event for input flux & geometry
  EventRecord * GenerateEvent (void);

  // info needed for computing the generated sample normalization
  double   GlobProbScale  (void) const { return fGlobPmax;                  }
  long int NFluxNeutrinos (void) const { return (long int) fNFluxNeutrinos; }
  map<int, double> SumFluxIntProbs(void) const { return fSumFluxIntProbs;   }

  // input flux and geometry drivers
  const GFluxI &        FluxDriver      (void) const { return *fFluxDriver;   }
  const GeomAnalyzerI & GeomAnalyzer    (void) const { return *fGeomAnalyzer; }
  GFluxI *              FluxDriverPtr   (void) const { return  fFluxDriver;   } 
  GeomAnalyzerI *       GeomAnalyzerPtr (void) const { return  fGeomAnalyzer; }

private:
 
  // private methods:
  void          InitJob                         (void);
  void          InitEventGeneration             (void);
  void          GetParticleLists                (void);
  void          GetMaxPathLengthList            (void);
  void          GetMaxFluxEnergy                (void);
  void          PopulateEventGenDriverPool      (void);
  void          BootstrapXSecSplines            (void);
  void          BootstrapXSecSplineSummation    (void);
  void          ComputeProbScales               (void);
  EventRecord * GenerateEvent1Try               (void);
  bool          GenerateFluxNeutrino            (void);
  bool          ComputePathLengths              (void);
  double	ComputeInteractionProbabilities (bool use_max_path_length);
  int           SelectTargetMaterial            (double R);
  void          GenerateEventKinematics         (void);
  void          GenerateVertexPosition          (void);
  void          ComputeEventProbability         (void);
  double        InteractionProbability          (double xsec, double pl, int A);
  double        PreGenFluxInteractionProbability(void);

  // private data members:
  GEVGPool *      fGPool;              ///< A pool of GEVGDrivers properly configured event generation drivers / one per init state
  GFluxI *        fFluxDriver;         ///< [input] neutrino flux driver
  GeomAnalyzerI * fGeomAnalyzer;       ///< [input] detector geometry analyzer
  double          fEmax;               ///< [declared by the flux driver] maximum neutrino energy 
  PDGCodeList     fNuList;             ///< [declared by the flux driver] list of neutrino codes 
  PDGCodeList     fTgtList;            ///< [declared by the geom driver] list of target codes 
  PathLengthList  fMaxPathLengths;     ///< [declared by the geom driver] maximum path length list 
  PathLengthList  fCurPathLengths;     ///< [current] path length list for current flux neutrino
  TLorentzVector  fCurVtx;             ///< [current] interaction vertex
  EventRecord *   fCurEvt;             ///< [current] generated event
  int             fSelTgtPdg;          ///< [current] selected target material PDG code
  map<int,double> fCurCumulProbMap;    ///< [current] cummulative interaction probabilities
  double          fNFluxNeutrinos;     ///< [current] number of flux nuetrinos fired by the flux driver so far 
  map<int,TH1D*>  fPmax;               ///< [computed at init] interaction probability scale /neutrino /energy for given geometry
  double          fGlobPmax;           ///< [computed at init] global interaction probability scale for given flux & geometry
  string          fEventGenList;       ///< [config] list of event generators loaded by this driver (what used to be the $GEVGL setting)
  TBits *         fUnphysEventMask;    ///< [config] controls whether unphysical events are returned (what used to be the $GUNPHYSMASK setting)
  string          fMaxPlXmlFilename;   ///< [config] input file with max density-weighted path lengths for all materials
  bool            fUseExtMaxPl;        ///< [config] using external max path length estimate?
  bool            fUseSplines;         ///< [config] compute all needed & not-loaded splines at init
  bool            fUseLogE;            ///< [config] build splines = f(logE) (rather than f(E)) ?
  bool            fKeepThrowingFluxNu; ///< [config] keep firing flux neutrinos till one of them interacts
  bool            fGenerateUnweighted; ///< [config] force single probability scale?
  bool            fPreSelect;          ///< [config] set whether to pre-select events using max interaction paths 
  TFile*          fFluxIntProbFile;    ///< [input] pre-generated flux interaction probability file
  TTree*          fFluxIntTree;        ///< [computed-or-loaded] pre-computed flux interaction probabilities (expected tree name is "gFlxIntProbs")
  double          fBrFluxIntProb;      ///< flux interaction probability (set to branch:"FluxIntProb")
  int             fBrFluxIndex;        ///< corresponding entry in flux input tree (set to address of branch:"FluxEntry")
  double          fBrFluxEnu;          ///< corresponding flux P4 (set to address of branch:"FluxP4") 
  double          fBrFluxWeight;       ///< corresponding flux weight (set to address of branch: "FluxWeight") 
  int             fBrFluxPDG;          ///< corresponding flux pdg code (set to address of branch: "FluxPDG") 
  string          fFluxIntFileName;    ///< whether to save pre-generated flux tree for use in later jobs
  string          fFluxIntTreeName;    ///< name for tree holding flux probabilities 
  map<int, double> fSumFluxIntProbs;   ///< map where the key is flux pdg code and the value is sum of fBrFluxWeight * fBrFluxIntProb for all these flux neutrinos 
};

}      // genie namespace
#endif // _GENIE_MC_JOB_DRIVER_H_
