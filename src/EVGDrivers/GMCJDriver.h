//____________________________________________________________________________
/*!

\class    genie::GMCJDriver

\brief    GENIE MC Job Driver (event generation for the input flux & geometry)

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  May 25, 2005

\cpright  Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _GENIE_MC_JOB_DRIVER_H_
#define _GENIE_MC_JOB_DRIVER_H_

#include <string>
#include <map>

#include <TBits.h>
#include <TH1D.h>

#include "EVGDrivers/PathLengthList.h"
#include "PDG/PDGCodeList.h"

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
  void UseFluxDriver               (GFluxI * flux);
  void UseGeomAnalyzer             (GeomAnalyzerI * geom);
  void UseSplines                  (bool useLogE = true);
  void UseMaxPathLengths           (string xml_filename);
  void KeepOnThrowingFluxNeutrinos (bool keep_on);
  void FilterUnphysical            (const TBits & unphysmask);
  void ForceSingleProbScale        (void);
  void Configure                   (void);

  // generate single neutrino event for input flux & geometry
  EventRecord * GenerateEvent (void);

  // methods for enquiring info needed for computing the generated sample normalization
  double   GlobProbScale  (void) const { return fGlobPmax;                  }
  long int NFluxNeutrinos (void) const { return (long int) fNFluxNeutrinos; }

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
  string          fMaxPlXmlFilename;   ///< [config] input file with max density-weighted path lengths for all materials
  bool            fUseExtMaxPl;        ///< [config] using external max path length estimate?
  bool            fUseSplines;         ///< [config] compute all needed & not-loaded splines at init
  bool            fUseLogE;            ///< [config] build splines = f(logE) (rather than f(E)) ?
  bool            fKeepThrowingFluxNu; ///< [config] keep firing flux neutrinos till one of them interacts
  bool            fGenerateUnweighted; ///< [config] force single probability scale?
  TBits           fUnphysMask;         ///< [config] unphysical events filtering mask
};

}      // genie namespace
#endif // _GENIE_MC_JOB_DRIVER_H_
