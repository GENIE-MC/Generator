//____________________________________________________________________________
/*!

\class    genie::flux::GSimpleNtpFlux

\brief    A GENIE flux driver using a simple ntuple format

\author   Robert Hatcher <rhatcher \at fnal.gov>
          Fermi National Accelerator Laboratory

\created  Jan 25, 2010

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _SIMPLE_NTP_FLUX_H_
#define _SIMPLE_NTP_FLUX_H_

#include <string>
#include <iostream>
#include <vector>
#include <set>

#include <TVector3.h>
#include <TLorentzVector.h>
#include <TLorentzRotation.h>

#include "Framework/EventGen/GFluxI.h"
#include "Tools/Flux/GFluxExposureI.h"
#include "Tools/Flux/GFluxFileConfigI.h"
#include "Framework/ParticleData/PDGUtils.h"

class TFile;
class TChain;
class TTree;
class TBranch;

using std::string;
using std::ostream;

namespace genie {
namespace flux  {

class GSimpleNtpEntry;    
ostream & operator << (ostream & stream, const GSimpleNtpEntry & info);

/// Small persistable C-struct -like classes that makes up the SimpleNtpFlux
/// ntuple.  This is only valid for a particular flux window (no reweighting,
/// no coordinate transformation available).
///
/// Order elements from largest to smallest for ROOT alignment purposes

/// GSimpleNtpEntry
/// =========================
/// This is the only required branch ("entry") of the "flux" tree
  class GSimpleNtpEntry {
  public:
    GSimpleNtpEntry();
    /* allow default copy constructor ... for now nothing special
       GSimpleNtpEntry(const GSimpleNtpEntry & info);
    */
    virtual ~GSimpleNtpEntry() { };
    void Reset();
    void Print(const Option_t* opt = "") const;
    friend ostream & operator << (ostream & stream, const GSimpleNtpEntry & info);
    
    Double_t   wgt;      ///< nu weight
    
    Double_t   vtxx;     ///< x position in lab frame
    Double_t   vtxy;     ///< y position in lab frame
    Double_t   vtxz;     ///< z position in lab frame
    Double_t   dist;     ///< distance from hadron decay
    
    Double_t   px;       ///< x momentum in lab frame
    Double_t   py;       ///< y momentum in lab frame
    Double_t   pz;       ///< z momentum in lab frame
    Double_t   E;        ///< energy in lab frame
    
    Int_t      pdg;      ///< nu pdg-code
    UInt_t     metakey;  ///< key to meta data

    ClassDef(GSimpleNtpEntry,1)
  };


  class GSimpleNtpNuMI;
  ostream & operator << (ostream & stream, const GSimpleNtpNuMI & info);

/// GSimpleNtpNuMI
/// =========================
/// Additional elements for NuMI (allow SKZP reweighting and reference
/// back to original GNuMI flux entries) as "numi" branch
  class GSimpleNtpNuMI {
  public:
    GSimpleNtpNuMI();
    virtual ~GSimpleNtpNuMI() { };
    void Reset();
    void Print(const Option_t* opt = "") const;
    friend ostream & operator << (ostream & stream, const GSimpleNtpNuMI & info);
    Double_t   tpx;      ///< parent particle px at target exit
    Double_t   tpy;
    Double_t   tpz;
    Double_t   vx;       ///< vertex position of hadron/muon decay
    Double_t   vy;
    Double_t   vz;
    Double_t   pdpx;     ///< nu parent momentum at time of decay
    Double_t   pdpy;
    Double_t   pdpz;
    Double_t   pppx;     ///< nu parent momentum at production point
    Double_t   pppy;
    Double_t   pppz;

    Int_t      ndecay;   ///< decay mode
    Int_t      ptype;    ///< parent type (PDG)
    Int_t      ppmedium; ///< tracking medium where parent was produced
    Int_t      tptype;   ///< parent particle type at target exit
    
    Int_t      run;      ///< 
    Int_t      evtno;    ///<
    Int_t      entryno;  ///<
    
    ClassDef(GSimpleNtpNuMI,3)
  };


  class GSimpleNtpAux;
  ostream & operator << (ostream & stream, const GSimpleNtpAux & info);

/// GSimpleNtpAux
/// =========================
/// Additional elements for expansion as "aux" branch
  class GSimpleNtpAux {
  public:
    GSimpleNtpAux();
    virtual ~GSimpleNtpAux() { };
    void Reset();
    void Print(const Option_t* opt = "") const;
    friend ostream & operator << (ostream & stream, const GSimpleNtpAux & info);

   std::vector<Int_t>    auxint;  ///< additional ints associated w/ entry
   std::vector<Double_t> auxdbl;  ///< additional doubles associated w/ entry
  
   ClassDef(GSimpleNtpAux,1)
  };


  class GSimpleNtpMeta;
  ostream & operator << (ostream & stream, const GSimpleNtpMeta & info);

/// GSimpleNtpMeta
/// =========================
/// A small persistable C-struct -like class that holds metadata 
/// about the the SimpleNtpFlux ntple.
///
  class GSimpleNtpMeta: public TObject {
  public:
    GSimpleNtpMeta();
    /* allow default copy constructor ... for now nothing special
       GSimpleNtpMeta(const GSimpleNtpMeta & info);
    */
    virtual ~GSimpleNtpMeta();

    void Reset();
    void AddFlavor(Int_t nupdg);
    void Print(const Option_t* opt = "") const;
    friend ostream & operator << (ostream & stream, const GSimpleNtpMeta & info);
    
    std::vector<Int_t>  pdglist; ///< list of neutrino flavors
   
    Double_t maxEnergy;   ///< maximum energy
    Double_t minWgt;      ///< minimum weight
    Double_t maxWgt;      ///< maximum weight
    Double_t protons;     ///< represented number of protons-on-target

    Double_t windowBase[3]; ///< x,y,z position of window base point
    Double_t windowDir1[3]; ///< dx,dy,dz of window direction 1
    Double_t windowDir2[3]; ///< dx,dy,dz of window direction 2
    
    std::vector<std::string>    auxintname;  ///< tagname of aux ints associated w/ entry
    std::vector<std::string>    auxdblname;  ///< tagname of aux doubles associated w/ entry
    std::vector<std::string>    infiles; ///< list of input files
    
    Int_t    seed;     ///< random seed used in generation
    UInt_t   metakey;  ///< index key to tie to individual entries

    static UInt_t mxfileprint;  ///< allow user to limit # of files to print

    ClassDef(GSimpleNtpMeta,1)
  };


/// GSimpleNtpFlux:
/// ==========
/// An implementation of the GFluxI interface that provides NuMI flux
///
class GSimpleNtpFlux 
  : public genie::GFluxI
  , public genie::flux::GFluxExposureI
  , public genie::flux::GFluxFileConfigI 
{

public :
  GSimpleNtpFlux();
 ~GSimpleNtpFlux();

  // Methods implementing the GENIE GFluxI interface, required for integrating
  // the NuMI neutrino flux simulations with the GENIE event generation drivers

  const PDGCodeList &    FluxParticles (void) { return *fPdgCList;            }
  double                 MaxEnergy     (void) { return  fMaxEv;               }
  bool                   GenerateNext  (void);
  int                    PdgCode       (void) { return  fCurEntry->pdg;       }
  double                 Weight        (void) { return  fWeight;              }
  const TLorentzVector & Momentum      (void) { return  fP4;  }
  const TLorentzVector & Position      (void) { return  fX4;  }
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
  const genie::flux::GSimpleNtpEntry *
    GetCurrentEntry(void) { return fCurEntry; } ///< GSimpleNtpEntry
  const genie::flux::GSimpleNtpNuMI *
    GetCurrentNuMI(void)  { return fCurNuMI; }  ///< GSimpleNtpNuMI
  const genie::flux::GSimpleNtpAux *
    GetCurrentAux(void)   { return fCurAux; }   ///< GSimpleNtpAux
  const genie::flux::GSimpleNtpMeta *
    GetCurrentMeta(void)  { return fCurMeta; }  ///< GSimpleNtpMeta

  // allow access to main tree so we can call Branch() to retrieve extra stuff
  TChain*
    GetFluxTChain(void) { return fNuFluxTree; } ///< 

  double    GetDecayDist() const; ///< dist (user units) from dk to current pos
  void      MoveToZ0(double z0);  ///< move ray origin to user coord Z0

  //
  // information about the current state
  //
  virtual double    GetTotalExposure() const;  ///< GFluxExposureI interface
  virtual long int  NFluxNeutrinos() const;    ///< # of rays generated

  double    UsedPOTs(void) const;       ///< # of protons-on-target used

  long int  NEntriesUsed(void) const { return fNEntriesUsed; } ///< number of entries read from files
  double    SumWeight(void) const { return fSumWeight;  } ///< integrated weight for flux neutrinos looped so far

  void      PrintCurrent(void);         ///< print current entry from leaves
  void      PrintConfig();              ///< print the current configuration

  std::vector<std::string> GetFileList();  ///< list of files currently part of chain

  // 
  // GFluxFileConfigI interface
  //
  virtual void  LoadBeamSimData(const std::vector<string>& filenames,
                                const std::string&         det_loc);
  using GFluxFileConfigI::LoadBeamSimData; // inherit the rest
  virtual void  GetBranchInfo(std::vector<std::string>& branchNames,
                              std::vector<std::string>& branchClassNames,
                              std::vector<void**>&      branchObjPointers);
  virtual TTree* GetMetaDataTree();

  //
  // configuration of GSimpleNtpFlux
  //

  void      SetRequestedBranchList(string blist="entry,numi,aux") { fNuFluxBranchRequest = blist; }

  void      SetMaxEnergy(double Ev);                              ///< specify maximum flx neutrino energy

  void      SetGenWeighted(bool genwgt=false) { fGenWeighted = genwgt; } ///< toggle whether GenerateNext() returns weight=1 flux (initial default false)

  void      SetEntryReuse(long int nuse=1);                       ///<  # of times to use entry before moving to next

  void      ProcessMeta(void);  ///< scan for max flux energy, weight

  void      GetFluxWindow(TVector3& p1, TVector3& p2, TVector3& p3) const; ///< 3 points define a plane in beam coordinate 

private:

  // Private methods
  //
  bool GenerateNext_weighted (void);
  void Initialize            (void);
  void SetDefaults           (void);
  void CleanUp               (void);
  void ResetCurrent          (void);
  void AddFile               (TTree* fluxtree, TTree* metatree, string fname);
  bool OptionalAttachBranch  (std::string bname);
  void CalcEffPOTsPerNu      (void);
  void ScanMeta              (void);

  // Private data members
  //
  double         fMaxEv;          ///< maximum energy
  bool           fEnd;            ///< end condition reached

  std::vector<string> fNuFluxFilePatterns;  ///< (potentially wildcarded) path(s)
  string    fNuFluxBranchRequest; ///< list of requested branches "entry,numi,au"
  TChain*   fNuFluxTree;          ///< TTree // REF ONLY
  TChain*   fNuMetaTree;          ///< TTree // REF ONLY

  int       fNFiles;              ///< number of files in chain
  Long64_t  fNEntries;            ///< number of flux ntuple entries
  Long64_t  fIEntry;              ///< current flux ntuple entry
  Int_t     fIFileNumber;         ///< which file for the current entry

  Double_t  fFilePOTs;            ///< # of protons-on-target represented by all files

  double    fWeight;              ///< current neutrino weight
  double    fMaxWeight;           ///< max flux neutrino weight in input file

  long int  fNUse;                ///< how often to use same entry in a row
  long int  fIUse;                ///< current # of times an entry has been used

  double    fSumWeight;           ///< sum of weights for nus thrown so far
  long int  fNNeutrinos;          ///< number of flux neutrinos thrown so far
  long int  fNEntriesUsed;        ///< number of entries read from files
  double    fEffPOTsPerNu;        ///< what a entry is worth ...
  double    fAccumPOTs;           ///< POTs used so far

  bool      fGenWeighted;         ///< does GenerateNext() give weights?
  bool      fAlreadyUnwgt;        ///< are input files already unweighted
                                  // i.e. are all entry "wgt" values = 1
  bool      fAllFilesMeta;        ///< do all files in chain have meta data

  GSimpleNtpEntry* fCurEntry;  ///< current entry
  GSimpleNtpNuMI*  fCurNuMI;   ///< current "numi" branch extra info
  GSimpleNtpAux*   fCurAux;    ///< current "aux" branch extra info
  TLorentzVector   fP4;        ///< reconstituted p4 vector
  TLorentzVector   fX4;        ///< reconstituted position vector
  GSimpleNtpMeta*  fCurMeta;   ///< current meta data 

  GSimpleNtpEntry* fCurEntryCopy;  ///< current entry
  GSimpleNtpNuMI*  fCurNuMICopy;   ///< current "numi" branch extra info
  GSimpleNtpAux*   fCurAuxCopy;    ///< current "aux" branch extra info

};

} // flux namespace
} // genie namespace

#endif // _SIMPLE_NTP_FLUX_H_
