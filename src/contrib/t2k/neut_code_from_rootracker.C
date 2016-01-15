//________________________________________________________________________________________
/*!

\macro   neut_code_from_rootracker.C

\brief   Macro to read a GENIE event tree in the t2k_rootracker format and calculate the 
         NEUT reaction code. Note that GENIE event files in t2k_rootracker format for 
         versions >= v2.5.1 already include the neut reaction code as a separate tree branch.

\usage   shell% root
         root[0] .L neut_code_from_rootracker.C++
         root[1] neut_code_from_rootracker("./your_rootracker_file.root");

\author  Costas Andreopoulos <costas.andreopoulos@stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

\created Nov 24, 2008

\cpright Copyright (c) 2003-2016, GENIE Neutrino MC Generator Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//_________________________________________________________________________________________

#include <cassert>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TBits.h>
#include <TObjString.h>
#include <TLorentzVector.h>

using std::cout;
using std::endl;
using std::ostringstream;
using std::ofstream;
using std::string;

void neut_code_from_rootracker(const char * filename)
{
  bool using_new_version = false; // StdHepReScat and G2NeutEvtCode branches available only for versions >= 2.5.1
  bool event_printout    = false; 

  //
  // constants
  //

  // set a max expected number of particles per event
  const int kNPmax = 100;

  // status codes
//const int kIStUndefined                  = -1;
//const int kIStInitialState               =  0;
  const int kIStStableFinalState           =  1;
//const int kIStIntermediateState          =  2;
  const int kIStDecayedState               =  3;
//const int kIStCorrelatedNucleon          = 10;
//const int kIStNucleonTarget              = 11;
//const int kIStDISPreFragmHadronicState   = 12;
//const int kIStPreDecayResonantState      = 13;
  const int kIStHadronInTheNucleus         = 14;

  // define some particle codes needed for converting GENIE -> NEUT reaction codes
  const int kPdgNuE          =   12;
  const int kPdgAntiNuE      =  -12;
  const int kPdgNuMu         =   14;
  const int kPdgAntiNuMu     =  -14;
  const int kPdgNuTau        =   16;
  const int kPdgAntiNuTau    =  -16;
  const int kPdgGamma        =    22; // photon
  const int kPdgProton       =  2212;
//const int kPdgAntiProton   = -2212;
  const int kPdgNeutron      =  2112;
//const int kPdgAntiNeutron  = -2112;
  const int kPdgPiP          =   211; // pi+
  const int kPdgPiM          =  -211; // pi-
  const int kPdgPi0          =   111; // pi0
  const int kPdgEta          =   221; // eta
  const int kPdgKP           =   321; // K+
  const int kPdgKM           =  -321; // K-
  const int kPdgK0           =   311; // K0
  const int kPdgAntiK0       =  -311; // \bar{K0}
  const int kPdgLambda       =  3122; // Lambda
  const int kPdgAntiLambda   = -3122; // \bar{Lambda}

  // stdhep momentum array indices
  const int kStdHepIdxPx = 0;
  const int kStdHepIdxPy = 1;
  const int kStdHepIdxPz = 2;
  const int kStdHepIdxE  = 3;

  //
  // get input tree
  //

  TFile file(filename, "READ");
  TTree * tree = (TTree *) file.Get("gRooTracker");
  assert(tree);

  //
  // event info in rootracker files
  //

  TBits*      EvtFlags = 0;             // generator-specific event flags
  TObjString* EvtCode = 0;              // generator-specific string with 'event code'
  int         EvtNum;                   // event num.
  double      EvtXSec;                  // cross section for selected event (1E-38 cm2)
  double      EvtDXSec;                 // cross section for selected event kinematics (1E-38 cm2 /{K^n})
  double      EvtWght;                  // weight for that event
  double      EvtProb;                  // probability for that event (given cross section, path lengths, etc)
  double      EvtVtx[4];                // event vertex position in detector coord syst (in geom units)
  int         StdHepN;                  // number of particles in particle array 
  int         StdHepPdg   [kNPmax];     // stdhep-like particle array: pdg codes (& generator specific codes for pseudoparticles)
  int         StdHepStatus[kNPmax];     // stdhep-like particle array: generator-specific status code
  int         StdHepRescat[kNPmax];     // stdhep-like particle array: intranuclear rescattering code [ >= v2.5.1 ]
  double      StdHepX4    [kNPmax][4];  // stdhep-like particle array: 4-x (x, y, z, t) of particle in hit nucleus frame (fm)
  double      StdHepP4    [kNPmax][4];  // stdhep-like particle array: 4-p (px,py,pz,E) of particle in LAB frame (GeV)
  double      StdHepPolz  [kNPmax][3];  // stdhep-like particle array: polarization vector
  int         StdHepFd    [kNPmax];     // stdhep-like particle array: first daughter
  int         StdHepLd    [kNPmax];     // stdhep-like particle array: last  daughter 
  int         StdHepFm    [kNPmax];     // stdhep-like particle array: first mother
  int         StdHepLm    [kNPmax];     // stdhep-like particle array: last  mother
  int         G2NeutEvtCode;            // NEUT code for the current GENIE event [ >= v2.5.1 ]
  int         NuParentPdg;              // parent hadron pdg code
  int         NuParentDecMode;          // parent hadron decay mode
  double      NuParentDecP4 [4];        // parent hadron 4-momentum at decay
  double      NuParentDecX4 [4];        // parent hadron 4-position at decay
  double      NuParentProP4 [4];        // parent hadron 4-momentum at production
  double      NuParentProX4 [4];        // parent hadron 4-position at production
  int         NuParentProNVtx;          // parent hadron vtx id

  //
  // get tre branches and set branch addresses
  //

  TBranch * brEvtFlags        = tree -> GetBranch ("EvtFlags");
  TBranch * brEvtCode         = tree -> GetBranch ("EvtCode");
  TBranch * brEvtNum          = tree -> GetBranch ("EvtNum");
  TBranch * brEvtXSec         = tree -> GetBranch ("EvtXSec");
  TBranch * brEvtDXSec        = tree -> GetBranch ("EvtDXSec");
  TBranch * brEvtWght         = tree -> GetBranch ("EvtWght");
  TBranch * brEvtProb         = tree -> GetBranch ("EvtProb");
  TBranch * brEvtVtx          = tree -> GetBranch ("EvtVtx");
  TBranch * brStdHepN         = tree -> GetBranch ("StdHepN");
  TBranch * brStdHepPdg       = tree -> GetBranch ("StdHepPdg");
  TBranch * brStdHepStatus    = tree -> GetBranch ("StdHepStatus");
  TBranch * brStdHepRescat    = (using_new_version) ? tree -> GetBranch ("StdHepRescat") : 0;
  TBranch * brStdHepX4        = tree -> GetBranch ("StdHepX4");
  TBranch * brStdHepP4        = tree -> GetBranch ("StdHepP4");
  TBranch * brStdHepPolz      = tree -> GetBranch ("StdHepPolz");
  TBranch * brStdHepFd        = tree -> GetBranch ("StdHepFd");
  TBranch * brStdHepLd        = tree -> GetBranch ("StdHepLd");
  TBranch * brStdHepFm        = tree -> GetBranch ("StdHepFm");
  TBranch * brStdHepLm        = tree -> GetBranch ("StdHepLm");
  TBranch * brG2NeutEvtCode   = (using_new_version) ? tree -> GetBranch ("G2NeutEvtCode") : 0;
  TBranch * brNuParentPdg     = tree -> GetBranch ("NuParentPdg");
  TBranch * brNuParentDecMode = tree -> GetBranch ("NuParentDecMode");
  TBranch * brNuParentDecP4   = tree -> GetBranch ("NuParentDecP4");
  TBranch * brNuParentDecX4   = tree -> GetBranch ("NuParentDecX4");
  TBranch * brNuParentProP4   = tree -> GetBranch ("NuParentProP4");     
  TBranch * brNuParentProX4   = tree -> GetBranch ("NuParentProX4");     
  TBranch * brNuParentProNVtx = tree -> GetBranch ("NuParentProNVtx");   

  brEvtFlags        -> SetAddress ( &EvtFlags        );
  brEvtCode         -> SetAddress ( &EvtCode         );
  brEvtNum          -> SetAddress ( &EvtNum          );
  brEvtXSec         -> SetAddress ( &EvtXSec         );
  brEvtDXSec        -> SetAddress ( &EvtDXSec        );
  brEvtWght         -> SetAddress ( &EvtWght         );
  brEvtProb         -> SetAddress ( &EvtProb         );
  brEvtVtx          -> SetAddress (  EvtVtx          );
  brStdHepN         -> SetAddress ( &StdHepN         );
  brStdHepPdg       -> SetAddress (  StdHepPdg       );
  brStdHepStatus    -> SetAddress (  StdHepStatus    );
  if(using_new_version) {
  brStdHepRescat    -> SetAddress (  StdHepRescat    ); 
  }
  brStdHepX4        -> SetAddress (  StdHepX4        );
  brStdHepP4        -> SetAddress (  StdHepP4        );
  brStdHepPolz      -> SetAddress (  StdHepPolz      );
  brStdHepFd        -> SetAddress (  StdHepFd        );
  brStdHepLd        -> SetAddress (  StdHepLd        );
  brStdHepFm        -> SetAddress (  StdHepFm        );
  brStdHepLm        -> SetAddress (  StdHepLm        );
  if(using_new_version) {
  brG2NeutEvtCode   -> SetAddress ( &G2NeutEvtCode   );
  }
  brNuParentPdg     -> SetAddress ( &NuParentPdg     );
  brNuParentDecMode -> SetAddress ( &NuParentDecMode );
  brNuParentDecP4   -> SetAddress (  NuParentDecP4   );
  brNuParentDecX4   -> SetAddress (  NuParentDecX4   );
  brNuParentProP4   -> SetAddress (  NuParentProP4   );     
  brNuParentProX4   -> SetAddress (  NuParentProX4   );     
  brNuParentProNVtx -> SetAddress ( &NuParentProNVtx );   


  //
  // open a text file to save the reaction codes
  //
  ostringstream outfilename;
  outfilename << filename << ".reaction_codes";

  ofstream outfile(outfilename.str().c_str());
  outfile << "#" << endl;
  outfile << "# NEUT reaction code for GENIE file: " << filename << endl;
  outfile << "#" << endl;
  outfile << "#" << endl;
  outfile << "# GENIE evt nu.	NEUT code" << endl;

  //
  // event loop
  //

  int n = tree->GetEntries(); 
  printf("Number of entries: %d", n);

  for(int i=0; i < tree->GetEntries(); i++) {

    //
    // get next event
    //

    printf("\n\n ** Current entry: %d \n", i);
    tree->GetEntry(i);

    //
    // print event
    //
    if(event_printout) {
       printf("\n -----------------------------------------------------------------------------------------------------------------");
       printf("\n Event code                 : %s", EvtCode->String().Data());
       printf("\n Event x-section            : %10.5f * 1E-38* cm^2",  EvtXSec);
       printf("\n Event kinematics x-section : %10.5f * 1E-38 * cm^2/{K^n}", EvtDXSec);
       printf("\n Event weight               : %10.8f", EvtWght);
       printf("\n Event vertex               : x = %8.2f mm, y = %8.2f mm, z = %8.2f mm", EvtVtx[0], EvtVtx[1], EvtVtx[2]);
       printf("\n * Particle list:");
       printf("\n --------------------------------------------------------------------------------------------------------------------------");
       printf("\n | Idx | Ist |    PDG     | Rescat |   Mother  |  Daughter |   Px   |   Py   |   Pz   |   E    |   x    |   y    |    z   |");
       printf("\n |     |     |            |        |           |           |(GeV/c) |(GeV/c) |(GeV/c) | (GeV)  |  (fm)  |  (fm)  |   (fm) |");
       printf("\n --------------------------------------------------------------------------------------------------------------------------");

       for(int ip=0; ip<StdHepN; ip++) {
          printf("\n | %3d | %3d | %10d | %6d | %3d | %3d | %3d | %3d | %6.3f | %6.3f | %6.3f | %6.3f | %6.3f | %6.3f | %6.3f |",
             ip, StdHepStatus[ip],  StdHepPdg[ip], StdHepRescat[ip], 
             StdHepFm[ip],  StdHepLm[ip], StdHepFd[ip],  StdHepLd[ip],
             StdHepP4[ip][0], StdHepP4[ip][1], StdHepP4[ip][2], StdHepP4[ip][3],
             StdHepX4[ip][0], StdHepX4[ip][1], StdHepX4[ip][2]); 
       }
       printf("\n --------------------------------------------------------------------------------------------------------------------------");
       printf("\n * Flux Info:");
       printf("\n Parent hadron pdg code    : %d", NuParentPdg);
       printf("\n Parent hadron decay mode  : %d", NuParentDecMode);
       printf("\n Parent hadron 4p at decay : Px = %6.3f GeV/c, Py = %6.3f GeV/c, Pz = %6.3f GeV/c, E = %6.3f GeV", 
                  NuParentDecP4[0], NuParentDecP4[1], NuParentDecP4[2], NuParentDecP4[3]);
       printf("\n Parent hadron 4p at prod. : Px = %6.3f GeV/c, Py = %6.3f GeV/c, Pz = %6.3f GeV/c, E = %6.3f GeV", 
                  NuParentProP4[0], NuParentProP4[1], NuParentProP4[2], NuParentProP4[3]);
       printf("\n Parent hadron 4x at decay : x = %6.3f m, y = %6.3f m, z = %6.3f m, t = %6.3f s", 
                  NuParentDecX4[0], NuParentDecX4[1], NuParentDecX4[2], NuParentDecX4[3]);
       printf("\n Parent hadron 4x at prod. : x = %6.3f m, y = %6.3f m, z = %6.3f m, t = %6.3f s", 
                  NuParentProX4[0], NuParentProX4[1], NuParentProX4[2], NuParentProX4[3]);
       printf("\n -------------------------------------------------------------------------------------------------------------------------- \n");
    }

    //
    // figure out the NEUT code for the current GENIE event
    //

    int evtype = 0;

    // basic GENIE reaction types
    string genie_code = EvtCode->GetString().Data();
    bool is_cc    = (genie_code.find("Weak[CC]") != string::npos);
    bool is_nc    = (genie_code.find("Weak[NC]") != string::npos);
    bool is_charm = (genie_code.find("charm")    != string::npos);
    bool is_qel   = (genie_code.find("QES")      != string::npos);
    bool is_dis   = (genie_code.find("DIS")      != string::npos);
    bool is_res   = (genie_code.find("RES")      != string::npos);
    bool is_cohpi = (genie_code.find("COH")      != string::npos);
    bool is_ve    = (genie_code.find("NuEEL")    != string::npos);
    bool is_imd   = (genie_code.find("IMD")      != string::npos);

    // basic initial state info
    int inu_pos    =  0;
    int fsl_pos    =  StdHepFd[inu_pos]; // final state primary lepton: neutrino daughter
    assert(fsl_pos>0);
    int tgt_pos    =  1;
    int hitnuc_pos = -1; // can be at slot 2 (for nuclear targets), slot 1 (for free nuc targets) or undefined (for ve-, IMD, coherent etc)

    int nu_code = StdHepPdg[inu_pos];
    bool is_nu    = (nu_code == kPdgNuE     || nu_code == kPdgNuMu     || nu_code == kPdgNuTau    );
    bool is_nubar = (nu_code == kPdgAntiNuE || nu_code == kPdgAntiNuMu || nu_code == kPdgAntiNuTau);

    int target_code = StdHepPdg[tgt_pos];
    bool nuclear_target = (target_code > 1000000000);

    bool hitnuc_set  = false;
    bool is_p        = false; 
    bool is_n        = false; 
    for(int ip = 1; ip <= 2; ip++) 
    {
       int ghep_pdgc = StdHepPdg[ip];
       if(ghep_pdgc == kPdgProton) {
         hitnuc_pos = ip;
         hitnuc_set = true;
         is_p       = true;
         break;
       }
       if(ghep_pdgc == kPdgNeutron) {
         hitnuc_pos = ip;
         hitnuc_set = true;
         is_n       = true;
         break;
       }
    }

    // calculate the true (at the hit nucleon rest frame) hadronic invariant mass
    // calculate only if the hit nucleon is set (calc not needed for coherent events etc)

    bool W_gt_2 = false;
    if(hitnuc_pos != -1) 
    {
      TLorentzVector p4v ( StdHepP4 [inu_pos]    [kStdHepIdxPx],
                           StdHepP4 [inu_pos]    [kStdHepIdxPy],
                           StdHepP4 [inu_pos]    [kStdHepIdxPz],
                           StdHepP4 [inu_pos]    [kStdHepIdxE]   );
      TLorentzVector p4l ( StdHepP4 [fsl_pos]    [kStdHepIdxPx],
                           StdHepP4 [fsl_pos]    [kStdHepIdxPy],
                           StdHepP4 [fsl_pos]    [kStdHepIdxPz],
                           StdHepP4 [fsl_pos]    [kStdHepIdxE]   );
      TLorentzVector p4n ( StdHepP4 [hitnuc_pos] [kStdHepIdxPx], 
                           StdHepP4 [hitnuc_pos] [kStdHepIdxPy],
                           StdHepP4 [hitnuc_pos] [kStdHepIdxPz],
                           StdHepP4 [hitnuc_pos] [kStdHepIdxE]   );

      TLorentzVector q = p4v - p4l;
      TLorentzVector w = p4n + q;

      double W = w.Mag();
      W_gt_2 = (W > 2.0);
    }
          
    // (quasi-)elastic, nc+cc, nu+nubar
    //
    if      (is_qel && !is_charm && is_cc && is_nu           ) evtype =   1;
    else if (is_qel && !is_charm && is_nc && is_nu && is_p   ) evtype =  51;   
    else if (is_qel && !is_charm && is_nc && is_nu && is_n   ) evtype =  52;
    else if (is_qel && !is_charm && is_cc && is_nubar        ) evtype =  -1;
    else if (is_qel && !is_charm && is_nc && is_nubar && is_p) evtype = -51;
    else if (is_qel && !is_charm && is_nc && is_nubar && is_n) evtype = -52;
           
    // quasi-elastic charm production
    //
    else if (is_qel && is_charm && is_cc && is_nu    ) evtype =   25;
    else if (is_qel && is_charm && is_cc && is_nubar ) evtype =  -25;

    // inverse mu- (tau-) decay and ve- elastic
    //
    else if ( is_imd ) evtype =  9;
    else if ( is_ve  ) evtype = 59;
               
    // coherent pi, nc+cc, nu+nubar
    //
    else if (is_cohpi && is_cc && is_nu   ) evtype =  16;
    else if (is_cohpi && is_cc && is_nubar) evtype = -16;
    else if (is_cohpi && is_nc && is_nu   ) evtype =  36;
    else if (is_cohpi && is_nc && is_nubar) evtype = -36;
                     
    // dis, W>2, nc+cc, nu+nubar
    // (charm DIS not simulated by NEUT, will bundle GENIE charm DIS into this category)
    //
    else if (is_dis && W_gt_2 && is_cc && is_nu   ) evtype =  26;
    else if (is_dis && W_gt_2 && is_nc && is_nu   ) evtype =  46;
    else if (is_dis && W_gt_2 && is_cc && is_nubar) evtype = -26; 
    else if (is_dis && W_gt_2 && is_nc && is_nubar) evtype = -46; 

    // resonance or dis with W < 2 GeV
    //
    else if ( is_res || (is_dis && !W_gt_2) ) {
        
       //cout << "Current event is RES or DIS with W<2" << endl;
        
       // check the number of pions and nucleons in the primary hadronic system
       // (_before_ intranuclear rescattering)
       //
       int nn=0, np=0, npi0=0, npip=0, npim=0, nKp=0, nKm=0, nK0=0, neta=0, nlambda=0, ngamma=0;

       for(int ip = 0; ip < StdHepN; ip++) 
       {
           int ghep_ist     = StdHepStatus[ip];
           int ghep_pdgc    = StdHepPdg[ip];
           int ghep_fm      = StdHepFm[ip];
           int ghep_fmpdgc  = (ghep_fm==-1) ? 0 : StdHepPdg[ghep_fm];
        
           // For nuclear targets use hadrons marked as 'hadron in the nucleus'
           // which are the ones passed in the intranuclear rescattering
           // For free nucleon targets use particles marked as 'final state'
           // but make an exception for decayed pi0's,eta's (count them and not their daughters)

           bool decayed         = (ghep_ist==kIStDecayedState && (ghep_pdgc==kPdgPi0 || ghep_pdgc==kPdgEta));
           bool parent_included = (ghep_fmpdgc==kPdgPi0 || ghep_fmpdgc==kPdgEta);

           bool count_it =
               ( nuclear_target && ghep_ist==kIStHadronInTheNucleus) ||
               (!nuclear_target && decayed) ||
               (!nuclear_target && ghep_ist==kIStStableFinalState && !parent_included);

           if(!count_it) continue;
                
           if(ghep_pdgc == kPdgProton )    np++;            // p
           if(ghep_pdgc == kPdgNeutron)    nn++;            // n
           if(ghep_pdgc == kPdgPiP)        npip++;          // pi+
           if(ghep_pdgc == kPdgPiM)        npim++;          // pi-
           if(ghep_pdgc == kPdgPi0)        npi0++;          // pi0
           if(ghep_pdgc == kPdgEta)        neta++;          // eta0
           if(ghep_pdgc == kPdgKP)         nKp++;           // K+
           if(ghep_pdgc == kPdgKM)         nKm++;           // K-
           if(ghep_pdgc == kPdgK0)         nK0++;           // K0
           if(ghep_pdgc == kPdgAntiK0)     nK0++;           // K0
           if(ghep_pdgc == kPdgLambda)     nlambda++;       // Lamda
           if(ghep_pdgc == kPdgAntiLambda) nlambda++;       // Lamda
           if(ghep_pdgc == kPdgGamma)      ngamma++;        // photon
       }
       if(event_printout) {
         cout  << "Num of primary particles: \n p = " << np << ", n = " << nn
               << ", pi+ = " << npip << ", pi- = " << npim << ", pi0 = " << npi0 
               << ", eta = " << neta 
               << ", K+ = " << nKp << ", K- = " << nKm << ", K0 = " << nK0 
               << ", Lambda's = " << nlambda
               << ", gamma's = " << ngamma
               << endl;
       }            
       int nnuc = np + nn;
       int npi  = npi0 + npip + npim;
       int nK   = nK0 + nKp + nKm;
       int neKL = neta + nK + nlambda;
              
       bool is_radiative_dec = (nnuc==1) && (npi==0) && (ngamma==1);

       //
       // single gamma from resonances
       //
              
       if      (is_res && is_nu    && is_cc && is_n && is_radiative_dec) evtype =  17;
       else if (is_res && is_nu    && is_nc && is_n && is_radiative_dec) evtype =  38;
       else if (is_res && is_nu    && is_nc && is_p && is_radiative_dec) evtype =  39;
               
       else if (is_res && is_nubar && is_cc && is_p && is_radiative_dec) evtype = -17;
       else if (is_res && is_nubar && is_nc && is_n && is_radiative_dec) evtype = -38;
       else if (is_res && is_nubar && is_nc && is_p && is_radiative_dec) evtype = -39;
               
       //
       // single pi (res + non-res bkg)
       //
            
       // nu CC
       else if (is_nu    && is_cc && is_p && np==1 && nn==0 && npip==1 && npim==0 && npi0==0 && neKL==0) evtype =  11;
       else if (is_nu    && is_cc && is_n && np==1 && nn==0 && npip==0 && npim==0 && npi0==1 && neKL==0) evtype =  12;
       else if (is_nu    && is_cc && is_n && np==0 && nn==1 && npip==1 && npim==0 && npi0==0 && neKL==0) evtype =  13;
           
       // nu NC
       else if (is_nu    && is_nc && is_n && np==0 && nn==1 && npip==0 && npim==0 && npi0==1 && neKL==0) evtype =  31;
       else if (is_nu    && is_nc && is_p && np==1 && nn==0 && npip==0 && npim==0 && npi0==1 && neKL==0) evtype =  32;
       else if (is_nu    && is_nc && is_n && np==1 && nn==0 && npip==0 && npim==1 && npi0==0 && neKL==0) evtype =  33;
       else if (is_nu    && is_nc && is_p && np==0 && nn==1 && npip==1 && npim==0 && npi0==0 && neKL==0) evtype =  34;
           
       //nubar CC
       else if (is_nubar && is_cc && is_n && np==0 && nn==1 && npip==0 && npim==1 && npi0==0 && neKL==0) evtype = -11;
       else if (is_nubar && is_cc && is_p && np==0 && nn==1 && npip==0 && npim==0 && npi0==1 && neKL==0) evtype = -12;
       else if (is_nubar && is_cc && is_p && np==1 && nn==0 && npip==0 && npim==1 && npi0==0 && neKL==0) evtype = -13;
                     
       //nubar NC
       else if (is_nubar && is_nc && is_n && np==0 && nn==1 && npip==0 && npim==0 && npi0==1 && neKL==0) evtype = -31;
       else if (is_nubar && is_nc && is_p && np==1 && nn==0 && npip==0 && npim==0 && npi0==1 && neKL==0) evtype = -32;
       else if (is_nubar && is_nc && is_n && np==1 && nn==0 && npip==0 && npim==1 && npi0==0 && neKL==0) evtype = -33;
       else if (is_nubar && is_nc && is_p && np==0 && nn==1 && npip==1 && npim==0 && npi0==0 && neKL==0) evtype = -34;
              
       //
       // single eta from res
       //
              
       else if (is_res &&  is_nu    && is_cc && is_n && np==1 && nn==0 && npi==0 && nK==0 && nlambda==0 && neta==1) evtype =  22;
       else if (is_res &&  is_nu    && is_nc && is_n && np==0 && nn==1 && npi==0 && nK==0 && nlambda==0 && neta==1) evtype =  42;
       else if (is_res &&  is_nu    && is_nc && is_p && np==1 && nn==0 && npi==0 && nK==0 && nlambda==0 && neta==1) evtype =  43;
              
       else if (is_res &&  is_nubar && is_cc && is_p && np==0 && nn==1 && npi==0 && nK==0 && nlambda==0 && neta==1) evtype = -22;
       else if (is_res &&  is_nubar && is_nc && is_n && np==0 && nn==1 && npi==0 && nK==0 && nlambda==0 && neta==1) evtype = -42;
       else if (is_res &&  is_nubar && is_nc && is_p && np==1 && nn==0 && npi==0 && nK==0 && nlambda==0 && neta==1) evtype = -43;
              
       //
       // single K from res
       //
              
       else if (is_res &&  is_nu    && is_cc && is_n && nnuc==0 && npi==0 && nK==1 && nlambda==1 && neta==0) evtype =  23;
       else if (is_res &&  is_nu    && is_nc && is_n && nnuc==0 && npi==0 && nK==1 && nlambda==1 && neta==0) evtype =  44;
       else if (is_res &&  is_nu    && is_nc && is_p && nnuc==0 && npi==0 && nK==1 && nlambda==1 && neta==0) evtype =  45;
              
       else if (is_res &&  is_nubar && is_cc && is_p && nnuc==0 && npi==0 && nK==1 && nlambda==1 && neta==0) evtype = -23;
       else if (is_res &&  is_nubar && is_nc && is_n && nnuc==0 && npi==0 && nK==1 && nlambda==1 && neta==0) evtype = -44;
       else if (is_res &&  is_nubar && is_nc && is_p && nnuc==0 && npi==0 && nK==1 && nlambda==1 && neta==0) evtype = -45;

       //
       // multi-pi (res or dis), W<2GeV
       //
              
       else if (is_nu    && is_cc && npi>1) evtype =  21;
       else if (is_nu    && is_nc && npi>1) evtype =  41;
       else if (is_nubar && is_cc && npi>1) evtype = -21;
       else if (is_nubar && is_nc && npi>1) evtype = -41;

       //
       // rare final state for RES or low-W (<2GeV) DIS events
       // (eg K0\bar{K0} final states, N0(1720) -> Sigma- K+ res decays, etc)
       // bundled-in with multi-pi
       //
       else {              
           
          cout << "Rare RES/low-W DIS final state: Bundled-in with multi-pi events" << endl;

               if (is_nu    && is_cc) evtype =  21;
          else if (is_nu    && is_nc) evtype =  41;
          else if (is_nubar && is_cc) evtype = -21;
          else if (is_nubar && is_nc) evtype = -41;
       }
    }

    cout << " *** GENIE event = " << i << " --> NEUT reaction code = " << evtype << endl;
    if(using_new_version) {
       // for validation, use a file generated + converted with the CVS head version of GENIE
       // where the NeutCode is stored at the t2k_rootracker tree
       cout << "NEUT reaction code stored at the rootracker file = " << G2NeutEvtCode << endl;
       assert(evtype == G2NeutEvtCode);
    }

    // save at the output file
    outfile << "\t" << i << "\t" << evtype << endl;

 } // event loop

 outfile.close();

 file.Close();
}
