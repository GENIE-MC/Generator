//________________________________________________________________________________________
/*!

\macro   read_t2k_rootracker.C

\brief   Test macro to read a GENIE event tree in the t2k_rootracker format
         and print-out the neutrino event and flux pass-through info.

\usage   shell% root
         root[0] .L read_t2k_rootracker.C++
         root[1] read_t2k_rootracker("./your_rootracker_file.root");

\author  Costas Andreopoulos <costas.andreopoulos@stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

\created Nov 24, 2008

\cpright Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//_________________________________________________________________________________________

#include <iostream>

#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TBits.h>
#include <TObjString.h>

using std::cout;
using std::endl;

void read_t2k_rootracker(const char * filename)
{
  bool using_new_version = false; // StdHepReScat and G2NeutEvtCode branches available only for versions >= 2.5.1

  const int kNPmax = 100;

  TFile file(filename, "READ");
  TTree * tree = (TTree *) file.Get("gRooTracker");
  assert(tree);

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

  // get branches
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

  // set address
  brEvtFlags        -> SetAddress ( &EvtFlags         );
  brEvtCode         -> SetAddress ( &EvtCode          );
  brEvtNum          -> SetAddress ( &EvtNum           );
  brEvtXSec         -> SetAddress ( &EvtXSec          );
  brEvtDXSec        -> SetAddress ( &EvtDXSec         );
  brEvtWght         -> SetAddress ( &EvtWght          );
  brEvtProb         -> SetAddress ( &EvtProb          );
  brEvtVtx          -> SetAddress (  EvtVtx           );
  brStdHepN         -> SetAddress ( &StdHepN          );
  brStdHepPdg       -> SetAddress (  StdHepPdg        );
  brStdHepStatus    -> SetAddress (  StdHepStatus     );
  if(using_new_version) {
  brStdHepRescat    -> SetAddress (  StdHepRescat     );
  }
  brStdHepX4        -> SetAddress (  StdHepX4         );
  brStdHepP4        -> SetAddress (  StdHepP4         );
  brStdHepPolz      -> SetAddress (  StdHepPolz       );
  brStdHepFd        -> SetAddress (  StdHepFd         );
  brStdHepLd        -> SetAddress (  StdHepLd         );
  brStdHepFm        -> SetAddress (  StdHepFm         );
  brStdHepLm        -> SetAddress (  StdHepLm         );
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

  int n = tree->GetEntries(); 
  printf("Number of entries: %d", n);

  // read tree & print some info for each entry

  for(int i=0; i < tree->GetEntries(); i++) {
    printf("\n\n ** Current entry: %d \n", i);
    tree->GetEntry(i);

    printf("\n -----------------------------------------------------------------------------------------------------------------");
    printf("\n Event code                 : %s", EvtCode->String().Data());
    printf("\n Event x-section            : %10.5f * 1E-38* cm^2",  EvtXSec);
    printf("\n Event kinematics x-section : %10.5f * 1E-38 * cm^2/{K^n}", EvtDXSec);
    printf("\n Event weight               : %10.8f", EvtWght);
    printf("\n Event vertex               : x = %8.2f mm, y = %8.2f mm, z = %8.2f mm", EvtVtx[0], EvtVtx[1], EvtVtx[2]);
    printf("\n *Particle list:");
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
    printf("\n *Flux Info:");
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

 } // tree entries

 file.Close();
}
