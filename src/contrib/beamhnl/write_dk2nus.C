//____________________________________________________________________________
/*!
  
\brief   Macro to flatten a dk2nu flux file into a format ready for use with the BeamHNL module

\author  John Plows <komninos-john.plows \at physics.ox.ac.uk>
         University of Oxford

\cpright Copyright (c) 2003-2025, The GENIE Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
 */
//____________________________________________________________________________

#include "write_dk2nus.h"

// ---- main method ----

int write_dk2nus(bool grid               = false,
                 bool debug              = false,
                 int job                 = 0,
                 std::string indir       = "",
		 std::string inString    = "",
                 std::string outdir      = "",
                 std::string tag         = "",
                 std::string playlist    = "minervame",   // these aren't
                 std::string beamconfig  = "me000z200i")  // used in ME
{

  //gROOT->ProcessLine("#include <vector>");

  //============================================================================
  // Print options
  //============================================================================
    std::cout << "  tag: " << tag << std::endl;
    std::cout << "  User is " << USER << std::endl;
  
  //============================================================================
  // Set I/O directories
  //============================================================================
    if (grid){
      indir  = INDIR_GRID;
      outdir = OUTDIR_GRID;
    }
    else{
      indir  = debug ? INDIR_DEBUG : indir;
      if (indir == ""){
        std::cout << "\nERROR: (extract_flux.C interactive mode) " << 
                     "Must specify input directory OR debug mode.\n\n";  
        exit(EXIT_FAILURE);
      }
      outdir = (outdir == "") ? string(std::getenv("HOME")) : outdir;
    }
    indir = gSystem->ExpandPathName(indir.c_str());
    outdir = gSystem->ExpandPathName(outdir.c_str());

    std::cout << "  indir: " << indir << std::endl;
    std::cout << "  outdir: " << outdir << std::endl;

  //============================================================================
  // Make output file
  //============================================================================
    std::string out_name = tag.empty() ? 
                              Form("dk2nu_flat_%d.root", job) : 
                                Form("dk2nu_flat_%s_%d.root", tag.c_str(), job);
    if(debug) out_name = "dtest1.root";

    out_name = outdir + "/" + out_name;
   
    TFile *fout = new TFile( out_name.c_str(), "RECREATE"); 
      
  //============================================================================
  // Get input and make chain
  //============================================================================
    //choose number of files by exension
    std::string extension = ".root";
    //std::string extension = debug ? "_999_0001.root" : "_0001.root";

    //std::string newIndir = indir + G4PREFIX + playlist + "_" + beamconfig + "_";
    //std::cout << "\n  Seaching for Ntuples in " << indir + "*" + extension << std::endl;
    //std::cout << "\n  Searching for Ntuples in " << newIndir + "*" + extension << std::endl;

    std::cout << "\t extract_flux.C sees this inString:\n\t " << inString.c_str() << std::endl;
    //now build an array with each entry from the inString
    std::string arrString[m_maxNFiles]; int ifile = 0;
    size_t pos = 0;
    while(true){
      pos = inString.find('|');
      std::string tmpString = pos != std::string::npos ? inString.substr(0,pos) : inString;
      arrString[ifile] = tmpString; ifile++;
      inString = inString.substr(pos+1,inString.length()-(pos+1));
      if(pos == std::string::npos){break;}
    }

    //Add to chain
    TChain* cflux = new TChain("dk2nuTree");
    TChain* dflux = new TChain("dkmetaTree");

    for(Int_t itmpfile = 0; itmpfile < ifile; itmpfile++){
      cflux->Add( arrString[itmpfile].c_str() );
      dflux->Add( arrString[itmpfile].c_str() );
      std::cout << "Adding file " << arrString[itmpfile].c_str() << std::endl;
    }

    //std::cout << indir.c_str() << "*" << extension.c_str() << std::endl;
    //cflux->Add( (indir + "*" + extension).c_str() );
    //cflux->Add( (newIndir + "*" + extension).c_str() );
  
    //for (int k = 686; k<691; ++k) cflux->Add( Form("%s/*_%d_*.root", findir.c_str(), k) );
    std::cout << "\n  Added files to TChain" << std::endl;
  //============================================================================

  //============================================================================
  // Initialize Trees
  //============================================================================
    dkMeta = new TTree("dkMeta", "meta");
    dkTree = new TTree("dkTree", "dk2nu");
    InitialiseMetaBranches(dkMeta);
    InitialiseTreeBranches(dkTree);

  //============================================================================
  // Loop entries and write trees
  //============================================================================
    LoopEntries (cflux, dflux, grid, debug);

    fout->Write();

    fout->Close();
    
    std::cout << "\n dk2nu written at: " << out_name <<"\n"<< std::endl;
  
    return 0;
}

//==============================================================================
// Loop Entries
//==============================================================================
void LoopEntries(TChain* cflux, TChain* dflux, bool grid, bool debug)
{
  std::cout << "\n  Beginning to loop over ntuples and write flat trees" << std::endl;

  //============================================================================
  // Load dk2nu ntuple trees
  //============================================================================
    std::cout << "Setting dk2nu branch address" << std::endl;
    bsim::Dk2Nu*  dk2nu  = new bsim::Dk2Nu;
    cflux->SetBranchAddress("dk2nu", &dk2nu);
    std::cout << "Setting dk2numeta branch address" << std::endl;
    bsim::DkMeta* dkmeta = new bsim::DkMeta;
    dflux->SetBranchAddress("dkmeta", &dkmeta);

  //============================================================================
  // The actual loop
  //============================================================================
    
    Long64_t nentries = debug ? 2500 : cflux->GetEntries();
    std::cout << "\n  Looping over " << nentries << " entries." << std::endl;
    
    for (Long64_t i=0; i < nentries; ++i ) {
      if (i % (nentries/100) == 0){
	std::cout << "entry " << i << " (" << i/(nentries/100) << "% done)" << std::endl;
      }
      cflux->GetEntry(i);
      dflux->GetEntry(i);
      // Sometimes g4numi ntuples are corrupt. Useful for debugging.
      if (i % (nentries/20) == 0){
	std::cout << (cflux->GetFile())->GetName() << std::endl;
      }

      // Fill the branches
      if( i==0 ) FillMetaBranches(dkmeta); //only 1 meta entry needed.
      FillTreeBranches(dk2nu);

    } // end entry loop
  
    std::cout << "\n Done writing trees" << std::endl;
}

//============================================================================
// Initialise branches here
//============================================================================
void InitialiseMetaBranches(TTree * meta)
{
  assert( meta );
  std::cout << "  Initialising meta branches\n";

  mJob = 0;
  mPots = 0.0;
  mBeam0x = 0.0;
  mBeam0y = 0.0;
  mBeam0z = 0.0;
  mBeamhwidth = 0.0;
  mBeamvwidth = 0.0;
  mBeamdxdz = 0.0;
  mBeamdydz = 0.0;
  for( unsigned int i = 0; i < maxC; i++ ){
    mBeamsim[i] = 0;
    mPhysics[i] = 0;
    mPhyscuts[i] = 0;
    mTgtcfg[i] = 0;
    mHorncfg[i] = 0;
    mDkvolcfg[i] = 0;
  }
  for( unsigned int i = 0; i < maxArray * maxC; i++ ){
    mLocationDotName[i] = 0;
  }
  for( unsigned int i = 0; i < maxArray; i++ ){
    mLocationDotX[i] = 0.0;
    mLocationDotY[i] = 0.0;
    mLocationDotZ[i] = 0.0;
  }
  //if( mLocationDotName.size() > 0 ) mLocationDotName.clear();
  //mLocationDotName = std::vector<std::string>();

  meta->Branch( "job",           &mJob,                 "job/I"                   );
  meta->Branch( "pots",          &mPots,                "pots/D"                  );
  meta->Branch( "beamsim",        mBeamsim,             "beamsim/C"               );
  meta->Branch( "physics",        mPhysics,             "physics/C"               );
  meta->Branch( "physcuts",       mPhyscuts,            "physcuts/C"              );
  meta->Branch( "tgtcfg",         mTgtcfg,              "tgtcfg/C"                );
  meta->Branch( "horncfg",        mHorncfg,             "horncfg/C"               );
  meta->Branch( "dkvolcfg",       mDkvolcfg,            "dkvolcfg/C"              );
  meta->Branch( "beam0x",        &mBeam0x,              "beam0x/D"                );
  meta->Branch( "beam0y",        &mBeam0y,              "beam0y/D"                );
  meta->Branch( "beam0z",        &mBeam0z,              "beam0z/D"                );
  meta->Branch( "beamhwidth",    &mBeamhwidth,          "beamhwidth/D"            );
  meta->Branch( "beamvwidth",    &mBeamvwidth,          "beamvwidth/D"            );
  meta->Branch( "beamdxdz",      &mBeamdxdz,            "beamdxdz/D"              );
  meta->Branch( "beamdydz",      &mBeamdydz,            "beamdydz/D"              );
  meta->Branch( "arSize",        &mArSize,              "arSize/I"                );
  meta->Branch( "location_x",     mLocationDotX,        "location_x[arSize]/D"    );
  meta->Branch( "location_y",     mLocationDotY,        "location_y[arSize]/D"    );
  meta->Branch( "location_z",     mLocationDotZ,        "location_z[arSize]/D"    );
  //meta->Branch( "location_name", &mLocationDotName                                );
  meta->Branch( "location_name",  mLocationDotName,     "location_name/C"         );
  //meta->Branch( "vintnames",      mVintnames                                      );
  //meta->Branch( "vdblnames",      mVdblnames                                      );
  
  std::cout << "  Meta branches initialised\n";
}

void InitialiseTreeBranches(TTree * tree)
{
  assert(tree);
  std::cout << "  Initialising dk2nu branches\n";

  dJob = 0;
  dPotnum = 0.0;
  dPpvx = 0.0;
  dPpvy = 0.0;
  dPpvz = 0.0;
  dDecayDotNorig = 0;
  dDecayDotNdecay = 0;
  dDecayDotNtype = 0;
  dDecayDotVx = 0.0;
  dDecayDotVy = 0.0;
  dDecayDotVz = 0.0;
  dDecayDotPdpx = 0.0;
  dDecayDotPdpy = 0.0;
  dDecayDotPdpz = 0.0;
  dDecayDotPpdxdz = 0.0;
  dDecayDotPpdydz = 0.0;
  dDecayDotPppz = 0.0;
  dDecayDotPpenergy = 0.0;
  dDecayDotPpmedium = 0;
  dDecayDotPtype = 0;
  dDecayDotMuparpx = 0.0;
  dDecayDotMuparpy = 0.0;
  dDecayDotMuparpz = 0.0;
  dDecayDotMupare = 0.0;
  dDecayDotNecm = 0.0;
  dDecayDotNimpwt = 0.0;
  dTgtexitDotTvx = 0.0;
  dTgtexitDotTvy = 0.0;
  dTgtexitDotTvz = 0.0;
  dTgtexitDotTpx = 0.0;
  dTgtexitDotTpy = 0.0;
  dTgtexitDotTpz = 0.0;
  dTgtexitDotTptype = 0;
  dTgtexitDotTgen = 0;

  for( unsigned int i = 0; i < maxArray; i++ ){
    dNurayDotPx[i] = 0.0;
    dNurayDotPy[i] = 0.0;
    dNurayDotPz[i] = 0.0;
    dNurayDotE[i] = 0.0;
    dNurayDotWgt[i] = 0.0;
    dAncestorDotPdg[i] = 0;
    dAncestorDotStartx[i] = 0.0;
    dAncestorDotStarty[i] = 0.0;
    dAncestorDotStartz[i] = 0.0;
    dAncestorDotStartt[i] = 0.0;
    dAncestorDotStartpx[i] = 0.0;
    dAncestorDotStartpy[i] = 0.0;
    dAncestorDotStartpz[i] = 0.0;
    dAncestorDotStoppx[i] = 0.0;
    dAncestorDotStoppy[i] = 0.0;
    dAncestorDotStoppz[i] = 0.0;
    dAncestorDotPolx[i] = 0.0;
    dAncestorDotPoly[i] = 0.0;
    dAncestorDotPolz[i] = 0.0;
    dAncestorDotPprodpx[i] = 0.0;
    dAncestorDotPprodpy[i] = 0.0;
    dAncestorDotPprodpz[i] = 0.0;
    dTrajDotTrkx[i] = 0.0;
    dTrajDotTrky[i] = 0.0;
    dTrajDotTrkz[i] = 0.0;
    dTrajDotTrkpx[i] = 0.0;
    dTrajDotTrkpy[i] = 0.0;
    dTrajDotTrkpz[i] = 0.0;
  }

  for( unsigned int i = 0; i < maxArray * maxC; i++ ){
    dAncestorDotProc[i] = 0;
    dAncestorDotIvol[i] = 0;
    dAncestorDotImat[i] = 0;
  }
  
  tree->Branch( "job",              &dJob,                  "job/I"                             );
  tree->Branch( "potnum",           &dPotnum,               "potnum/D"                          );
  tree->Branch( "decay_norig",      &dDecayDotNorig,        "decay_norig/I"                     );
  tree->Branch( "decay_ndecay",     &dDecayDotNdecay,       "decay_ndecay/I"                    );
  tree->Branch( "decay_ntype",      &dDecayDotNtype,        "decay_ntype/I"                     );
  tree->Branch( "decay_vx",         &dDecayDotVx,           "decay_vx/D"                        );
  tree->Branch( "decay_vy",         &dDecayDotVy,           "decay_vy/D"                        );
  tree->Branch( "decay_vz",         &dDecayDotVz,           "decay_vz/D"                        );
  tree->Branch( "decay_pdpx",       &dDecayDotPdpx,         "decay_pdpx/D"                      );
  tree->Branch( "decay_pdpy",       &dDecayDotPdpy,         "decay_pdpy/D"                      );
  tree->Branch( "decay_pdpz",       &dDecayDotPdpz,         "decay_pdpz/D"                      );
  tree->Branch( "decay_ppdxdz",     &dDecayDotPpdxdz,       "decay_ppdxdz/D"                    );
  tree->Branch( "decay_ppdydz",     &dDecayDotPpdydz,       "decay_ppdydz/D"                    );
  tree->Branch( "decay_pppz",       &dDecayDotPppz,         "decay_pppz/D"                      );
  tree->Branch( "decay_ppenergy",   &dDecayDotPpenergy,     "decay_ppenergy/D"                  );
  tree->Branch( "decay_ppmedium",   &dDecayDotPpmedium,     "decay_ppmedium/I"                  );
  tree->Branch( "decay_ptype",      &dDecayDotPtype,        "decay_ptype/I"                     );
  tree->Branch( "decay_muparpx",    &dDecayDotMuparpx,      "decay_muparpx/D"                   );
  tree->Branch( "decay_muparpy",    &dDecayDotMuparpy,      "decay_muparpy/D"                   );
  tree->Branch( "decay_muparpz",    &dDecayDotMuparpz,      "decay_muparpz/D"                   );
  tree->Branch( "decay_mupare",     &dDecayDotMupare,       "decay_muparpe/D"                   );
  tree->Branch( "decay_necm",       &dDecayDotNecm,         "decay_necm/D"                      );
  tree->Branch( "decay_nimpwt",     &dDecayDotNimpwt,       "decay_nimpwt/D"                    );
  tree->Branch( "nuray_size",       &dArSize,               "nuray_size/I"                      );
  tree->Branch( "nuray_px",          dNurayDotPx,           "nuray_px[nuray_size]/D"            );
  tree->Branch( "nuray_py",          dNurayDotPy,           "nuray_py[nuray_size]/D"            );
  tree->Branch( "nuray_pz",          dNurayDotPz,           "nuray_pz[nuray_size]/D"            );
  tree->Branch( "nuray_E",           dNurayDotE,            "nuray_E[nuray_size]/D"             );
  tree->Branch( "nuray_wgt",         dNurayDotWgt,          "nuray_wgt[nuray_size]/D"           );
  tree->Branch( "ancestor_size",    &dAnArSize,             "ancestor_size/I"                   );
  tree->Branch( "ancestor_pdg",      dAncestorDotPdg,       "ancestor_pdg[ancestor_size]/I"     );
  tree->Branch( "ancestor_startx",   dAncestorDotStartx,    "ancestor_startx[ancestor_size]/D"  );
  tree->Branch( "ancestor_starty",   dAncestorDotStarty,    "ancestor_starty[ancestor_size]/D"  );
  tree->Branch( "ancestor_startz",   dAncestorDotStartz,    "ancestor_startz[ancestor_size]/D"  );
  tree->Branch( "ancestor_startt",   dAncestorDotStartt,    "ancestor_startt[ancestor_size]/D"  );
  tree->Branch( "ancestor_startpx",  dAncestorDotStartpx,   "ancestor_startpx[ancestor_size]/D" );
  tree->Branch( "ancestor_startpy",  dAncestorDotStartpy,   "ancestor_startpy[ancestor_size]/D" );
  tree->Branch( "ancestor_startpz",  dAncestorDotStartpz,   "ancestor_startpz[ancestor_size]/D" );
  tree->Branch( "ancestor_stoppx",   dAncestorDotStoppx,    "ancestor_stoppx[ancestor_size]/D"  );
  tree->Branch( "ancestor_stoppy",   dAncestorDotStoppy,    "ancestor_stoppy[ancestor_size]/D"  );
  tree->Branch( "ancestor_stoppz",   dAncestorDotStoppz,    "ancestor_stoppz[ancestor_size]/D"  );
  tree->Branch( "ancestor_polx",     dAncestorDotPolx,      "ancestor_polx[ancestor_size]/D"    );
  tree->Branch( "ancestor_poly",     dAncestorDotPoly,      "ancestor_poly[ancestor_size]/D"    );
  tree->Branch( "ancestor_polz",     dAncestorDotPolz,      "ancestor_polz[ancestor_size]/D"    );
  tree->Branch( "ancestor_pprodpx",  dAncestorDotPprodpx,   "ancestor_pprodpx[ancestor_size]/D" );
  tree->Branch( "ancestor_pprodpy",  dAncestorDotPprodpy,   "ancestor_pprodpy[ancestor_size]/D" );
  tree->Branch( "ancestor_pprodpz",  dAncestorDotPprodpz,   "ancestor_pprodpz[ancestor_size]/D" );
  tree->Branch( "ancestor_nucleus",  dAncestorDotNucleus,   "ancestor_nucleus[ancestor_size]/I" );
  tree->Branch( "ancestor_proc",     dAncestorDotProc,      "ancestor_proc/C"                   );
  tree->Branch( "ancestor_ivol",     dAncestorDotIvol,      "ancestor_ivol/C"                   );
  tree->Branch( "ancestor_imat",     dAncestorDotImat,      "ancestor_imat/C"                   );
  tree->Branch( "ppvx",             &dPpvx,                 "ppvx/D"                            );
  tree->Branch( "ppvy",             &dPpvy,                 "ppvy/D"                            );
  tree->Branch( "ppvz",             &dPpvz,                 "ppvz/D"                            );
  tree->Branch( "tgtexit_tvx",      &dTgtexitDotTvx,        "tgtexit_tvx/D"                     );
  tree->Branch( "tgtexit_tvy",      &dTgtexitDotTvy,        "tgtexit_tvy/D"                     );
  tree->Branch( "tgtexit_tvz",      &dTgtexitDotTvz,        "tgtexit_tvz/D"                     );
  tree->Branch( "tgtexit_tpx",      &dTgtexitDotTpx,        "tgtexit_tpx/D"                     );
  tree->Branch( "tgtexit_tpy",      &dTgtexitDotTpy,        "tgtexit_tpy/D"                     );
  tree->Branch( "tgtexit_tpz",      &dTgtexitDotTpz,        "tgtexit_tpz/D"                     );
  tree->Branch( "tgtexit_tptype",   &dTgtexitDotTptype,     "tgtexit_tptype/I"                  );
  tree->Branch( "tgtexit_tgen",     &dTgtexitDotTgen,       "tgtexit_tgen/I"                    );
  tree->Branch( "traj_size",        &dTrArSize,             "traj_size/I"                       );
  tree->Branch( "traj_trkx",         dTrajDotTrkx,          "traj_trkx[traj_size]/D"            );
  tree->Branch( "traj_trky",         dTrajDotTrky,          "traj_trky[traj_size]/D"            );
  tree->Branch( "traj_trkz",         dTrajDotTrkz,          "traj_trkz[traj_size]/D"            );
  tree->Branch( "traj_trkpx",        dTrajDotTrkpx,         "traj_trkpx[traj_size]/D"           );
  tree->Branch( "traj_trkpy",        dTrajDotTrkpy,         "traj_trkpy[traj_size]/D"           );
  tree->Branch( "traj_trkpz",        dTrajDotTrkpz,         "traj_trkpz[traj_size]/D"           );
  //tree->Branch( "flagbits",         &dFlagbits,             "flagbits/I"                        );
  //tree->Branch( "vint",              dVint                                                      );
  //tree->Branch( "vdbl",              dVdbl                                                      );
  
  std::cout << "  Dk2nu branches initialised\n";
}

//============================================================================
// Fill branches here
//============================================================================

void FillMetaBranches(bsim::DkMeta * dkmeta)
{
  assert(dkmeta);

  mJob        = dkmeta->job;
  mPots       = dkmeta->pots;
  RootifyChar(dkmeta->beamsim,  mBeamsim );
  RootifyChar(dkmeta->physics,  mPhysics );
  RootifyChar(dkmeta->physcuts, mPhyscuts);
  RootifyChar(dkmeta->tgtcfg,   mTgtcfg  );
  RootifyChar(dkmeta->horncfg,  mHorncfg );
  RootifyChar(dkmeta->dkvolcfg, mDkvolcfg);
  mBeam0x     = dkmeta->beam0x;
  mBeam0y     = dkmeta->beam0y;
  mBeam0z     = dkmeta->beam0z;
  mBeamhwidth = dkmeta->beamhwidth;
  mBeamvwidth = dkmeta->beamvwidth;
  mBeamdxdz   = dkmeta->beamdxdz;
  mBeamdydz   = dkmeta->beamdydz;
  
  mArSize     = (dkmeta->location).size();
  std::string tmpString = std::string("");
  for( unsigned int i = 0; i < mArSize; i++ ){
    mLocationDotX[i]    = (dkmeta->location.at(i)).x;
    mLocationDotY[i]    = (dkmeta->location.at(i)).y;
    mLocationDotZ[i]    = (dkmeta->location.at(i)).z;
    //RootifyChar((dkmeta->location.at(i)).name, mLocationDotName[i]); //avoid char 2D arrays
    tmpString.append( ((dkmeta->location.at(i)).name).c_str() );
    tmpString.append( "||" );
  }
  RootifyChar( tmpString, mLocationDotName );

  //mVintnames  = dkmeta->vintnames;
  //mVdblnames  = dkmeta->vdblnames;

  dkMeta->Fill();
}

void FillTreeBranches(bsim::Dk2Nu * dk2nu)
{
  dJob              = dk2nu->job;
  dPotnum           = dk2nu->potnum;
  dPpvx             = dk2nu->ppvx;
  dPpvy             = dk2nu->ppvy;
  dPpvz             = dk2nu->ppvz;
  // - - -
  dDecayDotNorig    = dk2nu->decay.norig;
  dDecayDotNdecay   = dk2nu->decay.ndecay;
  dDecayDotNtype    = dk2nu->decay.ntype;
  dDecayDotVx       = dk2nu->decay.vx;
  dDecayDotVy       = dk2nu->decay.vy;
  dDecayDotVz       = dk2nu->decay.vz;
  dDecayDotPdpx     = dk2nu->decay.pdpx;
  dDecayDotPdpy     = dk2nu->decay.pdpy;
  dDecayDotPdpz     = dk2nu->decay.pdpz;
  dDecayDotPpdxdz   = dk2nu->decay.ppdxdz;
  dDecayDotPpdydz   = dk2nu->decay.ppdydz;
  dDecayDotPppz     = dk2nu->decay.pppz;
  dDecayDotPpenergy = dk2nu->decay.ppenergy;
  dDecayDotPpmedium = dk2nu->decay.ppmedium;
  dDecayDotPtype    = dk2nu->decay.ptype;
  dDecayDotMuparpx  = dk2nu->decay.muparpx;
  dDecayDotMuparpy  = dk2nu->decay.muparpy;
  dDecayDotMuparpz  = dk2nu->decay.muparpz;
  dDecayDotMupare   = dk2nu->decay.mupare;
  dDecayDotNecm     = dk2nu->decay.necm;
  dDecayDotNimpwt   = dk2nu->decay.nimpwt;
  // - - -
  dArSize           = (dk2nu->nuray).size();
  for( unsigned int i = 0; i < dArSize; i++ ){
    dNurayDotPx[i]  = (dk2nu->nuray.at(i)).px;
    dNurayDotPy[i]  = (dk2nu->nuray.at(i)).py;
    dNurayDotPz[i]  = (dk2nu->nuray.at(i)).pz;
    dNurayDotE[i]   = (dk2nu->nuray.at(i)).E;
    dNurayDotWgt[i] = (dk2nu->nuray.at(i)).wgt;
  }
  // - - -
  dAnArSize         = (dk2nu->ancestor).size();
  std::string tmpStringProc = std::string("");
  std::string tmpStringIvol = std::string("");
  std::string tmpStringImat = std::string("");
  for( unsigned int i = 0; i < dAnArSize; i++ ){
    dAncestorDotPdg[i]     = (dk2nu->ancestor.at(i)).pdg;
    dAncestorDotStartx[i]  = (dk2nu->ancestor.at(i)).startx;
    dAncestorDotStarty[i]  = (dk2nu->ancestor.at(i)).starty;
    dAncestorDotStartz[i]  = (dk2nu->ancestor.at(i)).startz;
    dAncestorDotStartt[i]  = (dk2nu->ancestor.at(i)).startt;
    dAncestorDotStartpx[i] = (dk2nu->ancestor.at(i)).startpx;
    dAncestorDotStartpy[i] = (dk2nu->ancestor.at(i)).startpy;
    dAncestorDotStartpz[i] = (dk2nu->ancestor.at(i)).startpz;
    dAncestorDotStoppx[i]  = (dk2nu->ancestor.at(i)).stoppx;
    dAncestorDotStoppy[i]  = (dk2nu->ancestor.at(i)).stoppy;
    dAncestorDotStoppz[i]  = (dk2nu->ancestor.at(i)).stoppz;
    dAncestorDotPolx[i]    = (dk2nu->ancestor.at(i)).polx;
    dAncestorDotPoly[i]    = (dk2nu->ancestor.at(i)).poly;
    dAncestorDotPolz[i]    = (dk2nu->ancestor.at(i)).polz;
    dAncestorDotPprodpx[i] = (dk2nu->ancestor.at(i)).pprodpx;
    dAncestorDotPprodpy[i] = (dk2nu->ancestor.at(i)).pprodpy;
    dAncestorDotPprodpz[i] = (dk2nu->ancestor.at(i)).pprodpz;
    dAncestorDotNucleus[i] = (dk2nu->ancestor.at(i)).nucleus;
    //RootifyChar((dk2nu->ancestor.at(i)).proc, dAncestorDotProc[i]);
    //RootifyChar((dk2nu->ancestor.at(i)).ivol, dAncestorDotIvol[i]);
    //RootifyChar((dk2nu->ancestor.at(i)).imat, dAncestorDotImat[i]);
    tmpStringProc.append( ((dk2nu->ancestor.at(i)).proc).c_str() );
    tmpStringIvol.append( ((dk2nu->ancestor.at(i)).ivol).c_str() );
    tmpStringImat.append( ((dk2nu->ancestor.at(i)).imat).c_str() );
    tmpStringProc.append( "||" );
    tmpStringIvol.append( "||" );
    tmpStringImat.append( "||" );
  }
  RootifyChar( tmpStringProc, dAncestorDotProc );
  RootifyChar( tmpStringIvol, dAncestorDotIvol );
  RootifyChar( tmpStringImat, dAncestorDotImat );
  // - - -
  dTgtexitDotTvx     = dk2nu->tgtexit.tvx;
  dTgtexitDotTvy     = dk2nu->tgtexit.tvy;
  dTgtexitDotTvz     = dk2nu->tgtexit.tvz;
  dTgtexitDotTpx     = dk2nu->tgtexit.tpx;
  dTgtexitDotTpy     = dk2nu->tgtexit.tpy;
  dTgtexitDotTpz     = dk2nu->tgtexit.tpz;
  dTgtexitDotTptype  = dk2nu->tgtexit.tptype;
  dTgtexitDotTgen    = dk2nu->tgtexit.tgen;
  // - - -
  dTrArSize          = (dk2nu->traj).size();
  for( unsigned int i = 0; i < dTrArSize; i++ ){
    dTrajDotTrkx[i]  = (dk2nu->traj.at(i)).trkx;
    dTrajDotTrky[i]  = (dk2nu->traj.at(i)).trky;
    dTrajDotTrkz[i]  = (dk2nu->traj.at(i)).trkz;
    dTrajDotTrkpx[i] = (dk2nu->traj.at(i)).trkpx;
    dTrajDotTrkpy[i] = (dk2nu->traj.at(i)).trkpy;
    dTrajDotTrkpz[i] = (dk2nu->traj.at(i)).trkpz;
  }
  // - - -
  //dFlagbits          = dk2nu->flagbits;
  //dVint              = dk2nu->vint;
  //dVdbl              = dk2nu->vdbl;

  dkTree->Fill();
}

//============================================================================
// Utility method for chars
//============================================================================

void RootifyChar( std::string rfch, char fdch[maxC] )
{
  //std::cout << "\n*\n" << rfch.c_str() << "\n";
  std::string tmpString(rfch.c_str());
  tmpString.append("\0");
  //std::cout << tmpString.c_str() << "\n*\n";
  sprintf( fdch, "%s", tmpString.c_str() );
}
