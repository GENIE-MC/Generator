//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Alfonso Garcia <alfonsog \at nikhef.nl>
         NIKHEF

 For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________

#include "Physics/HEDIS/XSection/HEDISStrucFunc.h"
#include "Physics/HEDIS/EventGen/HEDISInteractionListGenerator.h"
#include "Framework/EventGen/InteractionList.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Utils/RunOpt.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"

#include <TSystem.h>
#include <TMath.h>

#ifdef __GENIE_APFEL_ENABLED__
#include "APFEL/APFEL.h"
#endif
#include "LHAPDF/LHAPDF.h"
#ifdef __GENIE_LHAPDF6_ENABLED__
LHAPDF::PDF* pdf;
#endif

using namespace genie;
using namespace genie::constants;

// values from LHAPDF set
double xPDFmin;                  // Minimum values of x in grid from LHPADF set
double Q2PDFmin;                 // Minimum values of Q2 in grid from LHPADF set
double Q2PDFmax;                 // Maximum values of Q2 in grid from LHPADF set
std::map<int, double> mPDFQrk;   // Mass of the quark from LHAPDF set

//_________________________________________________________________________
HEDISStrucFunc * HEDISStrucFunc::fgInstance = 0;
//_________________________________________________________________________
HEDISStrucFunc::HEDISStrucFunc(SF_info sfinfo)
{

  fSF = sfinfo;

  string basedir = "";
  if ( gSystem->Getenv("HEDIS_SF_DATA_PATH")==NULL ) basedir = string(gSystem->Getenv("GENIE")) + "/data/evgen/hedis-sf";
  else                                               basedir = string(gSystem->Getenv("HEDIS_SF_DATA_PATH"));
  LOG("HEDISStrucFunc", pWARN) << "Base directory: " << basedir;

  if ( gSystem->AccessPathName( basedir.c_str(), kReadPermission ) ) {
      LOG("HEDISStrucFunc", pFATAL) << "Base directory doesnt exist or you dont have read permission.";
      LOG("HEDISStrucFunc", pFATAL) << "Remember!!!";
      LOG("HEDISStrucFunc", pFATAL) << "Path to base directory is defined with the enviroment variable HEDIS_SF_DATA_PATH.";
      LOG("HEDISStrucFunc", pFATAL) << "If not defined, default location is $GENIE/data/evgen/hedis-sf";
      assert(0);
  }

  // Name of the directory where SF tables are (or will be) saved. 
  string SFname = basedir + "/" + RunOpt::Instance()->Tune()->Name();

  // Check that the directory where SF tables are stored exists
  LOG("HEDISStrucFunc", pWARN) << "SF are (or will be) in following directory: " << SFname;
  if ( gSystem->AccessPathName( SFname.c_str(), kReadPermission ) ) {
    LOG("HEDISStrucFunc", pWARN) << "Bad news! Directory doesnt exists.";
    LOG("HEDISStrucFunc", pWARN) << "HEDIS package requires computation of SF";        
    LOG("HEDISStrucFunc", pWARN) << "This will be SLOW!!!!!";        
    if ( gSystem->mkdir(SFname.c_str())==0 ) {
      LOG("HEDISStrucFunc", pINFO) << "Creating Metafile with input information";
      std::ofstream meta_stream((SFname+"/Inputs.txt").c_str());
      meta_stream << fSF;
      meta_stream.close();      
    }
    else {
      LOG("HEDISStrucFunc", pFATAL) << "You dont have write permission in the following directory:";
      LOG("HEDISStrucFunc", pFATAL) << SFname;
      LOG("HEDISStrucFunc", pFATAL) << "Remember!!!";
      LOG("HEDISStrucFunc", pFATAL) << "Path to base directory is defined with the enviroment variable HEDIS_SF_DATA_PATH.";
      LOG("HEDISStrucFunc", pFATAL) << "If not defined, default location is $GENIE/data/evgen/hedis-sf";
      assert(0);
    }
  }
  else {
    LOG("HEDISStrucFunc", pWARN) << "Good news! Directory already exists.";
    LOG("HEDISStrucFunc", pWARN) << "Precomputed SF stored in that directory will be used";
    LOG("HEDISStrucFunc", pINFO) << "Comparing info from Metafile and selected Tune (CommonParam.xml)";
    SF_info cm;
    std::ifstream meta_stream((SFname+"/Inputs.txt").c_str(), std::ios::in);
    meta_stream >> cm;
    meta_stream.close();
    if (cm==fSF) {
      LOG("HEDISStrucFunc", pINFO) << "Info from MetaFile and Tune match";
    }
    else {
      LOG("HEDISStrucFunc", pFATAL) << "Info from MetaFile and Tune doesnt match";
      LOG("HEDISStrucFunc", pFATAL) << "MetaFile Path : " << SFname << "/Inputs.txt";        
      LOG("HEDISStrucFunc", pFATAL) << "From Tune : ";        
      LOG("HEDISStrucFunc", pFATAL) << cm;        
      LOG("HEDISStrucFunc", pFATAL) << "From MetaFile : ";        
      LOG("HEDISStrucFunc", pFATAL) << fSF;        
      assert(0);
    }
  }

  fSF.Sin2ThW = fSF.Rho==0. ? fSF.Sin2ThW : 1.-fSF.MassW*fSF.MassW/fSF.MassZ/fSF.MassZ/(1+fSF.Rho);
  LOG("HEDISStrucFunc", pINFO) << "Weingber angle computation:";
  LOG("HEDISStrucFunc", pINFO) << "Rho     = " << fSF.Rho;
  LOG("HEDISStrucFunc", pINFO) << "Sin2ThW = " << fSF.Sin2ThW;

  // initialising lhapdf
#ifdef __GENIE_LHAPDF6_ENABLED__
  LOG("HEDISStrucFunc", pINFO) << "Initialising LHAPDF6...";
  pdf = LHAPDF::mkPDF(fSF.LHAPDFset, fSF.LHAPDFmember);
  xPDFmin  = pdf->xMin();
  Q2PDFmin = pdf->q2Min();
  Q2PDFmax = pdf->q2Max();
  LOG("HEDISStrucFunc", pINFO) << "PDF info:";
  LOG("HEDISStrucFunc", pINFO) << "OrderQCD = " << pdf->orderQCD();
  LOG("HEDISStrucFunc", pINFO) << "FlavorScheme = " << pdf->info().get_entry("FlavorScheme");
  LOG("HEDISStrucFunc", pINFO) << "NumFlavors = " << pdf->info().get_entry("NumFlavors");
  LOG("HEDISStrucFunc", pINFO) << "Xmin = " << xPDFmin << "  Xmax = " << pdf->xMax() << "  Q2min = " << Q2PDFmin << "  Q2max = " << Q2PDFmax;
  LOG("HEDISStrucFunc", pINFO) << "MZ = " << pdf->info().get_entry("MZ");
  for (int i=1; i<7; i++) {
    mPDFQrk[i] = pdf->quarkMass(i);
    LOG("HEDISStrucFunc", pINFO) << "M" << i << " = " << mPDFQrk[i];
  }
#endif
#ifdef __GENIE_LHAPDF5_ENABLED__
  LOG("HEDISStrucFunc", pINFO) << "Initialising LHAPDF5...";
  LHAPDF::initPDFByName(fSF.LHAPDFset, LHAPDF::LHGRID, fSF.LHAPDFmember);
  xPDFmin  = LHAPDF::getXmin(0);
  Q2PDFmin = LHAPDF::getQ2min(0);
  Q2PDFmax = LHAPDF::getQ2max(0);
  LOG("HEDISStrucFunc", pINFO) << "PDF info:";
  LOG("HEDISStrucFunc", pINFO) << "Xmin = " << xPDFmin << "  Xmax = " << LHAPDF::getXmax(0) << "  Q2min = " << Q2PDFmin << "  Q2max = " << Q2PDFmax;
  for (int i=1; i<7; i++) {
    mPDFQrk[i] = LHAPDF::getQMass(i);
    LOG("HEDISStrucFunc", pINFO) << "M" << i << " = " << mPDFQrk[i];
  }
#endif

  // checking if LHAPDF boundaries matches boundaries defined by user
  if ( fSF.XGridMin < xPDFmin ) {
    LOG("HEDISStrucFunc", pWARN) << "Lower boundary in X is smaller than input PDF";
    LOG("HEDISStrucFunc", pWARN) << "xPDFmin = "  << xPDFmin;
    LOG("HEDISStrucFunc", pWARN) << "xGridMin = " << fSF.XGridMin;
    LOG("HEDISStrucFunc", pWARN) << "DANGER: An extrapolation will be done!!!!!";
  }
  else if ( fSF.Q2GridMin < Q2PDFmin ) {
    LOG("HEDISStrucFunc", pWARN) << "Lower boundary in Q2 is smaller than input PDF";
    LOG("HEDISStrucFunc", pWARN) << "Q2PDFmin = "  << Q2PDFmin;
    LOG("HEDISStrucFunc", pWARN) << "Q2GridMin = " << fSF.Q2GridMin;
    LOG("HEDISStrucFunc", pWARN) << "DANGER: An extrapolation will be done!!!!!";
  }
  else if ( fSF.Q2GridMax > Q2PDFmax ) {
    LOG("HEDISStrucFunc", pWARN) << "Upper boundary in Q2 is bigger than input PDF";
    LOG("HEDISStrucFunc", pWARN) << "Q2PDFmax = "  << Q2PDFmax;
    LOG("HEDISStrucFunc", pWARN) << "Q2GridMax = " << fSF.Q2GridMax;
    LOG("HEDISStrucFunc", pWARN) << "DANGER: An extrapolation will be done!!!!!";
  }

  // define arrays to fill from data files
  double dlogq2 = TMath::Abs( TMath::Log10(fSF.Q2GridMin)-TMath::Log10(fSF.Q2GridMax) ) / fSF.NGridQ2;
  double dlogx  = TMath::Abs( TMath::Log10(fSF.XGridMin)-TMath::Log10(1.) ) / fSF.NGridX;
  LOG("HEDISStrucFunc", pINFO) << "Grid x,Q2 :" << fSF.NGridX << " , " << fSF.NGridQ2;
  for ( double logq2 = TMath::Log10(fSF.Q2GridMin); logq2<TMath::Log10(fSF.Q2GridMax); logq2+= dlogq2 ) {
    double q2 = TMath::Power( 10, logq2 + 0.5*dlogq2 );
    if (q2>fSF.Q2GridMax) continue;
    sf_q2_array.push_back(q2);
    LOG("HEDISStrucFunc", pDEBUG) << "q2: " << sf_q2_array.back();
  }
  for ( double logx = TMath::Log10(fSF.XGridMin); logx<TMath::Log10(1.); logx+= dlogx ) {
    double x = TMath::Power( 10, logx + 0.5*dlogx );
    if ( x>1. ) continue;
    sf_x_array.push_back(x);
    LOG("HEDISStrucFunc", pDEBUG) << "x: " << sf_x_array.back();
  }

  // Change to variables that are suitable for BLI2DNonUnifGrid
  int nx = sf_q2_array.size();
  int ny = sf_x_array.size();
  double x[nx];
  double y[ny];
  double z[nx*ny];
  for (int i=0; i<nx; i++) x[i] = sf_q2_array[i];
  for (int j=0; j<ny; j++) y[j] = sf_x_array[j];
  
  vector<InteractionType_t> inttype;
  inttype.push_back(kIntWeakCC);
  inttype.push_back(kIntWeakNC);
  vector<InitialState> init_state;
  init_state.push_back(InitialState(1, 2, kPdgNuE));
  init_state.push_back(InitialState(1, 2, kPdgAntiNuE));

  HEDISInteractionListGenerator * helist = new HEDISInteractionListGenerator();
  InteractionList * ilist = helist->CreateHEDISlist(init_state,inttype);

  // Load structure functions for each quark at LO
  for(InteractionList::iterator in=ilist->begin(); in!=ilist->end(); ++in) {

    string sfFile = SFname + "/QrkSF_LO_" + QrkSFName(*in) + ".dat";
    // Make sure data files are available
    LOG("HEDISStrucFunc", pINFO) << "Checking if file " << sfFile << " exists...";        
    if ( gSystem->AccessPathName( sfFile.c_str()) ) {
      LOG("HEDISStrucFunc", pWARN) << "File doesnt exist. SF table will be computed.";        
      CreateQrkSF( *in, sfFile );
    }
    else if ( atoi(gSystem->GetFromPipe(("wc -w "+sfFile+" | awk '{print $1}'").c_str()))!=kSFT3*nx*ny ) {
      LOG("HEDISStrucFunc", pWARN) << "File does not contain all the need points. SF table will be recomputed.";        
      gSystem->Exec(("rm "+sfFile).c_str());
      CreateQrkSF( *in, sfFile );
    }
    std::ifstream sf_stream(sfFile.c_str(), std::ios::in);
    // Loop over F1,F2,F3
    for(int sf = 1; sf < kSFnumber; ++sf) {
      // Loop over x/Q2 bins
      for ( int ij=0; ij<nx*ny; ij++ ) sf_stream >> z[ij];
      // Create SF tables with BLI2DNonUnifGrid using x,Q2 binning
      fQrkSFLOTables[QrkSFCode(*in)].Table[(HEDISStrucFuncType_t)sf] = new genie::BLI2DNonUnifGrid( nx, ny, x, y, z );
    }
  }

  if (fSF.IsNLO) {
#ifdef __GENIE_APFEL_ENABLED__
    // initialising APFEL framework
    LOG("HEDISStrucFunc", pINFO) << "Initialising APFEL..." ; 
    APFEL::SetPDFSet(fSF.LHAPDFset);
    APFEL::SetReplica(fSF.LHAPDFmember);
    if (fSF.Scheme=="BGR") {
      APFEL::SetMassScheme("FONLL-B");
      APFEL::SetPoleMasses(mPDFQrk[4],mPDFQrk[5],mPDFQrk[6]);
    } 
    else if (fSF.Scheme=="CSMS") {
      APFEL::SetMassScheme("ZM-VFNS");
      APFEL::SetPoleMasses(mPDFQrk[4],mPDFQrk[5],mPDFQrk[5]+0.1);
    }
    else if (fSF.Scheme=="GGHR") {
      APFEL::SetMassScheme("FFNS5");
      APFEL::SetPoleMasses(mPDFQrk[4],mPDFQrk[5],mPDFQrk[6]);
      APFEL::SetMaxFlavourPDFs(5);
      APFEL::SetMaxFlavourAlpha(5);
    }
    else {
      LOG("HEDISStrucFunc", pERROR) << "Mass Scheme is not set properly";
      assert(0);
    }
    APFEL::SetQLimits(TMath::Sqrt(fSF.Q2GridMin),TMath::Sqrt(fSF.Q2GridMax));
    APFEL::SetMaxFlavourPDFs(6);
    APFEL::SetMaxFlavourAlpha(6);
    APFEL::SetNumberOfGrids(3);
    APFEL::SetGridParameters(1,90,3,fSF.XGridMin);
    APFEL::SetGridParameters(2,50,5,1e-1);
    APFEL::SetGridParameters(3,40,5,8e-1);
    APFEL::SetPerturbativeOrder(1);
    APFEL::SetAlphaQCDRef(pdf->alphasQ(fSF.MassZ),fSF.MassZ);    
    APFEL::SetProtonMass(kProtonMass);
    APFEL::SetWMass(fSF.MassW);
    APFEL::SetZMass(fSF.MassZ);
    APFEL::SetSin2ThetaW(fSF.Sin2ThW);
    APFEL::SetCKM(fSF.Vud, fSF.Vus, fSF.Vub,
                  fSF.Vcd, fSF.Vcs, fSF.Vcb,
                  fSF.Vtd, fSF.Vts, fSF.Vtb);
#endif

    // Load structure functions for each quark at LO
    int nch = -1;
    for(InteractionList::iterator in=ilist->begin(); in!=ilist->end(); ++in) {

      if ( nch==NucSFCode(*in) ) continue;
      nch = NucSFCode(*in);
      
      string sfFile = SFname + "/NucSF_NLO_" + NucSFName(*in) + ".dat";
      // Make sure data files are available
      LOG("HEDISStrucFunc", pINFO) << "Checking if file " << sfFile << " exists...";        
      if ( gSystem->AccessPathName( sfFile.c_str()) ) {
#ifdef __GENIE_APFEL_ENABLED__
        LOG("HEDISStrucFunc", pWARN) << "File doesnt exist. SF table will be computed.";        
        CreateNucSF( *in, sfFile );
#else
        LOG("HEDISStrucFunc", pERROR) << "File doesnt exist. APFEL is needed for NLO SF";        
        assert(0);
#endif
      }
      else if ( atoi(gSystem->GetFromPipe(("wc -w "+sfFile+" | awk '{print $1}'").c_str()))!=kSFT3*nx*ny ) {
#ifdef __GENIE_APFEL_ENABLED__
        LOG("HEDISStrucFunc", pWARN) << "File does not contain all the need points. SF table will be recomputed.";        
        gSystem->Exec(("rm "+sfFile).c_str());
        CreateNucSF( *in, sfFile );
#else
        LOG("HEDISStrucFunc", pERROR) << "File does not contain all the need points. APFEL is needed for NLO SF";        
        assert(0);
#endif
      }
      std::ifstream sf_stream(sfFile.c_str(), std::ios::in);
      // Loop over F1,F2,F3
      for(int sf = 1; sf < kSFnumber; ++sf) {
        // Loop over x/Q2 bins
        for ( int ij=0; ij<nx*ny; ij++ ) sf_stream >> z[ij];
        // Create SF tables with BLI2DNonUnifGrid using x,Q2 binning
        fNucSFNLOTables[nch].Table[(HEDISStrucFuncType_t)sf] = new genie::BLI2DNonUnifGrid( nx, ny, x, y, z );
      }        

      //compute structure functions for each nucleon at LO using quark grids
      LOG("HEDISStrucFunc", pDEBUG) << "Creating LO " << sfFile;              
      vector <int> qcodes;
      for(InteractionList::iterator in2=ilist->begin(); in2!=ilist->end(); ++in2) {
        if (NucSFCode(*in2)==nch) qcodes.push_back(QrkSFCode(*in2));
      }
      // Loop over F1,F2,F3
      for(int sf = 1; sf < kSFnumber; ++sf) {
        int ij = 0;
        // Loop over Q2 bins
        for (int i=0; i<nx; i++) {
          // Loop over x bins
          for (int j=0; j<ny; j++) {
            double sum = 0.;
            // NucSF = sum_qrks QrkSF
            for(vector<int>::const_iterator iq=qcodes.begin(); iq!=qcodes.end(); ++iq) 
              sum += fQrkSFLOTables[*iq].Table[(HEDISStrucFuncType_t)sf]->Evaluate(x[i],y[j]);
            z[ij] = sum;
            ij++;
          }
        }
        // Create SF tables with BLI2DNonUnifGrid using x,Q2 binning
        fNucSFLOTables[nch].Table[(HEDISStrucFuncType_t)sf] = new genie::BLI2DNonUnifGrid( nx, ny, x, y, z );
      }
    }
  }
  
  fgInstance = 0;

}
//_________________________________________________________________________
HEDISStrucFunc::~HEDISStrucFunc()
{

}
//_________________________________________________________________________
HEDISStrucFunc * HEDISStrucFunc::Instance(SF_info sfinfo)
{
  if(fgInstance == 0) {
    LOG("HEDISStrucFunc", pINFO) << "Late initialization";
    static HEDISStrucFunc::Cleaner cleaner;
    cleaner.DummyMethodAndSilentCompiler();
    fgInstance = new HEDISStrucFunc(sfinfo);
  }  
  return fgInstance;
}


//____________________________________________________________________________
void HEDISStrucFunc::CreateQrkSF( const Interaction * in, string sfFile ) 
{

  // variables used to tag the SF for particular channel
  bool iscc        = in->ProcInfo().IsWeakCC();
  bool isnu        = pdg::IsNeutrino(in->InitState().ProbePdg());
  bool ispr        = pdg::IsProton(in->InitState().Tgt().HitNucPdg());
  bool sea_iq      = in->InitState().Tgt().HitSeaQrk();
  int pdg_iq       = in->InitState().Tgt().HitQrkPdg();
  int pdg_fq       = in->ExclTag().FinalQuarkPdg();
  double mass_nucl = in->InitState().Tgt().HitNucMass();

  // up and down quark swicth depending on proton or neutron interaction  
  int qrkd = 0;
  int qrku = 0;
  if ( ispr ) { qrkd = 1 ; qrku = 2; }
  else        { qrkd = 2 ; qrku = 1; }

  // variables associated to the PDF and coupling of the quarks
  int qpdf1 = -999;                                     // number used to compute PDF
  int qpdf2 = -999;                                     // auxiliary number used to compute PDF for valence quarks
  double Cp2 = -999;                                    // couping for F1,F2
  double Cp3 = -999;                                    // couping for F3
  double sign3 = isnu ? +1. : -1.;  // sign change for nu/nubar in F3
  if ( iscc ) {
    if ( isnu ) {
      if      ( pdg_iq== 1 && !sea_iq && pdg_fq== 2 ) { qpdf1 =  qrkd; qpdf2 = -qrkd; Cp2 = 2*fSF.Vud*fSF.Vud; Cp3 =  2*fSF.Vud*fSF.Vud; }
      else if ( pdg_iq== 1 && !sea_iq && pdg_fq== 4 ) { qpdf1 =  qrkd; qpdf2 = -qrkd; Cp2 = 2*fSF.Vcd*fSF.Vcd; Cp3 =  2*fSF.Vcd*fSF.Vcd; }
      else if ( pdg_iq== 1 && !sea_iq && pdg_fq== 6 ) { qpdf1 =  qrkd; qpdf2 = -qrkd; Cp2 = 2*fSF.Vtd*fSF.Vtd; Cp3 =  2*fSF.Vtd*fSF.Vtd; }
      else if ( pdg_iq== 1 &&  sea_iq && pdg_fq== 2 ) { qpdf1 = -qrkd;                Cp2 = 2*fSF.Vud*fSF.Vud; Cp3 =  2*fSF.Vud*fSF.Vud; }
      else if ( pdg_iq== 1 &&  sea_iq && pdg_fq== 4 ) { qpdf1 = -qrkd;                Cp2 = 2*fSF.Vcd*fSF.Vcd; Cp3 =  2*fSF.Vcd*fSF.Vcd; }
      else if ( pdg_iq== 1 &&  sea_iq && pdg_fq== 6 ) { qpdf1 = -qrkd;                Cp2 = 2*fSF.Vtd*fSF.Vtd; Cp3 =  2*fSF.Vtd*fSF.Vtd; }
      else if ( pdg_iq== 3 &&  sea_iq && pdg_fq== 2 ) { qpdf1 =  3;                   Cp2 = 2*fSF.Vus*fSF.Vus; Cp3 =  2*fSF.Vus*fSF.Vus; }
      else if ( pdg_iq== 3 &&  sea_iq && pdg_fq== 4 ) { qpdf1 =  3;                   Cp2 = 2*fSF.Vcs*fSF.Vcs; Cp3 =  2*fSF.Vcs*fSF.Vcs; }
      else if ( pdg_iq== 3 &&  sea_iq && pdg_fq== 6 ) { qpdf1 =  3;                   Cp2 = 2*fSF.Vts*fSF.Vts; Cp3 =  2*fSF.Vts*fSF.Vts; }
      else if ( pdg_iq== 5 &&  sea_iq && pdg_fq== 2 ) { qpdf1 =  5;                   Cp2 = 2*fSF.Vub*fSF.Vub; Cp3 =  2*fSF.Vub*fSF.Vub; }
      else if ( pdg_iq== 5 &&  sea_iq && pdg_fq== 4 ) { qpdf1 =  5;                   Cp2 = 2*fSF.Vcb*fSF.Vcb; Cp3 =  2*fSF.Vcb*fSF.Vcb; }
      else if ( pdg_iq== 5 &&  sea_iq && pdg_fq== 6 ) { qpdf1 =  5;                   Cp2 = 2*fSF.Vtb*fSF.Vtb; Cp3 =  2*fSF.Vtb*fSF.Vtb; }
      else if ( pdg_iq==-2 &&  sea_iq && pdg_fq==-1 ) { qpdf1 = -qrku;                Cp2 = 2*fSF.Vud*fSF.Vud; Cp3 = -2*fSF.Vud*fSF.Vud; }
      else if ( pdg_iq==-2 &&  sea_iq && pdg_fq==-3 ) { qpdf1 = -qrku;                Cp2 = 2*fSF.Vus*fSF.Vus; Cp3 = -2*fSF.Vus*fSF.Vus; }
      else if ( pdg_iq==-2 &&  sea_iq && pdg_fq==-5 ) { qpdf1 = -qrku;                Cp2 = 2*fSF.Vub*fSF.Vub; Cp3 = -2*fSF.Vub*fSF.Vub; }
      else if ( pdg_iq==-4 &&  sea_iq && pdg_fq==-1 ) { qpdf1 = -4;                   Cp2 = 2*fSF.Vcd*fSF.Vcd; Cp3 = -2*fSF.Vcd*fSF.Vcd; }
      else if ( pdg_iq==-4 &&  sea_iq && pdg_fq==-3 ) { qpdf1 = -4;                   Cp2 = 2*fSF.Vcs*fSF.Vcs; Cp3 = -2*fSF.Vcs*fSF.Vcs; }
      else if ( pdg_iq==-4 &&  sea_iq && pdg_fq==-5 ) { qpdf1 = -4;                   Cp2 = 2*fSF.Vcb*fSF.Vcb; Cp3 = -2*fSF.Vcb*fSF.Vcb; }
    }
    else {
      if      ( pdg_iq== 2 && !sea_iq && pdg_fq== 1 ) { qpdf1 =  qrku; qpdf2 = -qrku; Cp2 = 2*fSF.Vud*fSF.Vud; Cp3 =  2*fSF.Vud*fSF.Vud; }
      else if ( pdg_iq== 2 && !sea_iq && pdg_fq== 3 ) { qpdf1 =  qrku; qpdf2 = -qrku; Cp2 = 2*fSF.Vus*fSF.Vus; Cp3 =  2*fSF.Vus*fSF.Vus; }
      else if ( pdg_iq== 2 && !sea_iq && pdg_fq== 5 ) { qpdf1 =  qrku; qpdf2 = -qrku; Cp2 = 2*fSF.Vub*fSF.Vub; Cp3 =  2*fSF.Vub*fSF.Vub; }
      else if ( pdg_iq== 2 &&  sea_iq && pdg_fq== 1 ) { qpdf1 = -qrku;                Cp2 = 2*fSF.Vud*fSF.Vud; Cp3 =  2*fSF.Vud*fSF.Vud; }
      else if ( pdg_iq== 2 &&  sea_iq && pdg_fq== 3 ) { qpdf1 = -qrku;                Cp2 = 2*fSF.Vus*fSF.Vus; Cp3 =  2*fSF.Vus*fSF.Vus; }
      else if ( pdg_iq== 2 &&  sea_iq && pdg_fq== 5 ) { qpdf1 = -qrku;                Cp2 = 2*fSF.Vub*fSF.Vub; Cp3 =  2*fSF.Vub*fSF.Vub; }
      else if ( pdg_iq== 4 &&  sea_iq && pdg_fq== 1 ) { qpdf1 =  4;                   Cp2 = 2*fSF.Vcd*fSF.Vcd; Cp3 =  2*fSF.Vcd*fSF.Vcd; }
      else if ( pdg_iq== 4 &&  sea_iq && pdg_fq== 3 ) { qpdf1 =  4;                   Cp2 = 2*fSF.Vcs*fSF.Vcs; Cp3 =  2*fSF.Vcs*fSF.Vcs; }
      else if ( pdg_iq== 4 &&  sea_iq && pdg_fq== 5 ) { qpdf1 =  4;                   Cp2 = 2*fSF.Vcb*fSF.Vcb; Cp3 =  2*fSF.Vcb*fSF.Vcb; }
      else if ( pdg_iq==-1 &&  sea_iq && pdg_fq==-2 ) { qpdf1 = -qrkd;                Cp2 = 2*fSF.Vud*fSF.Vud; Cp3 = -2*fSF.Vud*fSF.Vud; }
      else if ( pdg_iq==-1 &&  sea_iq && pdg_fq==-4 ) { qpdf1 = -qrkd;                Cp2 = 2*fSF.Vcd*fSF.Vcd; Cp3 = -2*fSF.Vcd*fSF.Vcd; }
      else if ( pdg_iq==-1 &&  sea_iq && pdg_fq==-6 ) { qpdf1 = -qrkd;                Cp2 = 2*fSF.Vtd*fSF.Vtd; Cp3 = -2*fSF.Vtd*fSF.Vtd; }
      else if ( pdg_iq==-3 &&  sea_iq && pdg_fq==-2 ) { qpdf1 = -3;                   Cp2 = 2*fSF.Vus*fSF.Vus; Cp3 = -2*fSF.Vus*fSF.Vus; }
      else if ( pdg_iq==-3 &&  sea_iq && pdg_fq==-4 ) { qpdf1 = -3;                   Cp2 = 2*fSF.Vcs*fSF.Vcs; Cp3 = -2*fSF.Vcs*fSF.Vcs; }
      else if ( pdg_iq==-3 &&  sea_iq && pdg_fq==-6 ) { qpdf1 = -3;                   Cp2 = 2*fSF.Vts*fSF.Vts; Cp3 = -2*fSF.Vts*fSF.Vts; }
      else if ( pdg_iq==-5 &&  sea_iq && pdg_fq==-2 ) { qpdf1 = -5;                   Cp2 = 2*fSF.Vub*fSF.Vub; Cp3 = -2*fSF.Vub*fSF.Vub; }
      else if ( pdg_iq==-5 &&  sea_iq && pdg_fq==-4 ) { qpdf1 = -5;                   Cp2 = 2*fSF.Vcb*fSF.Vcb; Cp3 = -2*fSF.Vcb*fSF.Vcb; }
      else if ( pdg_iq==-5 &&  sea_iq && pdg_fq==-6 ) { qpdf1 = -5;                   Cp2 = 2*fSF.Vtb*fSF.Vtb; Cp3 = -2*fSF.Vtb*fSF.Vtb; }
    }
  }
  else {           
    double c2u = TMath::Power( 1./2. - 4./3.*fSF.Sin2ThW,2) + 1./4.;
    double c2d = TMath::Power(-1./2. + 2./3.*fSF.Sin2ThW,2) + 1./4.;
    double c3u = 1./2. - 4./3.*fSF.Sin2ThW;
    double c3d = 1./2. - 2./3.*fSF.Sin2ThW;
    if      ( pdg_iq== 1 && !sea_iq && pdg_fq== 1 ) { qpdf1 =  qrkd; qpdf2 = -qrkd; Cp2 = c2d; Cp3 =  c3d; }
    else if ( pdg_iq== 2 && !sea_iq && pdg_fq== 2 ) { qpdf1 =  qrku; qpdf2 = -qrku; Cp2 = c2u; Cp3 =  c3u; }
    else if ( pdg_iq== 1 &&  sea_iq && pdg_fq== 1 ) { qpdf1 = -qrkd;                Cp2 = c2d; Cp3 =  c3d; }
    else if ( pdg_iq== 2 &&  sea_iq && pdg_fq== 2 ) { qpdf1 = -qrku;                Cp2 = c2u; Cp3 =  c3u; }
    else if ( pdg_iq== 3 &&  sea_iq && pdg_fq== 3 ) { qpdf1 =  3;                   Cp2 = c2d; Cp3 =  c3d; }
    else if ( pdg_iq== 4 &&  sea_iq && pdg_fq== 4 ) { qpdf1 =  4;                   Cp2 = c2u; Cp3 =  c3u; }
    else if ( pdg_iq== 5 &&  sea_iq && pdg_fq== 5 ) { qpdf1 =  5;                   Cp2 = c2d; Cp3 =  c3d; }
    else if ( pdg_iq==-1 &&  sea_iq && pdg_fq==-1 ) { qpdf1 = -qrkd;                Cp2 = c2d; Cp3 = -c3d; }
    else if ( pdg_iq==-2 &&  sea_iq && pdg_fq==-2 ) { qpdf1 = -qrku;                Cp2 = c2u; Cp3 = -c3u; }
    else if ( pdg_iq==-3 &&  sea_iq && pdg_fq==-3 ) { qpdf1 = -3;                   Cp2 = c2d; Cp3 = -c3d; }
    else if ( pdg_iq==-4 &&  sea_iq && pdg_fq==-4 ) { qpdf1 = -4;                   Cp2 = c2u; Cp3 = -c3u; }
    else if ( pdg_iq==-5 &&  sea_iq && pdg_fq==-5 ) { qpdf1 = -5;                   Cp2 = c2d; Cp3 = -c3d; }
  }   

  // open file in which SF will be stored
  std::ofstream sf_stream(sfFile.c_str());

  // loop over 3 different SF: F1,F2,F3
  for(int sf = 1; sf < 4; sf++) {
    for ( unsigned int i=0; i<sf_q2_array.size(); i++ ) {
      double Q2 = sf_q2_array[i];
      for ( unsigned int j=0; j<sf_x_array.size(); j++ ) {
        double x = sf_x_array[j];

        double z = x; // this variable is introduce in case you want to apply scaling

        // W threshold
        if      (fSF.QrkThrs==1) { 
            if ( Q2*(1/z-1)+mass_nucl*mass_nucl <= TMath::Power(mass_nucl+mPDFQrk[TMath::Abs(pdg_fq)],2) ) { sf_stream << 0. << "  "; continue; }
        } 
        // W threshold and slow rescaling
        else if (fSF.QrkThrs==2) {
            if ( Q2*(1/z-1)+mass_nucl*mass_nucl <= TMath::Power(mass_nucl+mPDFQrk[TMath::Abs(pdg_fq)],2) ) { sf_stream << 0. << "  "; continue; }
            z *= 1+mPDFQrk[TMath::Abs(pdg_fq)]*mPDFQrk[TMath::Abs(pdg_fq)]/Q2;
        }
        // Slow rescaling
        else if (fSF.QrkThrs==3) {
            z *= 1+mPDFQrk[TMath::Abs(pdg_fq)]*mPDFQrk[TMath::Abs(pdg_fq)]/Q2;
        }

        // Fill x,Q2 used to extract PDF. If values outside boundaries then freeze them.
        double xPDF = TMath::Max( z, xPDFmin );
        double Q2PDF = TMath::Max( Q2, Q2PDFmin );
        Q2PDF = TMath::Min( Q2PDF, Q2PDFmax  );

        // Extract PDF requiring then to be higher than zero
#ifdef __GENIE_LHAPDF6_ENABLED__
        double fPDF = fmax( pdf->xfxQ2(qpdf1, xPDF, Q2PDF)/z , 0.);
        if (qpdf2!= -999) fPDF -= fmax( pdf->xfxQ2(qpdf2, xPDF, Q2PDF)/z , 0.);
#endif
#ifdef __GENIE_LHAPDF5_ENABLED__
        double fPDF = fmax( LHAPDF::xfx(xPDF, TMath::Sqrt(Q2PDF), qpdf1)/z , 0.);
        if (qpdf2!= -999) fPDF -= fmax( LHAPDF::xfx(xPDF, TMath::Sqrt(Q2PDF), qpdf2)/z , 0.);
#endif

        // Compute SF
        double tmp = -999;
        if      ( sf==1 ) tmp = fPDF*Cp2/2;
        else if ( sf==2 ) tmp = fPDF*Cp2*z;
        else if ( sf==3 ) tmp = fPDF*Cp3*sign3;

        // Save SF for particular x and Q2 in file
        LOG("HEDISStrucFunc", pDEBUG) << "QrkSFLO" << sf << "[x=" << x << "," << Q2 << "] = " << tmp;
        sf_stream << tmp << "  ";
        
      }
    }
  }

  // Close file in which SF are stored
  sf_stream.close();

}
#ifdef __GENIE_APFEL_ENABLED__
//____________________________________________________________________________
void HEDISStrucFunc::CreateNucSF( const Interaction * in, string sfFile )
{

  // variables used to tag the SF for particular channel
  bool iscc = in->ProcInfo().IsWeakCC();
  bool isnu = pdg::IsNeutrino(in->InitState().ProbePdg());
  bool ispr = pdg::IsProton(in->InitState().Tgt().HitNucPdg());

  // Define the channel that is used in APFEL
  if ( isnu ) APFEL::SetProjectileDIS("neutrino");
  else        APFEL::SetProjectileDIS("antineutrino");
  if ( iscc ) APFEL::SetProcessDIS("CC");
  else        APFEL::SetProcessDIS("NC");
  if ( ispr ) APFEL::SetTargetDIS("proton");
  else        APFEL::SetTargetDIS("neutron");

  APFEL::InitializeAPFEL_DIS();

  // Using APFEL format to store the SF grid
  int nx  = sf_x_array.size();
  int nq2 = sf_q2_array.size();
  double xlist[nx*nq2];
  double q2list[nx*nq2];
  double F2list[nx*nq2];
  double FLlist[nx*nq2];
  double xF3list[nx*nq2];

  int nlist = 0;
  for ( unsigned int i=0; i<sf_q2_array.size(); i++ ) {
    double Q2 = sf_q2_array[i];
    double Q  = TMath::Sqrt(Q2);
    // SF from APFEL are multiplied by a prefactor in NC. We dont want that prefactor
    double norm = iscc ? 1. : 2./TMath::Power( Q2/(Q2 + TMath::Power(APFEL::GetZMass(),2))/4/APFEL::GetSin2ThetaW()/(1-APFEL::GetSin2ThetaW()), 2 );
    APFEL::SetAlphaQCDRef(pdf->alphasQ(Q),Q);
    APFEL::ComputeStructureFunctionsAPFEL(Q,Q);
    for ( unsigned int j=0; j<sf_x_array.size(); j++ ) {
      double x = sf_x_array[j];
      q2list[nlist]  = Q2;
      xlist[nlist]   = x;
      FLlist[nlist]  = norm*APFEL::FLtotal(x);
      F2list[nlist]  = norm*APFEL::F2total(x);
      xF3list[nlist] = norm*APFEL::F3total(x);
      nlist++;
    }
  }

  // open file in which SF will be stored
  std::ofstream sf_stream(sfFile.c_str());

  double sign3 = isnu ? +1. : -1.;  // sign change for nu/nubar in F3
  // loop over 3 different SF: F1,F2,F3
  for(int sf = 1; sf < 4; sf++) {
    for (int i=0; i<nx*nq2; i++) {
      double tmp = 0;
      if      ( sf==1 ) tmp = (F2list[i]-FLlist[i])/2/xlist[i];
      else if ( sf==2 ) tmp = F2list[i];
      else if ( sf==3 ) tmp = sign3 * xF3list[i] / xlist[i];
      // Save SF for particular x and Q2 in file
      LOG("HEDISStrucFunc", pDEBUG) << "NucSFNLO" << sf << "[x=" << xlist[i] << "," << q2list[i] << "] = " << tmp;
      sf_stream << tmp << "  ";
    }
  }
    
  // Close file in which SF are stored
  sf_stream.close();

}
#endif
//____________________________________________________________________________
string HEDISStrucFunc::QrkSFName( const Interaction * in) 
{
  string sin = pdg::IsNeutrino(in->InitState().ProbePdg()) ? "nu_" : "nubar_";
  sin += in->ProcInfo().IsWeakCC() ? "cc_" : "nc_";
  sin += pdg::IsProton(in->InitState().Tgt().HitNucPdg()) ? "p_" : "n_";
  sin += "iq"+to_string(in->InitState().Tgt().HitQrkPdg());
  sin += in->InitState().Tgt().HitSeaQrk() ? "sea_" : "val_";
  sin += "fq"+to_string(in->ExclTag().FinalQuarkPdg());
  return sin;
}
//____________________________________________________________________________
string HEDISStrucFunc::NucSFName( const Interaction * in) 
{
  string sin = pdg::IsNeutrino(in->InitState().ProbePdg()) ? "nu_" : "nubar_";
  sin += in->ProcInfo().IsWeakCC() ? "cc_" : "nc_";
  sin += pdg::IsProton(in->InitState().Tgt().HitNucPdg()) ? "p" : "n";
  return sin;
}
//____________________________________________________________________________
int HEDISStrucFunc::QrkSFCode( const Interaction * in) 
{
  int code = 10000000*pdg::IsNeutrino(in->InitState().ProbePdg());
  code    += 1000000*in->ProcInfo().IsWeakCC();
  code    += 100000*pdg::IsProton(in->InitState().Tgt().HitNucPdg());
  code    += 10000*in->InitState().Tgt().HitSeaQrk();
  code    += 100*(6+in->InitState().Tgt().HitQrkPdg());
  code    += 1*(6+in->ExclTag().FinalQuarkPdg());
  return code;
}
//____________________________________________________________________________
int HEDISStrucFunc::NucSFCode( const Interaction * in) 
{
  int code = 100*pdg::IsNeutrino(in->InitState().ProbePdg());
  code    += 10*in->ProcInfo().IsWeakCC();
  code    += 1*pdg::IsProton(in->InitState().Tgt().HitNucPdg());
  return code;
}
//____________________________________________________________________________
SF_xQ2 HEDISStrucFunc::EvalQrkSFLO( const Interaction * in, double x, double Q2 ) 
{
  int code = QrkSFCode(in);
  SF_xQ2 sf;
  sf.F1 = fQrkSFLOTables[code].Table[kSFT1]->Evaluate(Q2,x);
  sf.F2 = fQrkSFLOTables[code].Table[kSFT2]->Evaluate(Q2,x);
  sf.F3 = fQrkSFLOTables[code].Table[kSFT3]->Evaluate(Q2,x);
  return sf;
}
//____________________________________________________________________________
SF_xQ2 HEDISStrucFunc::EvalNucSFLO( const Interaction * in, double x, double Q2 ) 
{
  int code = NucSFCode(in);
  SF_xQ2 sf;
  sf.F1 = fNucSFLOTables[code].Table[kSFT1]->Evaluate(Q2,x);
  sf.F2 = fNucSFLOTables[code].Table[kSFT2]->Evaluate(Q2,x);
  sf.F3 = fNucSFLOTables[code].Table[kSFT3]->Evaluate(Q2,x);
  return sf;
}
//____________________________________________________________________________
SF_xQ2 HEDISStrucFunc::EvalNucSFNLO( const Interaction * in, double x, double Q2 ) 
{
  int code = NucSFCode(in);
  SF_xQ2 sf;
  sf.F1 = fNucSFNLOTables[code].Table[kSFT1]->Evaluate(Q2,x);
  sf.F2 = fNucSFNLOTables[code].Table[kSFT2]->Evaluate(Q2,x);
  sf.F3 = fNucSFNLOTables[code].Table[kSFT3]->Evaluate(Q2,x);
  return sf;
}

