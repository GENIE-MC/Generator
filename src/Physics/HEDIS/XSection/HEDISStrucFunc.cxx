//____________________________________________________________________________
/*
 Copyright (c) 2003-2018, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Alfonso Garcia <alfonsog \at nikhef.nl>
         NIKHEF

 For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________

#include "Physics/HEDIS/XSection/HEDISStrucFunc.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Utils/RunOpt.h"
#include "Framework/Conventions/Constants.h"

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
HEDISStrucFunc::HEDISStrucFunc(string basedir, SF_info sfinfo)
{

  fSF = sfinfo;

  if (basedir=="") basedir=string(gSystem->Getenv("GENIE")) + "/data/evgen/hedis-sf";
  LOG("HEDISStrucFunc", pERROR) << "Base diretory: " << basedir;

  if ( gSystem->AccessPathName( basedir.c_str(), kWritePermission ) ) {
      LOG("HEDISStrucFunc", pERROR) << "Base diretory doesnt exist or you dont have write permission.";
      assert(0);
  }

  // Name of the directory where SF tables are (or will be) saved. 
  string SFname = basedir + "/" + RunOpt::Instance()->Tune()->Name();

  // Check that the directory where SF tables are stored exists
  LOG("HEDISStrucFunc", pINFO) << "SF are (or will be) in following directory: " << SFname;
  if ( gSystem->mkdir(SFname.c_str())==0 ) {
    LOG("HEDISStrucFunc", pWARN) << "Bad news! Directory doesnt exists.";
    LOG("HEDISStrucFunc", pWARN) << "HEDIS package requires computation of SF";        
    LOG("HEDISStrucFunc", pWARN) << "This will be SLOW!!!!!";        
    LOG("HEDISStrucFunc", pINFO) << "Creating Metafile with input information";
    std::ofstream meta_stream((SFname+"/Inputs.txt").c_str());
    meta_stream << fSF;
    meta_stream.close();
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
      LOG("HEDISStrucFunc", pERROR) << "Info from MetaFile and Tune doesnt match";
      LOG("HEDISStrucFunc", pERROR) << "MetaFile Path : " << SFname << "/Inputs.txt";        
      LOG("HEDISStrucFunc", pERROR) << "From Tune : ";        
      LOG("HEDISStrucFunc", pERROR) << cm;        
      LOG("HEDISStrucFunc", pERROR) << "From MetaFile : ";        
      LOG("HEDISStrucFunc", pERROR) << fSF;        
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

  // Load structure functions for each quark at LO
  for( int ch=1; ch<kHEDISQrk_numofchannels; ch++ ) {
    string sfFile = SFname + "/QrkSF_LO_" + HEDISChannel::AsString((HEDISQrkChannel_t)ch) + ".dat";
    // Make sure data files are available
    LOG("HEDISStrucFunc", pINFO) << "Checking if file " << sfFile << " exists...";        
    if ( gSystem->AccessPathName( sfFile.c_str()) ) {
      LOG("HEDISStrucFunc", pWARN) << "File doesnt exist. SF table will be computed.";        
      CreateQrkSF( (HEDISQrkChannel_t)ch, sfFile );
    }
    else if ( atoi(gSystem->GetFromPipe(("wc -w "+sfFile+" | awk '{print $1}'").c_str()))!=kSFT3*nx*ny ) {
      LOG("HEDISStrucFunc", pWARN) << "File does not contain all the need points. SF table will be recomputed.";        
      gSystem->Exec(("rm "+sfFile).c_str());
      CreateQrkSF( (HEDISQrkChannel_t)ch, sfFile );
    }
    std::ifstream sf_stream(sfFile.c_str(), std::ios::in);
    // Loop over F1,F2,F3
    for(int sf = 1; sf <= kSFnumber; ++sf) {
      // Loop over x/Q2 bins
      for ( int ij=0; ij<nx*ny; ij++ ) sf_stream >> z[ij];
      // Create SF tables with BLI2DNonUnifGrid using x,Q2 binning
      fQrkSFLOTables[(HEDISQrkChannel_t)ch].Table[(HEDISStrucFuncType_t)sf] = new genie::BLI2DNonUnifGrid( nx, ny, x, y, z );
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
    //compute structure functions for each nucleon
    for( int ch=1; ch<kHEDISNuc_numofchannels; ch++ ) {
      // Load structure functions for each nucleon at NLO
      string sfFile = SFname + "/NucSF_NLO_" + HEDISChannel::AsString((HEDISNucChannel_t)ch) + ".dat";
      // Make sure data files are available
      LOG("HEDISStrucFunc", pINFO) << "Checking if file " << sfFile << " exists...";        
      if ( gSystem->AccessPathName( sfFile.c_str()) ) {
#ifdef __GENIE_APFEL_ENABLED__
        LOG("HEDISStrucFunc", pWARN) << "File doesnt exist. SF table will be computed.";        
        CreateNucSF( (HEDISNucChannel_t)ch, sfFile );
#else
        LOG("HEDISStrucFunc", pERROR) << "File doesnt exist. APFEL is needed for NLO SF";        
        assert(0);
#endif
      }
      else if ( atoi(gSystem->GetFromPipe(("wc -w "+sfFile+" | awk '{print $1}'").c_str()))!=kSFT3*nx*ny ) {
#ifdef __GENIE_APFEL_ENABLED__
        LOG("HEDISStrucFunc", pWARN) << "File does not contain all the need points. SF table will be recomputed.";        
        gSystem->Exec(("rm "+sfFile).c_str());
        CreateQrkSF( (HEDISQrkChannel_t)ch, sfFile );
#else
        LOG("HEDISStrucFunc", pERROR) << "File does not contain all the need points. APFEL is needed for NLO SF";        
        assert(0);
#endif
      }
      std::ifstream sf_stream(sfFile.c_str(), std::ios::in);
      // Loop over F1,F2,F3
      for(int sf = 1; sf <= kSFnumber; ++sf) {
        // Loop over x/Q2 bins
        for ( int ij=0; ij<nx*ny; ij++ ) sf_stream >> z[ij];
        // Create SF tables with BLI2DNonUnifGrid using x,Q2 binning
        fNucSFNLOTables[(HEDISNucChannel_t)ch].Table[(HEDISStrucFuncType_t)sf] = new genie::BLI2DNonUnifGrid( nx, ny, x, y, z );
      }        
      //compute structure functions for each nucleon at LO using quark grids
      LOG("HEDISStrucFunc", pDEBUG) << "Creating LO " << sfFile;        
      int frstqrkch = HEDISChannel::GetFirstHEDISQrkChannel((HEDISNucChannel_t)ch);
      int lastqrkch = HEDISChannel::GetLastHEDISQrkChannel((HEDISNucChannel_t)ch);
      // Loop over F1,F2,F3
      for(int sf = 1; sf <= kSFnumber; ++sf) {
        int ij = 0;
        // Loop over Q2 bins
        for (int i=0; i<nx; i++) {
          // Loop over x bins
          for (int j=0; j<ny; j++) {
            double sum = 0.;
            // NucSF = sum_qrks QrkSF
            for( int qch=frstqrkch; qch<=lastqrkch; qch++ ) sum += fQrkSFLOTables[(HEDISQrkChannel_t)qch].Table[(HEDISStrucFuncType_t)sf]->Evaluate(x[i],y[j]);
            z[ij] = sum;
            ij++;
          }
        }
        // Create SF tables with BLI2DNonUnifGrid using x,Q2 binning
        fNucSFLOTables[(HEDISNucChannel_t)ch].Table[(HEDISStrucFuncType_t)sf] = new genie::BLI2DNonUnifGrid( nx, ny, x, y, z );
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
HEDISStrucFunc * HEDISStrucFunc::Instance(string basedir, SF_info sfinfo)
{
  if(fgInstance == 0) {
    LOG("HEDISStrucFunc", pINFO) << "Late initialization";
    static HEDISStrucFunc::Cleaner cleaner;
    cleaner.DummyMethodAndSilentCompiler();
    fgInstance = new HEDISStrucFunc(basedir,sfinfo);
  }  
  return fgInstance;
}


//____________________________________________________________________________
void HEDISStrucFunc::CreateQrkSF( HEDISQrkChannel_t ch, string sfFile ) 
{

  // variables used to tag the SF for particular channel
  int pdg_nucl  = HEDISChannel::HitNuclPdg(ch);   // PDG of target nucleon
  bool sea_iq   = HEDISChannel::HitQuarkSea(ch);  // Valence/sea hit quark
  int pdg_iq    = HEDISChannel::HitQuarkPdg(ch);  // PDG hit quark
  int pdg_fq    = HEDISChannel::FnlQuarkPdg(ch);  // PDG produced quark

  // variables used for certain quark production threshold options
  double mass_fq   = mPDFQrk[TMath::Abs(pdg_fq)];                     // mass produced quark
  double mass_nucl = (pdg_nucl==2212) ? kProtonMass : kNeutronMass;   // mass target nucleon

  // up and down quark swicth depending on proton or neutron interaction  
  int qrkd = 0;
  int qrku = 0;
  if      ( pdg_nucl==2212 ) { qrkd = 1 ; qrku = 2; }
  else if ( pdg_nucl==2112 ) { qrkd = 2 ; qrku = 1; }

  // variables associated to the PDF and coupling of the quarks
  int qpdf1 = -999;                                     // number used to compute PDF
  int qpdf2 = -999;                                     // auxiliary number used to compute PDF for valence quarks
  double Cp2 = -999;                                    // couping for F1,F2
  double Cp3 = -999;                                    // couping for F3
  double sign3 = (HEDISChannel::IsNu(ch)) ? +1. : -1.;  // sign change for nu/nubar in F3
  if ( HEDISChannel::InteractionType(ch) == kIntWeakCC ) {
    if ( HEDISChannel::IsNu(ch) ) {
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
  else if ( HEDISChannel::InteractionType(ch) == kIntWeakNC ) {           
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
            if ( Q2*(1/z-1)+mass_nucl*mass_nucl <= TMath::Power(mass_nucl+mass_fq,2) ) { sf_stream << 0. << "  "; continue; }
        } 
        // W threshold and slow rescaling
        else if (fSF.QrkThrs==2) {
            if ( Q2*(1/z-1)+mass_nucl*mass_nucl <= TMath::Power(mass_nucl+mass_fq,2) ) { sf_stream << 0. << "  "; continue; }
            z *= 1+mass_fq*mass_fq/Q2;
        }
        // Slow rescaling
        else if (fSF.QrkThrs==3) {
            z *= 1+mass_fq*mass_fq/Q2;
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
void HEDISStrucFunc::CreateNucSF( HEDISNucChannel_t ch, string sfFile )
{

  // Define the channel that is used in APFEL
  if ( HEDISChannel::IsNu(ch) ) APFEL::SetProjectileDIS("neutrino");
  else                          APFEL::SetProjectileDIS("antineutrino");
  if      ( HEDISChannel::InteractionType(ch) == kIntWeakCC ) APFEL::SetProcessDIS("CC");
  else if ( HEDISChannel::InteractionType(ch) == kIntWeakNC ) APFEL::SetProcessDIS("NC");
  if      ( HEDISChannel::HitNuclPdg(ch)==2212 ) APFEL::SetTargetDIS("proton");
  else if ( HEDISChannel::HitNuclPdg(ch)==2112 ) APFEL::SetTargetDIS("neutron");

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
    double norm = (HEDISChannel::InteractionType(ch)==kIntWeakCC) ? 1. : 2./TMath::Power( Q2/(Q2 + TMath::Power(APFEL::GetZMass(),2))/4/APFEL::GetSin2ThetaW()/(1-APFEL::GetSin2ThetaW()), 2 );
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

  double sign3 = (HEDISChannel::IsNu(ch)) ? +1. : -1.;  // sign change for nu/nubar in F3
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
SF_xQ2 HEDISStrucFunc::EvalQrkSFLO( HEDISQrkChannel_t ch, double x, double Q2 ) 
{
  SF_xQ2 sf;
  sf.F1 = fQrkSFLOTables[ch].Table[kSFT1]->Evaluate(Q2,x);
  sf.F2 = fQrkSFLOTables[ch].Table[kSFT2]->Evaluate(Q2,x);
  sf.F3 = fQrkSFLOTables[ch].Table[kSFT3]->Evaluate(Q2,x);
  return sf;
}
//____________________________________________________________________________
SF_xQ2 HEDISStrucFunc::EvalNucSFLO( HEDISNucChannel_t ch, double x, double Q2 ) 
{
  SF_xQ2 sf;
  sf.F1 = fNucSFLOTables[ch].Table[kSFT1]->Evaluate(Q2,x);
  sf.F2 = fNucSFLOTables[ch].Table[kSFT2]->Evaluate(Q2,x);
  sf.F3 = fNucSFLOTables[ch].Table[kSFT3]->Evaluate(Q2,x);
  return sf;
}
//____________________________________________________________________________
SF_xQ2 HEDISStrucFunc::EvalNucSFNLO( HEDISNucChannel_t ch, double x, double Q2 ) 
{
  SF_xQ2 sf;
  sf.F1 = fNucSFNLOTables[ch].Table[kSFT1]->Evaluate(Q2,x);
  sf.F2 = fNucSFNLOTables[ch].Table[kSFT2]->Evaluate(Q2,x);
  sf.F3 = fNucSFNLOTables[ch].Table[kSFT3]->Evaluate(Q2,x);
  return sf;
}