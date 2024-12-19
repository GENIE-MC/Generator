//____________________________________________________________________________
/*!
\program gmkspl
\brief   GENIE utility program building Structure Functions needed for HEDIS
         package.
         Syntax :
           gmkphotonsf [-h]
         Note :
           [] marks optional arguments.
           <> marks a list of arguments out of which only one can be
              selected at any given time.
         Options :
        ***  See the User Manual for more details and examples. ***
\author  Alfonso Garcia <aagarciasoto \at km3net.de>
         Harvard University & IFIC
\created Dev 8, 2020
\cpright Copyright (c) 2003-2025, The GENIE Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#include "Framework/Messenger/Messenger.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"

#include <TSystem.h>
#include <TMath.h>
#include <Math/IFunction.h>
#include <Math/Integrator.h>

#include <string>
#include <iostream>
#include <fstream>

#ifdef __GENIE_APFEL_ENABLED__
#include "APFEL/APFEL.h"
#endif
#include "LHAPDF/LHAPDF.h"
#ifdef __GENIE_LHAPDF6_ENABLED__
const LHAPDF::PDF* pdf_cache;
#endif

using namespace std;

using namespace genie;
using namespace genie::constants;

int fNucPdg = 0;

#ifdef __GENIE_APFEL_ENABLED__
class PhotonConv: public ROOT::Math::IBaseFunctionOneDim
{
  public:
    PhotonConv(double x, double q) : ROOT::Math::IBaseFunctionOneDim(),xmin(x),Qin(q){};
   ~PhotonConv(){};
    ROOT::Math::IBaseFunctionOneDim * Clone  (void)     const { return new PhotonConv(xmin,Qin); }
    unsigned int                      NDim   (void)     const { return 1; }
    double                            DoEval (double y) const { return 2. * ( TMath::Power(xmin/y,2)+TMath::Power( 1.-xmin/y, 2) ) * pdf_cache->xfxQ(22, y, Qin); }
    void                              SetPar (double x, double q) { xmin=x; Qin=q; }

  private:
    double xmin,Qin;
};

extern "C" void externalsetapfellept_(double* x, double* q, int* irep, double* xl, double* xf);


// Requires a global cache of pdf_cache LHAPDF::PDF
void externalsetapfellept_(double* x, double* q, int* irep, double* xl, double* xf){
  if (*x >= 1 || *x < 0) {
    for ( int i=0; i<13; i++ ) xf[i] = 0.;
    for ( int i=0; i<7; i++ )  xl[i] = 0.;
    return;
  }
  else{
    for ( int i=0; i<13; i++ ) {
      xf[i] = pdf_cache->xfxQ(i-6, *x, *q);
      if ( pdg::IsNeutron(fNucPdg) ) {
        if     ( i==4 ) xf[i] = pdf_cache->xfxQ(-1, *x, *q); 
        else if( i==5 ) xf[i] = pdf_cache->xfxQ(-2, *x, *q);
        else if( i==7 ) xf[i] = pdf_cache->xfxQ( 2, *x, *q);
        else if( i==8 ) xf[i] = pdf_cache->xfxQ( 1, *x, *q);
      }
    }
    for ( int i=0; i<7; i++ ) {
      if ( pdg::IsProton (fNucPdg) ) {
        if      (i==0 || i==6) xl[i] = 0.;
        else if (i==3 )        xl[i] = pdf_cache->xfxQ(22, *x, *q);
        else{       
          double mlep;
          if      (i==1 || i==5) mlep = kMuonMass;
          else if (i==2 || i==4) mlep = kElectronMass;
          ROOT::Math::IBaseFunctionOneDim * func = new PhotonConv(*x,*q);
          ROOT::Math::IntegrationOneDim::Type ig_type = ROOT::Math::IntegrationOneDim::kADAPTIVE;
          ROOT::Math::Integrator ig(*func,ig_type,1,0.01,100000);
          double res = ig.Integral(*x,1.);
          xl[i] = APFEL::AlphaQED(*q)/2./kPi*TMath::Log( *q/mlep ) * res;
        }
      }
      else if ( pdg::IsNeutron(fNucPdg) ) xl[i] = 0.0;
    }
  }

  return;
}
#endif

//____________________________________________________________________________
int main(int argc, char ** argv)
{

  string basedir = "";
  if ( gSystem->Getenv("PHOTON_SF_DATA_PATH")==NULL ) basedir = string(gSystem->Getenv("GENIE")) + "/data/evgen/photon-sf";
  else                                                basedir = string(gSystem->Getenv("PHOTON_SF_DATA_PATH"));
  LOG("gmkphotonsf", pERROR) << "Base directory: " << basedir;

  if ( gSystem->AccessPathName( basedir.c_str(), kWritePermission ) ) {
      LOG("gmkphotonsf", pFATAL) << "Base directory doesnt exist or you dont have write permission.";
      LOG("gmkphotonsf", pFATAL) << "Remember!!!";
      LOG("gmkphotonsf", pFATAL) << "Path to base directory is defined with the enviroment variable PHOTON_SF_DATA_PATH.";
      LOG("gmkphotonsf", pFATAL) << "If not defined, default location is $GENIE/data/evgen/photon-sf";
      assert(0);
  }

  const int nx = 1000.;

  int nucs[2] = { kPdgProton, kPdgNeutron };
  int pdgs[6] = { kPdgNuE, kPdgAntiNuE, kPdgNuMu, kPdgAntiNuMu, kPdgNuTau, kPdgAntiNuTau };
  
#ifdef __GENIE_APFEL_ENABLED__
  for (int k=0; k<2; k++) {
    
    fNucPdg = nucs[k];

    // initialising APFEL framework
    LOG("gmkphotonsf", pINFO) << "Initialising APFEL..." ; 
    string pdfset;
    if      (pdg::IsProton (fNucPdg) ) pdfset = "NNPDF31_nnlo_as_0118_luxqed";
    else if (pdg::IsNeutron(fNucPdg) ) pdfset = "NNPDF31_nnlo_as_0118";

    const LHAPDF::PDFSet set(pdfset);
    pdf_cache = LHAPDF::mkPDF(pdfset,0);

    double xPDFmin,QPDFmin,QPDFmax,mc,mb,mt;
    stringstream( set.get_entry("MCharm") )  >> mc;
    stringstream( set.get_entry("MBottom") ) >> mb;
    stringstream( set.get_entry("MTop") )    >> mt;    
    stringstream( set.get_entry("QMin") )    >> QPDFmin;
    stringstream( set.get_entry("QMax") )    >> QPDFmax;  
    stringstream( set.get_entry("XMin") )    >> xPDFmin;
    LOG("gmkphotonsf", pINFO) << "xPDFmin = " << xPDFmin;
    LOG("gmkphotonsf", pINFO) << "QPDFmin = " << QPDFmin;
    LOG("gmkphotonsf", pINFO) << "QPDFmax = " << QPDFmax;
    LOG("gmkphotonsf", pINFO) << "mc = " << mc;
    LOG("gmkphotonsf", pINFO) << "mb = " << mb;
    LOG("gmkphotonsf", pINFO) << "mt = " << mt;

    APFEL::CleanUp();

    APFEL::SetPDFSet(pdfset);
    APFEL::SetReplica(0);
    APFEL::SetPerturbativeOrder(2);
    APFEL::SetQLimits(QPDFmin,QPDFmax);
    APFEL::SetMaxFlavourPDFs(6);
    APFEL::SetMaxFlavourAlpha(6);
    APFEL::SetNumberOfGrids(3);
    APFEL::SetGridParameters(1,100,5,1e-9);
    APFEL::SetGridParameters(2,40,3,1e-1);
    APFEL::SetGridParameters(3,60,3,8e-1);
    APFEL::SetGFermi(kGF);
    APFEL::SetPoleMasses(mc,mb,mt);
    APFEL::SetTheory("QUniD");
    APFEL::EnableLeptonEvolution(true);
    APFEL::SetFastEvolution(true);
    APFEL::SetPDFSet("leptexternal");
    APFEL::InitializeAPFEL();

    LOG("gmkphotonsf", pWARN) << "Init EvolveAPFEL";        
    APFEL::EvolveAPFEL(QPDFmin,kMw);
    LOG("gmkphotonsf", pWARN) << "End EvolveAPFEL";        

    // open file in which SF will be stored

    double x[nx];
    for ( int i=0; i<nx; i++ ) x[i] = TMath::Power( 10, TMath::Log10(xPDFmin) + i*(TMath::Log10(1.)-TMath::Log10(xPDFmin))/(1000.-1) );

    for(int j=0; j<6; j++) {

      string SFname = basedir + "/PhotonSF_hitnuc"+to_string(fNucPdg)+"_hitlep"+to_string(pdgs[j])+".dat";
      std::ofstream sf_stream(SFname);
      for ( int i=0; i<nx; i++ ) {
        double tmp = 0;
        if      ( pdg::IsNuE      (pdgs[j]) ) tmp = APFEL::xLepton( 1,x[i]);
        else if ( pdg::IsAntiNuE  (pdgs[j]) ) tmp = APFEL::xLepton(-1,x[i]);
        else if ( pdg::IsNuMu     (pdgs[j]) ) tmp = APFEL::xLepton( 2,x[i]);
        else if ( pdg::IsAntiNuMu (pdgs[j]) ) tmp = APFEL::xLepton(-2,x[i]);
        else if ( pdg::IsNuTau    (pdgs[j]) ) tmp = APFEL::xLepton( 3,x[i]);
        else if ( pdg::IsAntiNuTau(pdgs[j]) ) tmp = APFEL::xLepton(-3,x[i]);
        LOG("gmkphotonsf", pWARN) << "SF " << pdgs[j] << " [x=" << x[i] << "] = " << tmp;
        sf_stream << x[i] << " " << tmp << endl;
      }
      // Close file in which SF are stored
      sf_stream.close();
    }       
    
  }
#endif

}