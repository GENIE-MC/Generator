// --------------------------------------------------------------
// ROOT macro to load hadronization model tuning data to ntuples.
// See the data files for description of the ntuple fields.
//
// C.Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
// --------------------------------------------------------------
{
  //
  // Define input files and formats
  //

  TString genie    = gSystem->Getenv("GENIE");
  TString basedir  = genie + "/data/hadronization/";

  // -- backward - forward multiplicity correlation data
  TTree   bfc;
  TString bfcFile = basedir + "./hBkwFwdCorrelation.data";
  TString bfcFmt  = "nu/I:nuc/I:Wmin/F:Wmax/F:nF/F:nB/F:dnB/F";

  // -- multiplicity vs hadronic invariant mass data
  TTree   mult;
  TString multFile = basedir + "./hMultiplicityVsW2.data";
  TString multFmt  = "nu/I:nuc/I:tgt/I:W2/F:n/F:dn/F:nt/I:xfh/I";

  // -- pi0 - charged hadron multiplicity correlation data
  TTree   pi0c;
  TString pi0cFile = basedir + "./hPi0HCorrelation.data";
  TString pi0cFmt  = "nu/I:nuc/I:Wmin/F:Wmax/F:n/I:npi0/F:dnpi0/F:nt/I";

  // -- dispersion vs multiplicity data
  TTree   dn;
  TString dnFile = basedir + "./hDispersionVsMult.data";
  TString dnFmt  = "nu/I:nuc/I:tgt/I:n/F:D/F:dD/F:nt/I";

  // -- disperison over average multiplicity vs hadronic invariant mass
  TTree   dnw;
  TString dnwFile = basedir + "./hRDispersionVsW2.data";
  TString dnwFmt  = "nu/I:nuc/I:tgt/I:W2/F:DN/F:dDN/F:nt/I";

  // -- corrected normalized invariant xf distributions
  TTree   cnixf;
  TString cnixfFile = basedir + "./hCorrNormInvXFDistr.data";
  TString cnixfFmt  = "nu/I:nuc/I:Wmin/F:xBmin/F:xF/F:F/F:dF/F:nt/I";

  // -- corrected normalized xf distributions
  TTree   cnxf;
  TString cnxfFile = basedir + "./hCorrNormXFDistr.data";
  TString cnxfFmt  = "nu/I:nuc/I:Wmin/F:xBmin/F:xF/F:F/F:dF/F:nt/I";

  // -- normalized topological cross sections
  TTree   ntxs;
  TString ntxsFile = basedir + "./hNormTopologicalXSec.data";
  TString ntxsFmt  = "nu/I:nuc/I:n/I:W2/F:P/F:Pp/F:Pm";

  // -- nucleon xF
  TTree   nucxf;
  TString nucxfFile = basedir + "./hNucleonXF.data";
  TString nucxfFmt  = "xF/F:P/F:dxFm/F:dxFp/F:dPm/F:dPp/F";

  // -- nucleon pT2
  TTree   nucpt;
  TString nucptFile = basedir + "./hNucleonPT2.data";
  TString nucptFmt  = "pT2/F:P/F:dP/F";

  //
  // Load all ntuples
  //

  bfc.  ReadFile( bfcFile.  Data(), bfcFmt.  Data() );
  mult. ReadFile( multFile. Data(), multFmt. Data() );
  pi0c. ReadFile( pi0cFile. Data(), pi0cFmt. Data() );
  dn.   ReadFile( dnFile.   Data(), dnFmt.   Data() );
  dnw.  ReadFile( dnwFile.  Data(), dnwFmt.  Data() );
  cnixf.ReadFile( cnixfFile.Data(), cnixfFmt.Data() );
  cnxf. ReadFile( cnxfFile. Data(), cnxfFmt. Data() );
  ntxs. ReadFile( ntxsFile. Data(), ntxsFmt. Data() );
  nucxf.ReadFile( nucxfFile.Data(), nucxfFmt.Data() );
  nucpt.ReadFile( nucptFile.Data(), nucptFmt.Data() );
}

