//_____________________________________________________________________________
/*!

\class    genie::nuvld::GuiFitKernel

\brief

\author   Costas Andreopoulos (Rutherford Lab.)  <costas.andreopoulos \at stfc.ac.uk>

\created  January 25, 2005
*/
//_____________________________________________________________________________

#ifndef _GUI_FIT_KERNEL_H_
#define _GUI_FIT_KERNEL_H_
/*
class TF1;
class TH2F;
class TGraph;
class TGraphAsymmErrors;
class TMinuit;

namespace genie {

class Spline;

namespace nuvld {

class NeuGenFitParams;

class GuiFitKernel {

public:

  GuiFitKernel();
  ~GuiFitKernel();

   //-- configure the GuiFitKernel

   void SetFitParams       (NeuGenFitParams * nfp) { fNGFP = nfp; }
   void SetFitRange        (float xmin, float xmax);
   void SetScaleWithEnergy (bool tf) { fScaleWithE = tf; }

   void Reset       (void);
   void PrintConfig (void);

   //-- various fitters or parameter scanners

   void  DoSimpleFit       (bool fit_norm);
   void  DoFloatingNormFit (void);

   bool  MCParamScanning   (void);
   bool  ChisqScan1D       (void);
   bool  ChisqScan2D       (void);

   //-- get the fit function

   TF1 *    FitFunction         (void) const { return fFunc1d; }
   TGraph * GetResidualsAsGraph (void);

   //tmp
   TGraph * lowb;
   TGraph * highb;
   TGraph * chisq1d;
   TH2F *   chisq2d;

private:

   float CurrFitParamValue (int iparam);
   void  InitFitParameters (void);
   void  InitXSecNormFitParameters (void);
   void  InitMinuitFitParameters(TMinuit * minuit);

   bool  fScaleWithE;
   float fXmin;
   float fXmax;

   TF1 * fFunc1d;
   NeuGenFitParams * fNGFP;

};

} // nuvld namespace
} // genie namespace
*/
#endif

