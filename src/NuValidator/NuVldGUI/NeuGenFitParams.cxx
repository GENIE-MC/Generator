//_____________________________________________________________________________
/*!

\class    genie::nuvld::NeuGenFitParams

\brief    Holds the NeuGEN parameters that can be fitted and it is set by the
          NeuGenFitParamsDialog. It has a rather flat structure since this
          is more convenient for use with ROOT's fitting machinery. There is
          a correspondence with NeuGenConfig class - see this for more info
          on the fitted parameters.

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>

\created  January 28, 2005
*/
//_____________________________________________________________________________

#include "Messenger/Messenger.h"
#include "NuVldGUI/NeuGenFitParams.h"

using namespace genie::nuvld;

//_____________________________________________________________________________
NeuGenFitParams::NeuGenFitParams ()
{
  this->Init();
}
//_____________________________________________________________________________
NeuGenFitParams::NeuGenFitParams(const NeuGenFitParams & ngfp)
{
  this->Init();

  for (int ip = 0; ip < kNNGFitParams; ip++)
  {
    fIsFitted[ip] = ngfp.fIsFitted[ip];
    fRangeMin[ip] = ngfp.fRangeMin[ip];
    fRangeMax[ip] = ngfp.fRangeMax[ip];
    fStep[ip]     = ngfp.fStep[ip];
  }
}
//_____________________________________________________________________________
NeuGenFitParams::~NeuGenFitParams()
{
  delete [] fIsFitted;
  delete [] fRangeMin;
  delete [] fRangeMax;
  delete [] fStep;
}
//_____________________________________________________________________________
int NeuGenFitParams::NFittedParams(void) const
{
  int nf = 0;

  for (int ip = 0; ip < kNNGFitParams; ip++) if(this->IsFitted(ip)) nf++;

  return nf;
}
//_____________________________________________________________________________
bool NeuGenFitParams::IsFitted(int iparam) const
{
  if ( this->InRange(iparam) ) return fIsFitted[iparam];
  else 
  {
    LOG("NuVld", pERROR) << "Fit param: " << iparam << " out of range";
  }
  return false;
}
//_____________________________________________________________________________
double NeuGenFitParams::RangeMin(int iparam) const
{
  if ( this->InRange(iparam) ) return fRangeMin[iparam];
  else 
  {
    LOG("NuVld", pERROR) << "Fit param: " << iparam << " out of range";
  }
  return 0.;
}
//_____________________________________________________________________________
double NeuGenFitParams::RangeMax(int iparam) const
{
  if ( this->InRange(iparam) ) return fRangeMax[iparam];
  else 
  {
    LOG("NuVld", pERROR) << "Fit param: " << iparam << " out of range";
  }
  return 0.;
}
//_____________________________________________________________________________
double NeuGenFitParams::Step(int iparam) const
{
  if ( this->InRange(iparam) ) return fStep[iparam];
  else 
  {
    LOG("NuVld", pERROR) << "Fit param: " << iparam << " out of range";
  }
  return 0.;
}
//_____________________________________________________________________________
string NeuGenFitParams::ParamAsString(int iparam) const
{
  if ( this->InRange(iparam) ) {

     NeuGenFitParam_t param = (NeuGenFitParam_t) iparam;
  
     switch(param) {

     case ( kNgfMaQel        ): return "Ma-QEL (GeV)"; break;
     case ( kNgfMaRes        ): return "Ma-RES (GeV)"; break;
     case ( kNgfMaCoh        ): return "Ma-COH (GeV)"; break;
     case ( kNgfQelFa0       ): return "QEL-FA(Q2=0)"; break;
     case ( kNgfQelEta       ): return "QEL-ETA     "; break;
     case ( kNgfResOmega     ): return "RES-OMEGA   "; break;
     case ( kNgfResZ         ): return "RES-Z       "; break;
     case ( kNgfCohR0        ): return "COH Nucl. R0"; break;
     case ( kNgfCohREI       ): return "COH pi Re/Im"; break;
     case ( kNgfKnoB         ): return "KNO B       "; break;
     case ( kNgfKnoAvp       ): return "KNO A (v+p) "; break;
     case ( kNgfKnoAvn       ): return "KNO A (v+n) "; break;
     case ( kNgfKnoAvbp      ): return "KNO A (vb+p)"; break;
     case ( kNgfKnoAvbn      ): return "KNO A (vb+n)"; break;
     case ( kNgfDisResM2vp   ): return "D/R (v+p /2)"; break;
     case ( kNgfDisResM3vp   ): return "D/R (v+p /3)"; break;
     case ( kNgfDisResM2vn   ): return "D/R (v+n /2)"; break;
     case ( kNgfDisResM3vn   ): return "D/R (v+n /3)"; break;
     case ( kNgfDisResM2vbp  ): return "D/R (vb+p/2)"; break;
     case ( kNgfDisResM3vbp  ): return "D/R (vb+p/3)"; break;
     case ( kNgfDisResM2vbn  ): return "D/R (vb+n/2)"; break;
     case ( kNgfDisResM3vbn  ): return "D/R (vb+n/3)"; break;
     case ( kNgfNuDisScale   ): return "Neutrino DIS scale factor";      break;
     case ( kNgfNuBarDisScale): return "Anti-Neutrino DIS scale factor"; break;
     }
     
  } else {
    LOG("NuVld", pERROR) << "Fit param: " << iparam << " out of range";
  }
  return "";
}
//_____________________________________________________________________________
bool NeuGenFitParams::InRange(int iparam) const
{
  return (iparam >= 0 && iparam < kNNGFitParams);
}
//_____________________________________________________________________________
void NeuGenFitParams::Init(void)
{
  fIsFitted = new bool   [kNNGFitParams];
  fRangeMin = new double [kNNGFitParams];
  fRangeMax = new double [kNNGFitParams];
  fStep     = new double [kNNGFitParams];

  for (int ip = 0; ip < kNNGFitParams; ip++)
  {
    fIsFitted[ip] = false;
    fRangeMin[ip] = 0.;
    fRangeMax[ip] = 0.;
    fStep[ip]     = 0.01;
  }    
}
//_____________________________________________________________________________
