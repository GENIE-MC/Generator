//____________________________________________________________________________
/*!

\class    genie::FKR

\brief    Rein-Seghal package utility class for computing and holding the
          Feynmann-Kislinger-Ravndall (FKR) baryon excitation model parameters.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 03, 2004

*/
//____________________________________________________________________________

#include "BaryonResonance/BaryonResUtils.h"
#include "Conventions/Constants.h"
#include "ReinSeghal/FKR.h"
#include "Messenger/Messenger.h"

using std::endl;

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
namespace genie
{
  ostream & operator << (ostream & stream, const FKR & parameters)
  {
     parameters.Print(stream);
     return stream;
  }
}
//____________________________________________________________________________
FKR::FKR()
{
  this->Initialize();
}
//____________________________________________________________________________
FKR::~FKR()
{

}
//____________________________________________________________________________
void FKR::Calculate(double q2, double W, double mN, int n)
{
// q2   : momentum transfer < 0
// W    : hadronic invariant mass
// nqres: resonance index

  double mN2 = TMath::Power(mN, 2);
  double W2  = TMath::Power(W,  2);
  double k   = 0.5 * (W2 - mN2) / mN;
  double v   = k - 0.5*q2 / mN;
  double Q2  = TMath::Power(v, 2) - q2; // note: in RS Q2 does not denote -q2
  double Q   = TMath::Sqrt(TMath::Abs(Q2));

  LOG("FKR", pDEBUG)
         << ENDL
         <<  "q2 = " << q2 << ENDL
         <<  "W  = " << W  << ENDL
         <<  "mN = " << mN << ENDL
         <<  "k  = " << k  << ENDL
         <<  "v  = " << v  << ENDL
         <<  "Q2 = " << Q2 << ENDL
         <<  "Q  = " << Q  << ENDL;

  double omega = fOmega;
  double zeta  = fZeta;
  double ma2   = fMa2;
  double mv2   = fMv2;
  double d     = TMath::Power(W+mN,2.) - q2;
  double osc   = TMath::Power(1 - 0.25 * q2/mN2, 0.5-n);
  double GV    = osc * TMath::Power( 1./(1-q2/mv2), 2);
  double GA    = osc * TMath::Power( 1./(1-q2/ma2), 2);

  fLamda = TMath::Sqrt(2./omega) * mN * Q / W;
  fTv    = TMath::Sqrt(omega/2.) * GV / (3*W);
  fRv    = kSqrt2 * (mN/W)*(W+mN)*Q*GV / d;
  fS     = (-q2/Q2) * (3*W*mN + q2 - mN2) * GV / (6*mN2);
  fTa    = (2./3.) *zeta * TMath::Sqrt(omega/2.) * (mN/W) * Q * GA / d;
  fRa    = (kSqrt2/6.) * zeta * (GA/W) * ( (W+mN) + 2*n*omega*W/ d );
  fB     = (1./3.) * (zeta/W) * TMath::Sqrt(omega/2.) * ( 1 + (W2-mN2+q2) / d ) * GA;
  fC     = (1./6.) * (zeta/Q) * ( W2 - mN2 + n*omega*(W2-mN2+q2) /d) * (GA/mN);

  //-- compute frequently used combinations of FKR params

  fR      = fRv;
  fRplus  = - (fRv + fRa);
  fRminus = - (fRv - fRa);
  fT      = fTv;
  fTplus  = - (fTv + fTa);
  fTminus = - (fTv - fTa);

  double w = - kSin8w2; // 8w = theta-weinberg

  fLamdaRminus  = fLamda  * fRminus;
  fLamdaRplus   = fLamda  * fRplus;
  fLamdaTminus  = fLamda  * fTminus;
  fLamdaTplus   = fLamda  * fTplus;
  fLC           = fLamda  * fC;
  fLS           = fLamda  * fS;
  fLR           = fLamda  * fR;
  fLT           = fLamda  * fT;
  fLC_2B        = fLamda  * fC - 2 * fB;
  fLC_3B        = fLamda  * fC - 3 * fB;
  fLC_5B        = fLamda  * fC - 5 * fB;
  fLamda2       = fLamda  * fLamda;
  fLamda2Rminus = fLamda2 * fRminus;
  fLamda2Rplus  = fLamda2 * fRplus;
  fLamda2C      = fLamda2 * fC;
  fLamda2S      = fLamda2 * fS;
  fLamda2R      = fLamda2 * fR;
  fRminus_wR    = fRminus - 1 * w * fRv;
  fRminus_2wR   = fRminus - 2 * w * fRv;
  fRminus_3wR   = fRminus - 3 * w * fRv;
  fRminus_4wR   = fRminus - 4 * w * fRv;
  fRplus_wR     = fRplus  - 1 * w * fRv;
  fRplus_2wR    = fRplus  - 2 * w * fRv;
  fRplus_3wR    = fRplus  - 3 * w * fRv;
  fRplus_4wR    = fRplus  - 4 * w * fRv;
  fTminus_wTv   = fTminus - 1 * w * fTv;
  fTminus_2wTv  = fTminus - 2 * w * fTv;
  fTminus_3wTv  = fTminus - 3 * w * fTv;
  fTminus_4wTv  = fTminus - 4 * w * fTv;
  fTplus_wTv    = fTplus  - 1 * w * fTv;
  fTplus_2wTv   = fTplus  - 2 * w * fTv;
  fTplus_3wTv   = fTplus  - 3 * w * fTv;
  fTplus_4wTv   = fTplus  - 4 * w * fTv;
}
//____________________________________________________________________________
void FKR::Print(ostream & stream) const
{
  stream << endl;
  stream << " lamda = " << fLamda   << endl;
  stream << " Tv    = " << fTv      << endl;
  stream << " Rv    = " << fRv      << endl;
  stream << " S     = " << fS       << endl;
  stream << " Ta    = " << fTa      << endl;
  stream << " Ra    = " << fRa      << endl;
  stream << " B     = " << fB       << endl;
  stream << " C     = " << fC       << endl;
  stream << " R     = " << fR       << endl;
  stream << " T     = " << fT       << endl;
  stream << " T+    = " << fTplus   << endl;
  stream << " T-    = " << fTminus  << endl;
  stream << " R+    = " << fRplus   << endl;
  stream << " R-    = " << fRminus  << endl;
}
//____________________________________________________________________________
void FKR::Initialize(void)
{
  fZeta      = 0;
  fOmega     = 0;
  fMa2       = 0;
  fMv2       = 0;

  //-- initialize FKR options

  fLamda        =  0.0;
  fTv           =  0.0;
  fRv           =  0.0;
  fS            =  0.0;
  fTa           =  0.0;
  fRa           =  0.0;
  fB            =  0.0;
  fC            =  0.0;
  fR            =  0.0;
  fT            =  0.0;
  fTplus        =  0.0;
  fTminus       =  0.0;
  fRplus        =  0.0;
  fRminus       =  0.0;

  fLamdaRminus  = 0.0;
  fLamdaRplus   = 0.0;
  fLamdaTminus  = 0.0;
  fLamdaTplus   = 0.0;
  fLC           = 0.0;
  fLS           = 0.0;
  fLR           = 0.0;
  fLT           = 0.0;
  fLC_2B        = 0.0;
  fLC_3B        = 0.0;
  fLC_5B        = 0.0;
  fLamda2       = 0.0;
  fLamda2Rminus = 0.0;
  fLamda2Rplus  = 0.0;
  fLamda2C      = 0.0;
  fLamda2S      = 0.0;
  fLamda2R      = 0.0;
  fRminus_wR    = 0.0;
  fRminus_2wR   = 0.0;
  fRminus_3wR   = 0.0;
  fRminus_4wR   = 0.0;
  fRplus_wR     = 0.0;
  fRplus_2wR    = 0.0;
  fRplus_3wR    = 0.0;
  fRplus_4wR    = 0.0;
  fTminus_wTv   = 0.0;
  fTminus_2wTv  = 0.0;
  fTminus_3wTv  = 0.0;
  fTminus_4wTv  = 0.0;
  fTplus_wTv    = 0.0;
  fTplus_2wTv   = 0.0;
  fTplus_3wTv   = 0.0;
  fTplus_4wTv   = 0.0;
}
//____________________________________________________________________________
void FKR::SetZeta(double zeta )
{
  fZeta = zeta;
}
//____________________________________________________________________________
void FKR::SetOmega(double omega)
{
  fOmega = omega;
}
//____________________________________________________________________________
void FKR::SetMa2(double ma2)
{
  fMa2 = ma2;
}
//____________________________________________________________________________
void FKR::SetMv2(double mv2)
{
  fMv2 = mv2;
}
//____________________________________________________________________________

