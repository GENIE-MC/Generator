//____________________________________________________________________________
/*
 Copyright (c) 2003-2022, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Authors: Igor Kakorin <kakorin@jinr.ru>, Joint Institute for Nuclear Research
          adapted from fortran code provided by Toru Sato <tsato@rcnp.osaka-u.ac.jp>

 For the class documentation see the corresponding header file.
*/
//____________________________________________________________________________

#include <TSystem.h>
#include <TMath.h>
#include <Math/SpecFuncMathMore.h>

#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/RefFrame.h"
#include "Framework/Utils/Range1.h"
#include "Framework/Conventions/KineVar.h"
#include "Framework/Conventions/Units.h"
#include <Framework/Conventions/KinePhaseSpace.h>
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/Utils/KineUtils.h"
#include "Physics/XSectionIntegration/XSecIntegratorI.h"
#include "Physics/Resonance/XSection/DCCSPPPXSec.h"
#include "Physics/NuclearState/FermiMomentumTablePool.h"
#include "Physics/NuclearState/FermiMomentumTable.h"
#include "Physics/NuclearState/NuclearUtils.h"
#include "Physics/Resonance/XSection/DCCSPPPXSec.h"

#include <fstream>
#include <limits>

using namespace genie;
using namespace genie::constants;
using namespace genie::units;
using namespace std::complex_literals;

//____________________________________________________________________________
DCCSPPPXSec::DCCSPPPXSec() :
XSecAlgorithmI("genie::DCCSPPPXSec")
{

}
//____________________________________________________________________________
DCCSPPPXSec::DCCSPPPXSec(std::string config) :
XSecAlgorithmI("genie::DCCSPPPXSec", config)
{

}
//____________________________________________________________________________
DCCSPPPXSec::~DCCSPPPXSec()
{

}
//____________________________________________________________________________
double DCCSPPPXSec::XSec(const Interaction * interaction, KinePhaseSpace_t kps) const
{

    // dimension of kine phase space
    std::string s = KinePhaseSpace::AsString(kps);
    int kpsdim = s!="<|E>"?1 + std::count(s.begin(), s.begin()+s.find('}'), ','):0;
    if (kpsdim < 2 || kpsdim > 4) return 0.;
    // TODO: check 3d-case by integration 4d-case

    if(! this -> ValidProcess    (interaction) ) return 0.;
    if(! this -> ValidKinematics (interaction) ) return 0.;

    double xsec = 0;

    const InitialState & init_state = interaction -> InitState();
    const ProcessInfo &  proc_info  = interaction -> ProcInfo();
    const Target & target = init_state.Tgt();

    // Get kinematical parameters
    const Kinematics & kinematics = interaction -> Kine();
    double Q2    = kinematics.Q2();
    double W     = kinematics.W();
    double Wsq   = W*W;

    //-- Get 1pi exclusive channel
    SppChannel_t spp_channel = SppChannel::FromInteraction(interaction);

    int  probepdgc = init_state.ProbePdg();
    int  nucpdgc   = target.HitNucPdg();
    int  helicity  = init_state.ProbeHelicity();
    bool is_nubar  = pdg::IsAntiNeutrino(probepdgc);
    bool is_EM     = proc_info.IsEM();
    bool is_CC     = proc_info.IsWeakCC();
    bool is_NC     = proc_info.IsWeakNC();
    bool is_EM0    = is_EM && (helicity != 0);
    int  inubnu    = is_nubar?-1:1;
    bool is_p      = pdg::IsProton  (nucpdgc);


    PDGLibrary * pdglib = PDGLibrary::Instance();
    double flepi = 0, flepf = 0;
    if ( !is_EM0 )
    {
       // mass of initial lepton
       flepi    = pdglib->Find(probepdgc)->Mass();
       // mass of final lepton
       flepf    = interaction->FSPrimLepton()->Mass();
    }
    double flepi2   = flepi*flepi;
    double flepf2   = flepf*flepf;
    
    // mass of isoscalar pion
    double fpio     = (pdglib->Find(kPdgPiP)->Mass() + pdglib->Find(kPdgPi0)->Mass() + pdglib->Find(kPdgPiM)->Mass())/3;
    // mass of isoscalar nucleon
    double fnuc     = (pdglib->Find(kPdgProton)->Mass() + pdglib->Find(kPdgNeutron)->Mass())/2;
    double fnuc2    = fnuc*fnuc;

    // energy transfer in the frame of center of mass \piN
    double omegc    = (Wsq - fnuc2 - Q2)/2/W;
    // vec{qc} - momentum transfer in the frame of center of mass \piN
    // vec{qc}^2 - square of momentum transfer in the frame of center of mass \piN
    double qc2sp    = Q2 + omegc*omegc;
    // |vec{qc}|
    double qgamc    = TMath::Sqrt(qc2sp);
    // vec{q} - momentum transfer in the LAB-frame
    // |vec{q}|
    double qgam     = qgamc*W/fnuc;
    // vec{q}^2 - square of momentum transfer in the LAB-frame
    double q2sp     = qgam*qgam;
    // energy transfer in the LAB-frame
    double omeg     = TMath::Sqrt(q2sp - Q2);

    // energy of pion in the frame of center of mass \piN
    double epioc    = (Wsq - fnuc2 + fpio*fpio)/2/W;
    //  vec{k} - momentum of pion in the frame of center of mass \piN
    // |vec{k}|
    double qpioc    = TMath::Sqrt(TMath::Max(epioc*epioc - fpio*fpio, 0.));


    // energy of initial lepton
    double elepi    = init_state.ProbeE(kRfHitNucRest);
    // magnitude of momentum of initial lepton
    double plepi    = TMath::Sqrt(elepi*elepi - flepi2);
    // energy of final lepton
    double elepf    = elepi - omeg;
    // magnitude of momentum of final lepton
    double plepf    = TMath::Sqrt(elepf*elepf - flepf2);
    // theta - scattering angle of final lepton
    // cos(theta)
    double costhl   = (2*elepi*elepf - Q2 - flepi2 - flepf2)/2/plepi/plepf;
    costhl = costhl >  1?1:costhl;
    costhl = costhl < -1?-1:costhl;
    //  cosine of polar angle of pion in the frame of center of mass \piN
    double costhpi = 0;
    if (kpsdim > 2)
        costhpi = kinematics.GetKV(kKVctp);
    // azimuthal angle of pion in the frame of center of mass \piN
    double phipi = 0;
    if (kpsdim == 4)
        phipi = kinematics.GetKV(kKVphip);

    double fcrsc = 0;
    if (is_EM)
        // a factor for EM scattering
        fcrsc = 16*kPi2*kAem2/Q2/Q2/2;
    else if (is_CC)
        // a factor for CC scattering
        fcrsc = kGF2*fVud*fVud*kMw2*kMw2/(kMw2 + Q2)/(kMw2 + Q2);
    else if (is_NC)
        // a factor for NC scattering
        fcrsc = kGF2*kMz2*kMz2/(kMz2 + Q2)/(kMz2 + Q2);

    double fcrsa = plepf/plepi/4/kPi2*fcrsc;


    if (kpsdim == 2)
    {

        double ss2, cc2;
        // speed of initial lepton
        double betai = plepi/elepi;
        // speed of final lepton
        double betaf = plepf/elepf;
        // Further cos(chi) = pf*cos(theta)/Ef,
        // cos(chi)->cos(theta), when mf_0->0
        // ss2=sin^2(chi/2)-mi*mf/Ei/Ef
        // cc2=cos^2(chi/2)+mi*mf/Ei/Ef

        // for EM scattering
        if(is_EM)
        {
            ss2 = (1 - betai*betaf*costhl)/2 - flepf2/elepi/elepf;
            cc2 = (1 + betai*betaf*costhl)/2 + flepf2/elepi/elepf/2;
        }

        // for CC scattering
        if(is_CC)
        {
            ss2 = (1 - betaf*costhl)/2;
            cc2 = (1 + betaf*costhl)/2;
        }

        // for NC scattering
        if(is_NC)
        {
            ss2 = (1 - costhl)/2;
            cc2 = (1 + costhl)/2;
        }
        std::complex<double> zvep, zvem, zvmp, zvmm, zvsp, zvsm, zaep, zaem, zamp, zamm, zalp, zalm, zasp, zasm, zrhp, zrhm, zaxp, zaxm;
        double fact   = 4*W*qpioc/fnuc;
        double facl   = Q2/qgamc/qgamc;
        double rt = 0, rl = 0, rtp = 0, rrh = 0, rrh0 = 0, rrh0i = 0;
        for (int il = 0; il < maxl; il++)
        {
            double l = static_cast<double>(il);
            // vector current
            zvep  = MultipoleV(interaction, 0, l-1);  // E^+_{l-1}
            zvem  = MultipoleV(interaction, 1, l+1);  // E^-_{l+1}
            zvmp  = MultipoleV(interaction, 2, l  );  // M^+_{l}
            zvmm  = MultipoleV(interaction, 3, l  );  // M^-_{l}
            zvsp  = MultipoleV(interaction, 6, l-1);  // S^+_{l-1}
            zvsm  = MultipoleV(interaction, 7, l+1);  // S^-_{l+1}

            // axial vector current
            zaep  = MultipoleA(interaction, 0, l  );  // E^+_{l}
            zaem  = MultipoleA(interaction, 1, l  );  // E^-_{l}
            zamp  = MultipoleA(interaction, 2, l-1);  // M^+_{l-1}
            zamm  = MultipoleA(interaction, 3, l+1);  // M^-_{l+1}
            zalp  = MultipoleA(interaction, 4, l  );  // L^+_{l}
            zalm  = MultipoleA(interaction, 5, l  );  // L^-_{l}
            zasp  = MultipoleA(interaction, 6, l  );  // S^+_{l}
            zasm  = MultipoleA(interaction, 7, l  );  // S^-_{l}

            zrhp  = omegc*zasp - qgamc*zalp;
            zrhm  = omegc*zasm - qgamc*zalm;
            zaxp  = zasp + zrhp*omegc/Q2;
            zaxm  = zasm + zrhm*omegc/Q2;

            rt    += std::real((l+1)*(l+1)*l*(zvmp*std::conj(zvmp) + zvem*std::conj(zvem) + zamm*std::conj(zamm) + zaep*std::conj(zaep)) +
                            (l+1)*l*l*    (zvmm*std::conj(zvmm) + zvep*std::conj(zvep) + zamp*std::conj(zamp) + zaem*std::conj(zaem)));
            rl    += std::real((l+1)*(l+1)*(l+1)*(zvsm*std::conj(zvsm) + zaxp*std::conj(zaxp)) + l*l*l*(zvsp*std::conj(zvsp) + zaxm*std::conj(zaxm)));
            rtp   -= std::real((l+1)*(l+1)*l*(zvmp*std::conj(zaep) + zvem*std::conj(zamm)) - l*l*(l+1)*(zvmm*std::conj(zaem) + zvep*std::conj(zamp)));
            rrh   += std::real((l+1)*(l+1)*(l+1)*zrhp*std::conj(zrhp) + l*l*l*zrhm*std::conj(zrhm));
            rrh0  += std::real((l+1)*(l+1)*(l+1)*zaxp*std::conj(zrhp) + l*l*l*zaxm*std::conj(zrhm));
            rrh0i += std::imag((l+1)*(l+1)*(l+1)*zaxp*std::conj(zrhp) + l*l*l*zaxm*std::conj(zrhm));
        }

        double W1  =  rt*fact/2;
        double W2  =  Q2/q2sp*(rt/2 + rl*facl)*fact;
        double W3  = -2*fnuc/qgam*rtp*fact;
        double W4  =  fnuc2/Q2/Q2*rrh*fact;
        double W5  = -W/q2sp*rrh0*fact;

        xsec = 2*ss2*W1 + cc2*W2;

        if(is_CC || is_NC)
        {
            double pq  = -fnuc*omeg/Q2;
            double W4x = W4 - fnuc2/Q2*W1 + pq*pq*W2 - 2*pq*W5;
            double W5x = W5 - pq*W2;
            xsec += inubnu*W3*((elepi+elepf)*ss2 - flepf2/2/elepf)/fnuc + flepf2/fnuc2*ss2*W4x - flepf2/fnuc/elepf*W5x;
        }

        // a factor for d\sigma3/dE_f dO_f
        double factor = fcrsa*elepi*elepf*2;
        xsec *= factor;

    }
    else
    {
        double conv    = omegc/qgamc;
        std::complex<double> zv1 = 0, zv2 = 0, zv3 = 0, zv4 = 0;
        std::complex<double> zv5 = 0, zv6 = 0, zv7 = 0, zv8 = 0;
        std::complex<double> za1 = 0, za2 = 0, za3 = 0, za4 = 0;
        std::complex<double> za5 = 0, za6 = 0, za7 = 0, za8 = 0;
        std::complex<double> zelp, zelm, zmlp, zmlm, zllp, zllm, zslp, zslm;

        double leg[3][maxl+1];
        CalculateLegendre(costhpi, leg);
        double * dleg  = &leg[1][1];
        double * ddleg = &leg[2][1];


        for (int il = 0; il < maxl; il++)
        {
            double l = static_cast<double>(il);
            double faz11 = -(l+1)*dleg[il];
            double faz12 =  l*dleg[il];
            double faz21 =  (l+1)*dleg[il+1];
            double faz22 = -l*dleg[il-1];
            zelp     = MultipoleV(interaction, 0, il);
            zelm     = MultipoleV(interaction, 1, il);
            zmlp     = MultipoleV(interaction, 2, il);
            zmlm     = MultipoleV(interaction, 3, il);
            zslp     = MultipoleV(interaction, 6, il);
            zslm     = MultipoleV(interaction, 7, il);
            zllp     = zslp*conv;                       // cvc
            zllm     = zslm*conv;                       // cvc
            zv1     += dleg[il+1] *(zelp + l*zmlp) + dleg[il-1]*(zelm + (l+1)*zmlm);
            zv2     += dleg[il]   *((l+1)*zmlp+ l*zmlm);
            zv3     += ddleg[il+1]*(zelp - zmlp) + ddleg[il-1]*(zelm + zmlm);
            zv4     += ddleg[il]  *(zmlp - zelp - zelm - zmlm);
            zv5     += faz21*zllp + faz22*zllm;
            zv6     += faz11*zllp + faz12*zllm;
            zv7     += faz11*zslp + faz12*zslm;
            zv8     += faz21*zslp + faz22*zslm;
        }

        double del   = 0;
        if(is_CC || is_NC)
        {
            del         = flepf*flepf;
            for (int il = 0; il < maxl; il++)
            {
                double l = static_cast<double>(il);
                double faz11 = -(l+1)*dleg[il];
                double faz12 =  l*dleg[il];
                double faz21 =  (l+1)*dleg[il+1];
                double faz22 = -l*dleg[il-1];
                zelp     = MultipoleA(interaction, 0, il);
                zelm     = MultipoleA(interaction, 1, il);
                zmlp     = MultipoleA(interaction, 2, il);
                zmlm     = MultipoleA(interaction, 3, il);
                zllp     = MultipoleA(interaction, 4, il);
                zllm     = MultipoleA(interaction, 5, il);
                zslp     = MultipoleA(interaction, 6, il);
                zslm     = MultipoleA(interaction, 7, il);
                za1     += dleg[il]    *(zelp + zelm + (l+2)*zmlp + (l-1)*zmlm);
                za2     += dleg[il+1]  *(l+1)*zmlp + dleg[il-1]*l*zmlm;
                za3     += ddleg[il]   *(zelp + zelm + zmlp - zmlm);
                za4     += ddleg[il-1] *(zmlm - zelm) - ddleg[il+1]*(zelp + zmlp);
                za5     += faz11*zllp + faz12*zllm;
                za6     += faz21*zllp + faz22*zllm;
                za7     += faz21*zslp + faz22*zslm;
                za8     += faz11*zslp + faz12*zslm;
            }
        }

        double sinthpi   = TMath::Sqrt(1 - costhpi*costhpi);
        std::complex<double> zff[4][4];
        zff[0][1] = -sinthpi*za8;
        zff[0][3] = -1i*(costhpi*za8 + za7);
        zff[3][1] = -sinthpi*za5;
        zff[3][3] = -1i*(costhpi*za5 + za6);
        zff[1][1] = costhpi*za1 - za2 - sinthpi*sinthpi*za3;
        zff[1][3] = -1i*sinthpi*(za1 + za4 + costhpi*za3);
        zff[2][0] = za2 - costhpi*za1;
        zff[2][2] = sinthpi*za1;
        zff[0][0] = 1i*sinthpi*zv7;
        zff[0][2] = 1i*(costhpi*zv7 + zv8);
        zff[3][0] = 1i*sinthpi*zv6;
        zff[3][2] = 1i*(costhpi*zv6 + zv5);
        zff[1][0] = 1i*(zv1 - costhpi*zv2 + sinthpi*sinthpi*zv4);
        zff[1][2] = 1i*sinthpi*(zv2 + zv3 + costhpi*zv4);
        zff[2][1] = 1i*(zv1 - costhpi*zv2);
        zff[2][3] = -sinthpi*zv2;

        std::complex<double> zsum;
        std::complex<double> ztr[4][4];
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
            {
                zsum = 0;
                for (int k = 0; k < 4; k++)
                    zsum += zff[i][k]*std::conj(zff[j][k]);
                ztr[i][j] = zsum;
            }

        double sinthl  = TMath::Sqrt(1 - costhl*costhl);
        double ch      = (fnuc + omeg)/W;
        double sh      = qgam/W;
        double q       = TMath::Sqrt(plepi*plepi + plepf*plepf - 2*plepf*plepi*costhl);
        double plkv0   = (elepi + elepf)/2;
        double plkv3   = (plepi*plepi - plepf*plepf)/2/q;
        
        double plkcv0  = ch*plkv0 - sh*plkv3;
        double plkcv1  = plepi*plepf/q*sinthl;  // plkcv1 = plkv1
        // plkcv2 = plkv2 = 0;
        double plkcv3  = ch*plkv3 - sh*plkv0;
        double xkbq    = plkcv0*qgamc - plkcv3*omegc;
        double xkzc    = plkcv3;
        double xk0c    = plkcv0;
        double xkxc    = plkcv1;
        
        double xbj0j0 = 0, xqjqj = 0, xkjkj = 0;
        std::complex<double> zqj, zbqj, zbj0, zkj;
        std::complex<double> zkjjx = 0, zkjjy  = 0, zbqjjx = 0, zbj0jx = 0, zbqjjy = 0;
        for (int i = 0; i < 4; i++)
        {
            zqj      = omegc*zff[0][i] - qgamc*zff[3][i];
            zbqj     = qgamc*zff[0][i] - omegc*zff[3][i];
            zkj      = xk0c* zff[0][i] - xkzc *zff[3][i];
            zbj0     = zff[0][i] + omegc*zqj/Q2; 
            
            xbj0j0  += std::real(zbj0*std::conj(zbj0));
            xqjqj   += std::real(zqj*std::conj(zqj));
            xkjkj   += std::real(zkj*std::conj(zkj));
            if (kpsdim == 4)
            {
                zbqjjx  += zbqj*std::conj(zff[1][i]);
                zbqjjy  += zbqj*std::conj(zff[2][i]);
                zbj0jx  += zbj0*std::conj(zff[1][i]);
                zkjjx   += zkj *std::conj(zff[1][i]);
                zkjjy   += zkj *std::conj(zff[2][i]);
            }

        }

        double r0     = 0;                                   // RT + RL
        double rc1    = 0;                                   // RLT
        double rs1    = 0;                                   // RLT'
        double rc2    = 0;                                   // RTT
        double rs2    = 0;                                   // RTT'
        
        double eps    = 0;
        if (costhl > -1)
        {
            double tanthl2  = (1 - costhl)/(1 + costhl);
                       eps  = 1/(1 + 2*qgam*qgam*tanthl2/Q2);
        }
        
        double rr1    = std::real(ztr[1][1] + ztr[2][2])/2;
        if (is_EM0)
        {
            double rr3    =  Q2/qgamc/qgamc*xbj0j0;
                   r0     =  rr1 + eps*rr3; 
        }
        else
        {
            double vt1     =  2*xkxc*xkxc + Q2 + del;
            double vt2     = -2*xkbq;
            double rr2     = std::imag(ztr[1][2]);
            double rr11    = xkjkj;
            double rr12    = xqjqj;
            double rr13    = std::real(ztr[0][0] - ztr[3][3]);
            double rta     = vt1*rr1;
            double rtb     = vt2*rr2;
            double rl      = 2*rr11 - (rr12 + (Q2 + del)*rr13)/2;
                   r0      = rta + rl + inubnu*rtb;
        }

        // a factor for d\sigma5/dE_f dO_f dcos_theta_pi
        double factor;
        if (is_EM0)
        {
            double qgamL  = (Wsq - fnuc2)/2/fnuc;
            double Gamma  = kAem*qgamL*elepf/2/kPi2/Q2/elepi/(1 - eps);
            factor = Gamma*qpioc*W/fnuc/qgamL;
        }
        else
        {
            factor = 2*fcrsa*W*qpioc/fnuc;
        }
        

        if (kpsdim == 4)
        {
            double rr8    = std::real(ztr[1][1] - ztr[2][2])/2;
            if (is_EM0)
            {
                double facl = TMath::Sqrt(Q2)/qgamc;
                double rr4  = facl*std::real(zbj0jx);
                double rr7  = facl*std::imag(zbj0jx);
                       rc1  = -TMath::Sqrt(2*eps*(1 + eps))*rr4;
                       rc2  =  eps*rr8;
                       rs1  = -helicity*TMath::Sqrt(2*eps*(1 - eps))*rr7;
            }
            else
            {
                double vrc1a  = -4*xkxc;
                double vrc1b  =  2*xkxc;
                double vrs1a  = -vrc1a;
                double vrs1b  =  vrc1b;
                double vrc2   =  2*xkxc*xkxc;
                double vrs2   =  vrc2;
                double rr9    =-std::real(ztr[1][2]);
                double rr14   = std::real(zkjjx);
                double rr15   = std::real(zkjjy);
                double rr16   = std::imag(zbqjjx);
                double rr17   = std::imag(zbqjjy);
                double rc1a   = vrc1a*rr14;
                double rc1b   = vrc1b*rr17;
                double rs1a   = vrs1a*rr15;
                double rs1b   = vrs1b*rr16;
                       rc2    = vrc2*rr8;
                       rs2    = vrs2*rr9;
                       rc1    = rc1a + inubnu*rc1b;
                       rs1    = rs1a + inubnu*rs1b;
            }
            // a factor for d\sigma5/dE_f dO_f dO_pi
            factor /= 2*kPi;
        }
        
        xsec = r0 + rc1*TMath::Cos(phipi) + rs1*TMath::Sin(phipi) + rc2*TMath::Cos(2*phipi) + rs2*TMath::Sin(2*phipi);
        xsec *= factor;
    }
    // a factor to convert d/dE_fdO_f -> d/dWdQ2
    double factor = kPi*W/fnuc/plepf/plepi;
    xsec *= factor;


  // The algorithm computes d^4xsec/dWdQ2dCosTheta_pidPhi_pi or d^3xsec/dWdQ2dCosTheta_pi or d^2xsec/dWdQ2
  // Check whether variable tranformation is needed
  if ( kps != kPSWQ2ctpphipfE && kps != kPSWQ2ctpfE && kps != kPSWQ2fE)
  {
     double J = 1.;
     if (kpsdim == 2)
       J = utils::kinematics::Jacobian(interaction, kPSWQ2fE, kps);
     else if (kpsdim == 3)
       J = utils::kinematics::Jacobian(interaction, kPSWQ2ctpfE, kps);
     else
       J = utils::kinematics::Jacobian(interaction, kPSWQ2ctpphipfE, kps);
     xsec *= J;
  }

  // Apply given scaling factor
  if (is_EM)
    xsec *= fXSecScaleEM;
  if (is_CC)
    xsec *= fXSecScaleCC;
  if (is_NC)
    xsec *= fXSecScaleNC;

  // If requested return the free nucleon xsec even for input nuclear tgt
  if ( interaction->TestBit(kIAssumeFreeNucleon) ) return xsec;

  int Z = target.Z();
  int A = target.A();
  int N = A-Z;

  // Take into account the number of scattering centers in the target
  int NNucl = (SppChannel::InitStateNucleon(spp_channel) == kPdgProton) ? Z : N;
  xsec*=NNucl; // nuclear xsec (no nuclear suppression symmetry_factor)

  if ( fUsePauliBlocking && A!=1 && kps == kPSWQ2ctpfE )
  {
    // Calculation of Pauli blocking according to refs. 9-11
    double P_Fermi = 0.0;

    // Maximum value of Fermi momentum of target nucleon (GeV)
    if ( A<6 || ! fUseRFGParametrization )
    {
        // look up the Fermi momentum for this target
        FermiMomentumTablePool * kftp = FermiMomentumTablePool::Instance();
        const FermiMomentumTable * kft = kftp->GetTable(fKFTable);
        P_Fermi = kft->FindClosestKF(pdg::IonPdgCode(A, Z), nucpdgc);
     }
     else {
        // define the Fermi momentum for this target
        P_Fermi = utils::nuclear::FermiMomentumForIsoscalarNucleonParametrization(target);
        // correct the Fermi momentum for the struck nucleon
        if(is_p) { P_Fermi *= TMath::Power( 2.*Z/A, 1./3); }
        else     { P_Fermi *= TMath::Power( 2.*N/A, 1./3); }
     }

     double FactorPauli_RES = 1.0;

     if (P_Fermi > 0.)
     {
        if (2*P_Fermi < qgamc - qpioc)
           FactorPauli_RES = 1.0;
        if (2*P_Fermi >= qgamc + qpioc)
           FactorPauli_RES = ((3*qgamc*qgamc + qpioc*qpioc)/2/P_Fermi - (5*TMath::Power(qgamc,4) + TMath::Power(qpioc,4) + 10*qgamc*qgamc*qpioc*qpioc)/(40*TMath::Power(P_Fermi,3)))/2/qgamc;
        if (2*P_Fermi >= qgamc - qpioc && 2*P_Fermi <= qgamc + qpioc)
           FactorPauli_RES = ((qpioc + qgamc)*(qpioc + qgamc) - 4*P_Fermi*P_Fermi/5 - TMath::Power(qgamc - qpioc, 3)/(2*P_Fermi) + TMath::Power(qgamc - qpioc, 5)/(40*TMath::Power(P_Fermi, 3)))/4/qpioc/qgamc;
     }

     xsec *= FactorPauli_RES;
  }
  return xsec;

}
//____________________________________________________________________________
void DCCSPPPXSec::CalculateLegendre(double x, double la[3][maxl+1]) const
{
    double * leg   = &la[0][1];
    double * dleg  = &la[1][1];
    double * ddleg = &la[2][1];

    leg[0]   = 1;
    leg[1]   = x;
    dleg[-1] = 0;
    dleg[0]  = 0;
    dleg[1]  = 1;
    ddleg[-1]= 0;
    ddleg[0] = 0;
    ddleg[1] = 0;

    for (int l = 2; l < maxl; l++)
    {
        leg[l]  = (x*(2*l-1)*leg[l-1] - (l-1)*leg[l-2])/l;
        dleg[l] = l*leg[l-1] + x*dleg[l-1];
        ddleg[l]= (l+1)*dleg[l-1]+ x*ddleg[l-1];
    }

}
//____________________________________________________________________________
std::string DCCSPPPXSec::FindDataTableFile(const std::string &basename, bool &ok) const
{

  ok = true;
  for (size_t p = 0; p < fDataPaths.size(); ++p)
  {
    const std::string& path = fDataPaths.at( p );
    std::string full_name = path + '/' + basename;
    // is file exist?
    if ( std::ifstream(full_name.c_str()).good() )
      return full_name;
  }

  // A matching file could not be found
  ok = false;
  return std::string();

}
//____________________________________________________________________________
void DCCSPPPXSec::ReadMultipoleTable(void)
{

    bool table_ok;
    std::string full_file_name = FindDataTableFile(fDataFileName, table_ok);

    if ( table_ok )
    {
        LOG("DCCSPPPXSec", pINFO) <<
            "Loading the table with ANL-Osaka multipole amplitudes from file " << full_file_name;

        std::ifstream file(full_file_name);
        char cdum[20];
        int mxw, mxq, mxl;
        file >> mxw >> mxq >> mxl >> cdum;
        for (int iW = 0; iW < maxW; iW++)
            for (int iQ2 = 0; iQ2 < maxQ2; iQ2++)
            {
                file >> Wnodes[iW] >> Q2nodes[iQ2] >> cdum >> cdum;
                for (int imu = 0; imu < maxmu; imu++)
                    for (int iso = 0; iso < 3; iso++)
                    {
                        file.ignore(11);
                        for (int il = 0; il < maxl; il++)
                            file >> MultipoleTbl[MultipoleTblIndx(imu, il, iso, 0, iQ2, 0, iW)] >> MultipoleTbl[MultipoleTblIndx(imu, il, iso, 1, iQ2, 0, iW)];
                    }
                for (int imu = 0; imu < maxmu; imu++)
                    for (int iso = 3; iso < maxiso; iso++)
                    {
                        file.ignore(11);
                        for (int il = 0; il < maxl; il++)
                            file >> MultipoleTbl[MultipoleTblIndx(imu, il, iso, 0, iQ2, 0, iW)] >> MultipoleTbl[MultipoleTblIndx(imu, il, iso, 1, iQ2, 0, iW)];
                    }
            }
        file.close();
    }
    else
    {
        LOG("DCCSPPPXSec", pERROR) <<
            "Couldn't load the table with ANL-Osaka multipole amplitudes from file " << full_file_name;
        std::exit(EXIT_FAILURE);
    }

}
//____________________________________________________________________________
void DCCSPPPXSec::Spline (const int n, const double * x, const double * y, double * b, double * c, double * d) const
// Calculate spline coefficients bi, ci, di for cubic spline Si(x) = yi + bi(x - xi) + ci(x-xi)^2 + di(x - xi)^3 for a function given in points xi, yi
/*!
    @param[in]  n    - number of points
    @param[in]  x    - pointer to array with \f$x_i\f$
    @param[in]  y    - pointer to array with \f$y_i\f$
    @param[out] b    - pointer to array with \f$b_i\f$
    @param[out] c    - pointer to array with \f$c_i\f$
    @param[out] d    - pointer to array with \f$d_i\f$
*/
{
    int nm1 = n-1;
    //  set up tridiagonal system
    //  b = diagonal, d = offdiagonal, c = right hand side.
    d[0] = x[1] - x[0];
    c[1] = (y[1] - y[0])/d[0];
    for (int i = 1; i < nm1; i++)
    {
        d[i]   = x[i+1] - x[i];
        b[i]   = 2*(d[i-1] + d[i]);
        c[i+1] = (y[i+1] - y[i])/d[i];
        c[i]   = c[i+1] - c[i];
    }

    // end conditions, third derivatives at x[0] and x[n]
    // obtained from divided differences
      b[0]   = -d[0];
      b[n-1] = -d[n-2];
      c[0]   = c[2]/(x[3]-x[1]) - c[1]/(x[2]-x[0]);
      c[n-1] = c[n-2]/(x[n-1]-x[n-3]) - c[n-3]/(x[n-2]-x[n-4]);
      c[0]   = c[0]*d[0]*d[0]/(x[3]-x[0]);
      c[n-1] = -c[n-1]*d[n-2]*d[n-2]/(x[n-1]-x[n-4]);

    // forward elimination
    for (int i = 1; i < n; i++)
    {
        double t = d[i-1]/b[i-1];
        b[i] -= t*d[i-1];
        c[i] -= t*c[i-1];
    }

    // back substitution
    c[n-1] = c[n-1]/b[n-1];
    for (int i = nm1-1; i >= 0; i--)
        c[i] = (c[i] - d[i]*c[i+1])/b[i];

    // compute polynomial coefficients
    b[n-1] = (y[n-1] - y[nm1-1])/d[nm1-1] + d[nm1-1]*(c[nm1-1] + 2*c[n-1]);
    for (int i = 0; i < nm1; i++)
    {
        b[i] = (y[i+1] - y[i])/d[i] - d[i]*(c[i+1] + 2*c[i]);
        d[i] = (c[i+1] - c[i])/d[i];
        c[i] *= 3;
    }
    c[n-1] *= 3;
    d[n-1] = d[n-2];
}
//____________________________________________________________________________
void DCCSPPPXSec::InitializeSplineCoefficients()
{
    for (int imu = 0; imu < maxmu; imu++)
        for (int il = 0; il < maxl; il++)
            for (int iso = 0; iso < maxiso; iso++)
                for (int izp = 0; izp < maxzp; izp++)
                    for (int iQ2 = 0; iQ2 < maxQ2; iQ2++)
                        Spline (maxW, &Wnodes[0],
                                        &MultipoleTbl[MultipoleTblIndx(imu, il, iso, izp, iQ2, 0, 0)],
                                        &MultipoleTbl[MultipoleTblIndx(imu, il, iso, izp, iQ2, 1, 0)],
                                        &MultipoleTbl[MultipoleTblIndx(imu, il, iso, izp, iQ2, 2, 0)],
                                        &MultipoleTbl[MultipoleTblIndx(imu, il, iso, izp, iQ2, 3, 0)]);
}
//____________________________________________________________________________
int DCCSPPPXSec::MultipoleTblIndx (int imu, int il, int iso, int izpart, int iQ2, int icf, int iW) const
// Calculate position in table with ANL-Osaka multipole amplitudes
/*!
    @param[in] imu    - multipole index: 0 -\f$E_{l+}\f$, 1 -\f$E_{l-}\f$, 2 -\f$M_{l+}\f$, 3 -\f$M_{l-}\f$, 4 -\f$S_{l+}\f$, 5 -\f$S_{l-}\f$, 6 -\f$L_{l+}\f$, 7 -\f$L_{l-}\f$
    @param[in] il     - \f$l\f$ (orbital momentum from 0 to DCCSPPPXSec::maxl - 1)
    @param[in] iso    - isospin index: 0 -\f$a^V_{\frac{3}{2}}\f$, 1 -\f$a^V_{\frac{1}{2}}\f$, 2 -\f$a^V_{0}\f$, 3 -\f$a^A_{\frac{3}{2}}\f$, 4 -\f$a^A_{\frac{1}{2}}\f$
    @param[in] izpart - 0 - Re, 1 - Im
    @param[in] iQ2    - index of \f$Q^2\f$-node
    @param[in] icf    - cubic spline (\f$S_\textrm{iW}(W) = y_\textrm{iW}+b_\textrm{iW}(W-W_\textrm{iW})+{c_\textrm{iW}}(W-W_\textrm{iW})^2+d_\textrm{iW}(W-W_\textrm{iW})^3\f$) coefficient for <b>iW</b>-node: 0 -\f$y_\textrm{iW}\f$, 1 -\f$b_\textrm{iW}\f$, 2 -\f$c_\textrm{iW}\f$, 3 -\f$d_\textrm{iW}\f$
    @param[in] iW     - index of \f$W\f$-node
*/
{
    return (((((imu*maxl + il)*maxiso + iso)*maxzp + izpart)*maxQ2 + iQ2)*maxcf + icf)*maxW + iW;
}
//____________________________________________________________________________
int DCCSPPPXSec::Wnode (double W) const
{
    const int N = 11;
    double pos[N][2] = { {1.07701, 0}, {1.08, 1}, {1.16, 9}, {1.21, 19}, {1.22, 23}, {1.225, 24},
                         {1.2255, 25}, {1.23, 26}, {1.24, 28}, {2.0, 104}, {2.1, 109} };
    if (W < pos[0][0])   return pos[0][1];
    if (W > pos[N-1][0]) return pos[N-1][1];
    int low = 0;
    int up = N;
    while (up-low > 1)
    {
        int inx = (low + up)/2;
        if (W<=pos[inx][0])
            up = inx;
        else
            low =inx;
    }
    return TMath::FloorNint(pos[low][1] + (W-pos[low][0])*(pos[up][1]-pos[low][1])/(pos[up][0]-pos[low][0]));
}
//____________________________________________________________________________
int DCCSPPPXSec::Q2node (double Q2) const
{
    const int N = 6;
    double pos[N][2] = { {1E-7, 0}, {0.02, 1}, {0.2, 10}, {1.0, 18}, {2.0, 23}, {3.0, 27}};
    if (Q2 < pos[0][0])   return pos[0][1];
    if (Q2 > pos[N-1][0]) return pos[N-1][1];
    int low = 0;
    int up = N;
    while (up-low > 1)
    {
        int inx = (low + up)/2;
        if (Q2<=pos[inx][0])
            up = inx;
        else
            low =inx;
    }
    return TMath::FloorNint(pos[low][1] + (Q2-pos[low][0])*(pos[up][1]-pos[low][1])/(pos[up][0]-pos[low][0]));
}
//____________________________________________________________________________
double DCCSPPPXSec::IsospinAmplitude(const Interaction * interaction, int iso) const
/*!
    @param[in] iso - isospin index: 0 -\f$a^V_{\frac{3}{2}}\f$, 1 -\f$a^V_{\frac{1}{2}}\f$, 2 -\f$a^V_{0}\f$, 3 -\f$a^A_{\frac{3}{2}}\f$, 4 -\f$a^A_{\frac{1}{2}}\f$
*/
{

     //-- Get 1pi exclusive channel
    SppChannel_t spp_chn = SppChannel::FromInteraction(interaction);

    double val = 0;

    if (iso == 0 || iso == 3)
        switch (spp_chn)
        {
            // EM, NC
            case kSpp_lp_em_10010:  // ich = 1
            case kSpp_vp_nc_10010:  // ich = 1
            case kSpp_vbp_nc_10010: // ich = 1
            case kSpp_ln_em_01010:  // ich = 3
            case kSpp_vn_nc_01010:  // ich = 3
            case kSpp_vbn_nc_01010: // ich = 3
                val = 2./3;
                break;
            case kSpp_lp_em_01100:  // ich = 2
            case kSpp_vp_nc_01100:  // ich = 2
            case kSpp_vbp_nc_01100: // ich = 2
            case kSpp_ln_em_10001:  // ich = 4
            case kSpp_vn_nc_10001:  // ich = 4
            case kSpp_vbn_nc_10001: // ich = 4
                val = TMath::Sqrt(2)/3;
                break;
            // CC
            case kSpp_vp_cc_10100:   // ich = 1
            case kSpp_vbn_cc_01001:  // ich = 1
                val = TMath::Sqrt(2);
                break;
            case kSpp_vn_cc_01100:   // ich = 2
            case kSpp_vbp_cc_10001:  // ich = 2
                val = TMath::Sqrt(2)/3;
                break;
            case kSpp_vn_cc_10010:   // ich = 3
            case kSpp_vbp_cc_01010:  // ich = 3
                val = 2./3;
                break;
            default:
                val = 0;
        }
    else if (iso == 1 || iso == 4)
        switch (spp_chn)
        {
            // EM, NC
            case kSpp_lp_em_10010:  // ich = 1
            case kSpp_vp_nc_10010:  // ich = 1
            case kSpp_vbp_nc_10010: // ich = 1
            case kSpp_ln_em_01010:  // ich = 3
            case kSpp_vn_nc_01010:  // ich = 3
            case kSpp_vbn_nc_01010: // ich = 3
                val = 1./3;
                break;
            case kSpp_lp_em_01100:  // ich = 2
            case kSpp_vp_nc_01100:  // ich = 2
            case kSpp_vbp_nc_01100: // ich = 2
            case kSpp_ln_em_10001:  // ich = 4
            case kSpp_vn_nc_10001:  // ich = 4
            case kSpp_vbn_nc_10001: // ich = 4
                val = -TMath::Sqrt(2)/3;
                break;
            // CC
            case kSpp_vn_cc_01100:   // ich = 2
            case kSpp_vbp_cc_10001:  // ich = 2
                val = 2*TMath::Sqrt(2)/3;
                break;
            case kSpp_vn_cc_10010:   // ich = 3
            case kSpp_vbp_cc_01010:  // ich = 3
                val = -2./3;
                break;
            default:
                val = 0;
        }
    else if (iso == 2)
        switch (spp_chn)
        {
            // EM, NC
            case kSpp_lp_em_10010:  // ich = 1
            case kSpp_vp_nc_10010:  // ich = 1
            case kSpp_vbp_nc_10010: // ich = 1
                val = 1;
                break;
            case kSpp_lp_em_01100:  // ich = 2
            case kSpp_vp_nc_01100:  // ich = 2
            case kSpp_vbp_nc_01100: // ich = 2
                val = -TMath::Sqrt(2);
                break;
            case kSpp_ln_em_01010:  // ich = 3
            case kSpp_vn_nc_01010:  // ich = 3
            case kSpp_vbn_nc_01010: // ich = 3
                val = -1;
                break;
            case kSpp_ln_em_10001:  // ich = 4
            case kSpp_vn_nc_10001:  // ich = 4
            case kSpp_vbn_nc_10001: // ich = 4
                val = TMath::Sqrt(2);
                break;
            default:
                val = 0;
        }

    double sgn = -1;
    const InitialState & init_state = interaction -> InitState();
    if (pdg::IsAntiNeutrino (init_state.ProbePdg()))
        sgn = 1;

    const ProcessInfo &  proc_info  = interaction -> ProcInfo();
    bool is_EM     = proc_info.IsEM();
    bool is_CC     = proc_info.IsWeakCC();
    bool is_NC     = proc_info.IsWeakNC();

    if (is_CC)
        val *= sgn; //account for neutrino/antineutrino case

    if (is_EM && iso > 2)
        val = 0;   // no axial-vector part for EM-processes

    if (is_NC)
    {
        double fx = 1;
        if (iso < 2)
            fx = (1 - 2*fSin2Wein);
        else if (iso == 2)
            fx =  -2*fSin2Wein;
        val *= fx;
    }

    return val;

}
//____________________________________________________________________________
std::complex<double> DCCSPPPXSec::MultipoleV(const Interaction * interaction, int mult, int l) const
/*!
    @param[in] mult  - multipole index: 0 -\f$E_{l+}\f$, 1 -\f$E_{l-}\f$, 2 -\f$M_{l+}\f$, 3 -\f$M_{l-}\f$, 4 -\f$S_{l+}\f$, 5 -\f$S_{l-}\f$, 6 -\f$L_{l+}\f$, 7 -\f$L_{l-}\f$
    @param[in] l     - \f$l\f$ (orbital momentum from 0 to DCCSPPPXSec::maxl)
*/
{
    if (l < 0 || l >= maxl)
        return std::complex<double>(0, 0);

    // Get kinematical parameters
    const Kinematics & kinematics = interaction -> Kine();
    double Q2    = kinematics.Q2();
    double W     = kinematics.W();

    int iW = Wnode(W);

    double re = 0, im = 0;
    for (int iso = 0; iso < 3;iso++)
    {
        double arQ2[4][maxQ2];
        for (int iQ2 = 0; iQ2 < maxQ2; iQ2++)
        {
            double a = MultipoleTbl[MultipoleTblIndx(mult, l, iso, 0, iQ2, 0, iW)];
            double b = MultipoleTbl[MultipoleTblIndx(mult, l, iso, 0, iQ2, 1, iW)];
            double c = MultipoleTbl[MultipoleTblIndx(mult, l, iso, 0, iQ2, 2, iW)];
            double d = MultipoleTbl[MultipoleTblIndx(mult, l, iso, 0, iQ2, 3, iW)];
            arQ2[0][iQ2] = a + ( b + ( c + d*(W - Wnodes[iW]) )*(W - Wnodes[iW]) )*(W - Wnodes[iW]);
        }
        Spline (maxQ2, &Q2nodes[0], &arQ2[0][0], &arQ2[1][0], &arQ2[2][0], &arQ2[3][0]);
        int iQ2 = Q2node(Q2);
        double val = arQ2[0][iQ2] + ( arQ2[1][iQ2] + ( arQ2[2][iQ2] + arQ2[3][iQ2]*(Q2 - Q2nodes[iQ2]) )*(Q2 - Q2nodes[iQ2]) )*(Q2 - Q2nodes[iQ2]);
        re += val*IsospinAmplitude(interaction, iso);
    }
    for (int iso = 0; iso < 3;iso++)
    {
        double arQ2[4][maxQ2];
        for (int iQ2 = 0; iQ2 < maxQ2; iQ2++)
        {
            double a = MultipoleTbl[MultipoleTblIndx(mult, l, iso, 1, iQ2, 0, iW)];
            double b = MultipoleTbl[MultipoleTblIndx(mult, l, iso, 1, iQ2, 1, iW)];
            double c = MultipoleTbl[MultipoleTblIndx(mult, l, iso, 1, iQ2, 2, iW)];
            double d = MultipoleTbl[MultipoleTblIndx(mult, l, iso, 1, iQ2, 3, iW)];
            arQ2[0][iQ2] = a + ( b + ( c + d*(W - Wnodes[iW]) )*(W - Wnodes[iW]) )*(W - Wnodes[iW]);
        }
        Spline (maxQ2, &Q2nodes[0], &arQ2[0][0], &arQ2[1][0], &arQ2[2][0], &arQ2[3][0]);
        int iQ2 = Q2node(Q2);
        double val = arQ2[0][iQ2] + ( arQ2[1][iQ2] + ( arQ2[2][iQ2] + arQ2[3][iQ2]*(Q2 - Q2nodes[iQ2]) )*(Q2 - Q2nodes[iQ2]) )*(Q2 - Q2nodes[iQ2]);
        im += val*IsospinAmplitude(interaction, iso);
    }

    return std::complex<double>(re, im);
}
//____________________________________________________________________________
std::complex<double> DCCSPPPXSec::MultipoleA(const Interaction * interaction, int mult, int l) const
/*!
    @param[in] mult  - multipole index: 0 -\f$E_{l+}\f$, 1 -\f$E_{l-}\f$, 2 -\f$M_{l+}\f$, 3 -\f$M_{l-}\f$, 4 -\f$S_{l+}\f$, 5 -\f$S_{l-}\f$, 6 -\f$L_{l+}\f$, 7 -\f$L_{l-}\f$
    @param[in] l     - \f$l\f$ (orbital momentum from 0 to DCCSPPPXSec::maxl)
*/
{
    if (l < 0 || l >= maxl)
        return std::complex<double>(0, 0);

    // Get kinematical parameters
    const Kinematics & kinematics = interaction -> Kine();
    double Q2    = kinematics.Q2();
    double W     = kinematics.W();

    int iW = Wnode(W);

    double re = 0, im = 0;
    for (int iso = 3; iso < maxiso;iso++)
    {
        double arQ2[4][maxQ2];
        for (int iQ2 = 0; iQ2 < maxQ2; iQ2++)
        {
            double a = MultipoleTbl[MultipoleTblIndx(mult, l, iso, 0, iQ2, 0, iW)];
            double b = MultipoleTbl[MultipoleTblIndx(mult, l, iso, 0, iQ2, 1, iW)];
            double c = MultipoleTbl[MultipoleTblIndx(mult, l, iso, 0, iQ2, 2, iW)];
            double d = MultipoleTbl[MultipoleTblIndx(mult, l, iso, 0, iQ2, 3, iW)];
            arQ2[0][iQ2] = a + ( b + ( c + d*(W - Wnodes[iW]) )*(W - Wnodes[iW]) )*(W - Wnodes[iW]);
        }
        Spline (maxQ2, &Q2nodes[0], &arQ2[0][0], &arQ2[1][0], &arQ2[2][0], &arQ2[3][0]);
        int iQ2 = Q2node(Q2);
        double val = arQ2[0][iQ2] + ( arQ2[1][iQ2] + ( arQ2[2][iQ2] + arQ2[3][iQ2]*(Q2 - Q2nodes[iQ2]) )*(Q2 - Q2nodes[iQ2]) )*(Q2 - Q2nodes[iQ2]);
        re += val*IsospinAmplitude(interaction, iso);
    }
    for (int iso = 3; iso < maxiso;iso++)
    {
        double arQ2[4][maxQ2];
        for (int iQ2 = 0; iQ2 < maxQ2; iQ2++)
        {
            double a = MultipoleTbl[MultipoleTblIndx(mult, l, iso, 1, iQ2, 0, iW)];
            double b = MultipoleTbl[MultipoleTblIndx(mult, l, iso, 1, iQ2, 1, iW)];
            double c = MultipoleTbl[MultipoleTblIndx(mult, l, iso, 1, iQ2, 2, iW)];
            double d = MultipoleTbl[MultipoleTblIndx(mult, l, iso, 1, iQ2, 3, iW)];
            arQ2[0][iQ2] = a + ( b + ( c + d*(W - Wnodes[iW]) )*(W - Wnodes[iW]) )*(W - Wnodes[iW]);
        }
        Spline (maxQ2, &Q2nodes[0], &arQ2[0][0], &arQ2[1][0], &arQ2[2][0], &arQ2[3][0]);
        int iQ2 = Q2node(Q2);
        double val = arQ2[0][iQ2] + ( arQ2[1][iQ2] + ( arQ2[2][iQ2] + arQ2[3][iQ2]*(Q2 - Q2nodes[iQ2]) )*(Q2 - Q2nodes[iQ2]) )*(Q2 - Q2nodes[iQ2]);
        im += val*IsospinAmplitude(interaction, iso);
    }

    return std::complex<double>(re, im);
}
//____________________________________________________________________________
double DCCSPPPXSec::Integral(const Interaction * interaction) const
{
  double xsec = fXSecIntegrator->Integrate(this,interaction);
  return xsec;
}
//____________________________________________________________________________
bool DCCSPPPXSec::ValidProcess(const Interaction * interaction) const
{

  if(interaction->TestBit(kISkipProcessChk)) return true;

  //-- Get the requested SPP channel
  SppChannel_t spp_channel = SppChannel::FromInteraction(interaction);
  if( spp_channel == kSppNull ) {
    return false;
  }

  return true;

}
//____________________________________________________________________________
bool DCCSPPPXSec::ValidKinematics(const Interaction * interaction) const
{
    // call only after ValidProcess
    if ( interaction->TestBit(kISkipKinematicChk) ) return true;

    const KPhaseSpace& kps = interaction->PhaseSpace();

    // Get kinematical parameters
    const InitialState & init_state = interaction -> InitState();
    const Kinematics & kinematics   = interaction -> Kine();
    const ProcessInfo &  proc_info  = interaction -> ProcInfo();
    
    bool is_EM     = proc_info.IsEM();
    int  helicity  = init_state.ProbeHelicity();
    bool is_EM0    = is_EM && (helicity != 0);
    
    double Enu  = init_state.ProbeE(kRfHitNucRest);
    double W    = kinematics.W();
    double Q2   = kinematics.Q2();

    if (Enu < kps.Threshold_SPP_iso(is_EM0))
        return false;

    Range1D_t Wl  = kps.WLim_SPP_iso(is_EM0);
    Range1D_t Q2l = kps.Q2Lim_W_SPP_iso(is_EM0);

    // model restrictions
    Wl.min  = TMath::Max (Wl.min,  1.077);
    Wl.max  = TMath::Max (Wl.max,  1.077);
    Wl.max  = TMath::Min (Wl.max,  2.1);
    Q2l.min = TMath::Max (Q2l.min, 0.0);
    Q2l.max = TMath::Min (Q2l.max, 3.0);

    if (W < Wl.min || W > Wl.max)
        return false;

    if (Q2 < Q2l.min || Q2 > Q2l.max)
        return false;

    return true;
}
//____________________________________________________________________________
void DCCSPPPXSec::Configure(const Registry & config)
{
    Algorithm::Configure(config);
    this->LoadConfig();
}
//____________________________________________________________________________
void DCCSPPPXSec::Configure(std::string config)
{
    Algorithm::Configure(config);
    this->LoadConfig();
}
//____________________________________________________________________________
void DCCSPPPXSec::LoadConfig(void)
{
    // Cross section scaling symmetry_factors
    this->GetParam( "DCC-EM-XSecScale", fXSecScaleEM );
    this->GetParam( "DCC-CC-XSecScale", fXSecScaleCC );
    this->GetParam( "DCC-NC-XSecScale", fXSecScaleNC );

    double thw;
    this->GetParam( "WeinbergAngle", thw );
    fSin2Wein = TMath::Power( TMath::Sin(thw), 2 );

    this->GetParam("CKM-Vud", fVud );


    // Either a data path relative to the root GENIE folder
    // or an absolute path can be used. Find out which
    // option was chosen.
    std::string path_type;
    GetParamDef( "DataPathType", path_type, std::string("relative") );

    // Right now, there can only be a single data path
    // specified. We use a vector of paths to allow for
    // easy expansion later.
    std::string data_path;
    GetParam( "DataPath", data_path );

    // Convert the relative path to an absolute one if needed
    if ( path_type == "relative" ) {
    data_path = std::string( gSystem->Getenv("GENIE") ) + '/' + data_path;
    }

    fDataPaths.push_back( data_path );

    GetParam( "DataFileName", fDataFileName );

    GetParam("FermiMomentumTable", fKFTable);
    GetParam("RFG-UseParametrization", fUseRFGParametrization);
    GetParam("UsePauliBlockingForRES", fUsePauliBlocking);

    // Load the differential cross section integrator
    fXSecIntegrator = dynamic_cast<const XSecIntegratorI *> (this->SubAlg("XSec-Integrator"));
    assert(fXSecIntegrator);
    
    // Read table with amplitudes and initialize spline coefficients
    ReadMultipoleTable();
    InitializeSplineCoefficients();

}
