#ifndef _PARAMETERS_H_
#define _PARAMETERS_H_
#include <cmath>
#include <complex>

namespace NCgamma_Param
{

    const std::complex<double> c_i(0.0,1.0);
    const double pi= 4.0*atan(1.0);
    const double alpha= 1.0/137.0;
    const double hc=0.19733; // (GeV x fm)
    const double hccm=197.33e-16; // (GeV x cm)
    const double Gf=1.16639e-5; // GeV^-2
    const double cosw2=0.23122;    // Weak angle squared
    const double sinw2=0.2223;
    const double ga=1.2723; //from the PDG value for gA/gV (measured in beta decay: http://pdglive.lbl.gov/DataBlock.action?node=S017AV&init=0) using gV=1 as expected.
    const double gas=0.15;
    const double MA = 1.049; // GeV Axial mass

    //nuclear matter density (in fm^-3 * hc^3)
    const double rho0=0.17 * hc*hc*hc;

    //lepton mases
    const double me= 0.51099892e-3; //GeV
    const double mmu= 0.1056583715; //GeV

    //messon mass
    const double mpi=0.13957; //GeV

    //nucleon
    const double mp= 0.93827203; //GeV
    const double mn= 0.93827203; //GeV
    const double mup = 1.7928; //proton anomalous magnetic moment
    const double mun = -1.913; //neutron anomalous magnetic moment
        //for the form factors:
    const double lambdanN=5.6;
    const double MD = 0.843; // GeV

    //barion mass
    const double mLambda_1115= 1.115; //GeV
    const double mSigma_1190= 1.197449; //GeV
    const double mDelta= 1.232; //GeV
    const double mN1440= 1.440; //GeV
    const double mN1520= 1.515; //GeV
    const double mN1535= 1.530; //GeV

    const double coupling_DeltaNpi= 2.14;
    const double N_Delta_MA= 0.93; //GeV
    const double Delta_V0= 0.08; //GeV

    //Carbon nucleus
    const double AC= 12.0;
    const double ZC= 6.0;
    const double MC= 12.0*mp;

    //Powers

    const double hccm2= hccm*hccm;
    const double Gf2=Gf*Gf;
    const double mpi2=mpi*mpi;
    const double mp2=mp*mp;

    //Delta FF parameters
    const double parameter_001 = 0.01; //GeV^-2
    const double parameter_023 = 0.23; //GeV^-2
    const double parameter_071 = 0.71; //GeV^2
    const double parameter_093 = 0.93; //GeV^2
    const double parameter_03 = 0.3; //GeV^(-1/2)


}

#endif
