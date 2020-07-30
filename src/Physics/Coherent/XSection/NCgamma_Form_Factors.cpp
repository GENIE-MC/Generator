//
// Created by edusaul on 8/04/19.
//

#include "NCgamma_Form_Factors.h"
#include <cmath>
#include "NCgamma_Parameters_GeV.h"
#include <iostream>

using namespace NC_gamma;

const std::complex<double> &Nuclear_FF::getFfP() const {
    return FF_p;
}

const std::complex<double> &Nuclear_FF::getFfN() const {
    return FF_n;
}


/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

Nucleus_FF_DeVries::Nucleus_FF_DeVries(const std::string &n) : nucleus(n) {
    this->pi= NCgamma_Param::pi;
    this->hc=NCgamma_Param::hc;

    if(this->nucleus == "12C") this->set_12C();
    else if(this->nucleus == "40Ar") this->set_40Ar();
    else {
        std::cout<<"WARNING: wrong nucleus, setting 12C as default..."<<std::endl;
        this->set_12C();
    }

}


void Nucleus_FF_DeVries::setFF(double Q) {

    double aux,nu;
    int i;

    aux=0.0;

    for (i=1;i<16; i++){
        nu=double(i);
        aux += pow(-1.0,nu+1.0)*R*R*R*a[i]*hc*hc*hc/(pi*pi*nu*nu-Q*Q*R*R);
    }

    this->FF_p = 4.0*pi*sin(Q*R)/(Q*R)*aux;// /6.0;
    this->FF_n = this->FF_p;
}

void Nucleus_FF_DeVries::set_12C() {
    //a(i) en fm^-3
    this->a[1]=0.15737e-1;
    this->a[2]=0.38897e-1;
    this->a[3]=0.37085e-1;
    this->a[4]=0.14795e-1;
    this->a[5]=-0.44834e-2;

    this->a[6]=-0.10057e-1;
    this->a[7]=-0.68696e-2;
    this->a[8]=-0.28813e-2;
    this->a[9]=-0.77229e-3;
    this->a[10]=0.66908e-4;

    this->a[11]=0.10636e-3;
    this->a[12]=-0.36864e-4;
    this->a[13]=-0.50135e-5;
    this->a[14]=0.9455e-5;
    this->a[15]=-0.47687e-5;

    this->R=8.0/hc; // [fm]/hc
}

void Nucleus_FF_DeVries::set_40Ar() {
    //a(i) en fm^-3
    this->a[1]=0.30451e-1;
    this->a[2]=0.55337e-1;
    this->a[3]=0.20203e-1;
    this->a[4]=-0.16765e-1;
    this->a[5]=-0.13578e-1;

    this->a[6]=-0.43204e-4;
    this->a[7]=0.91988e-3;
    this->a[8]=-0.41205e-3;
    this->a[9]=0.11971e-3;
    this->a[10]=-0.19801e-4;

    this->a[11]=-0.43204e-5;
    this->a[12]=0.61205e-5;
    this->a[13]=-0.37803e-5;
    this->a[14]=0.18001e-5;
    this->a[15]=-0.77407e-6;

    this->R=9.0/hc; // [fm]/hc
}

Nucleus_FF_DeVries::~Nucleus_FF_DeVries() {}



/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////




Form_Factors_Delta::Form_Factors_Delta() {
    this->mn = NCgamma_Param::mp;
    this->mn2 = NCgamma_Param::mp2;
    this->mDelta = NCgamma_Param::mDelta;
    this->mDelta2 = this->mDelta*this->mDelta;
    this->N_Delta_MA = NCgamma_Param::N_Delta_MA;
    this->aNC = 1.0 - 2.0*NCgamma_Param::sinw2;

    this->kgcm0 = (this->mDelta2-this->mn2)/2.0/this->mDelta;
    this->mpw2 = (this->mn+this->mDelta)*(this->mn+this->mDelta);
    this->mmw2 = (this->mn - this->mDelta)*(this->mn - this->mDelta);

    this->C4v = 0.0;
    this->C5v = 0.0;
    this->C3a = 0.0;
    this->C4a = 0.0;
    this->C5a = 0.0;
    this->C6a = 0.0;

    this->C4vNC = 0.0;
    this->C5vNC = 0.0;
    this->C3aNC = 0.0;
    this->C4aNC = 0.0;
    this->C6aNC = 0.0;

}

void Form_Factors_Delta::setFF(double Q2) {

    double egcm = (this->mDelta2-Q2-this->mn2)/2.0/this->mDelta;
    double qcm = sqrt(egcm*egcm + Q2);

    double Fq = 1.0/((1.0 + Q2/NCgamma_Param::parameter_071)*(1.0 + Q2/NCgamma_Param::parameter_071));
    double AM = NCgamma_Param::parameter_03*(1.0 + NCgamma_Param::parameter_001*Q2)*exp(-NCgamma_Param::parameter_023*Q2)*(qcm/this->kgcm0)*Fq;
    double A32= sqrt(3.0)/2.0*(-AM);

    double r = sqrt(2.0*this->kgcm0/NCgamma_Param::pi/NCgamma_Param::alpha*this->mn*this->mDelta/(this->mmw2+Q2));

    this->C3v=-r*this->mn*this->mDelta/(this->mpw2+Q2)*(2.0*A32);

    double Fd=pow((1.0 + Q2/(this->N_Delta_MA*this->N_Delta_MA)),-2);
    this->C5aNC=1.2*Fd;
    this->C3vNC = this->C3v * this->aNC;

}

double Form_Factors_Delta::getC3v() const {
    return C3v;
}

double Form_Factors_Delta::getC4v() const {
    return C4v;
}

double Form_Factors_Delta::getC5v() const {
    return C5v;
}

double Form_Factors_Delta::getC3a() const {
    return C3a;
}

double Form_Factors_Delta::getC4a() const {
    return C4a;
}

double Form_Factors_Delta::getC5a() const {
    return C5a;
}

double Form_Factors_Delta::getC6a() const {
    return C6a;
}

double Form_Factors_Delta::getC3vNC() const {
    return C3vNC;
}

double Form_Factors_Delta::getC4vNC() const {
    return C4vNC;
}

double Form_Factors_Delta::getC5vNC() const {
    return C5vNC;
}

double Form_Factors_Delta::getC3aNC() const {
    return C3aNC;
}

double Form_Factors_Delta::getC4aNC() const {
    return C4aNC;
}

double Form_Factors_Delta::getC5aNC() const {
    return C5aNC;
}

double Form_Factors_Delta::getC6aNC() const {
    return C6aNC;
}

