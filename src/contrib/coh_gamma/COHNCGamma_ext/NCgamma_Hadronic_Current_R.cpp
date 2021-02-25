//
// Created by edusaul on 11/03/19.
//

#include "NCgamma_Parameters_GeV.h"
#include "NCgamma_Hadronic_Current_R.h"
#include <algorithm>
#include <iostream>

using namespace NC_gamma;

Delta_in_medium::Delta_in_medium() {
    this->mDelta = NCgamma_Param::mDelta;
    this->mDelta2 = NCgamma_Param::mDelta*NCgamma_Param::mDelta;
    this->mn = NCgamma_Param::mp;
    this->mn2 = NCgamma_Param::mp2;
    this->mpi = NCgamma_Param::mpi;
    this->V0 = NCgamma_Param::Delta_V0;
    this->coupling_DeltaNpi = NCgamma_Param::coupling_DeltaNpi;
    this->c_i = NCgamma_Param::c_i;
    this->hc3 = NCgamma_Param::hc * NCgamma_Param::hc * NCgamma_Param::hc;
    this->dir_or_crs = "dir";
    this->rho0 = NCgamma_Param::rho0;

    this->rho_avg = 3.0/(4.0*NCgamma_Param::pi*(1.2*1.2*1.2)) * this->hc3; // average density in fm^-3 * hc^3
    this->rho_r = this->rho_avg ;
}


Delta_in_medium::~Delta_in_medium() {

}


std::complex<double> Delta_in_medium::Delta_propagator_dir(double r) {

    std::complex<double> sigma = this->Sigma();

    std::complex<double> D_Delta_med = 1.0/(this->p2 - (this->mDelta + sigma.real() )*(this->mDelta + sigma.real() )
                                            + this->c_i * (this->mDelta + sigma.real()) * (Gamma_tilde_Delta(this->rho_r) - 2.0 * sigma.imag() ) );
    return D_Delta_med;
}

std::complex<double> Delta_in_medium::Delta_propagator_crs() {
    return 1.0/(this->p2 - this->mDelta2 + this->c_i * this->mDelta * this->Gamma_vacuum(this->p2));
}

void Delta_in_medium::change_other_parameters(std::vector<std::complex<double> > param) {
    this->p2 = param[0].real();
    this->q_minus_kg = param[1].real();
    this->p = sqrt(this->p2);
    this->qcm = sqrt(this->lambda_func(p2, this->mpi*this->mpi, this->mn*this->mn))/(2.0 * this->p);
    this->gamma_vac = this->Gamma_vacuum(this->p2);

    this->propagator_crs = this->Delta_propagator_crs();
    this->propagator_avg = this->Delta_propagator_avg_dir();
}

std::complex<double> Delta_in_medium::Sigma() {
    double reSigma = 0.0;

    double imSigma = - 1.0/2.0 * this->V0 * this->rho_r/this->rho0;

    std::complex<double> Sigma(reSigma, imSigma);
    return Sigma;
}

double Delta_in_medium::Gamma_tilde_Delta(double rho_r_in) {

    double q_tilde = this->qcm/this->kf(rho_r_in);
    return this->gamma_vac * this->I_series(q_tilde);
}

double Delta_in_medium::Gamma_vacuum(double p2_in) {
    if(p2_in > (this->mn + this->mpi)*(this->mn + this->mpi)) {

        return 1.0/(6.0*NCgamma_Param::pi) * (this->coupling_DeltaNpi/this->mpi)*(this->coupling_DeltaNpi/this->mpi)
               *this->mn/sqrt(p2_in) * this->qcm*qcm*qcm;
    }else return 0;
}

double Delta_in_medium::I_series(double q) {
    double res = 1.0;

    if (q != 0) {
        if (q > 1.0) res += -2.0 / (5.0 * q * q) + 9.0 / (35.0 * pow(q, 4)) - 2.0 / (21.0 * pow(q, 6));
        else if (q < 1.0) res += 34.0 / 35.0 * q - 22.0 / 105.0 * q * q * q - 1.0;
    }

    return res;
}

void Delta_in_medium::setNucleon(const std::string &n) {
    Delta_in_medium::nucleon = n;
}

double Delta_in_medium::lambda_func(double x, double y, double z) {
    return x*x + y*y + z*z - 2.0*x*y - 2.0*y*z - 2.0*x*z;
}

double Delta_in_medium::kf(double rho_r_in) {
    return pow(3.0 * NCgamma_Param::pi*NCgamma_Param::pi * rho_r_in/2.0, 1.0/3.0);
}

void Delta_in_medium::setDirOrCrs(const std::string &dirOrCrs) {
    this->dir_or_crs = dirOrCrs;
}

std::complex<double> Delta_in_medium::Delta_propagator_avg_dir() {
    std::complex<double> sigma = this->Sigma();

    std::complex<double> D_Delta_med = 1.0/(this->p2 - (this->mDelta + sigma.real() )*(this->mDelta + sigma.real() )
                                            + this->c_i * (this->mDelta + sigma.real()) * (Gamma_tilde_Delta(this->rho_r) - 2.0 * sigma.imag() ) );
    return D_Delta_med;
}






//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////


Hadronic_Current_R::Hadronic_Current_R() {
    this->c_i = NCgamma_Param::c_i;
}

void Hadronic_Current_R::set_tr_dir(Array4x4 &tr) {}

void Hadronic_Current_R::set_tr_crs(Array4x4 &tr) {}

void Hadronic_Current_R::setQ(const std::vector<double> &q_in) {
    Hadronic_Current_R::q = q_in;
    double q2 = q_in[0] * q_in[0] - q_in[1] * q_in[1] -q_in[2] * q_in[2] -q_in[3] * q_in[3];
    this->setFF(q2);
}

void Hadronic_Current_R::setKg(const std::vector<double> &kg_in) {
    Hadronic_Current_R::kg = kg_in;
}

void Hadronic_Current_R::setP(std::vector<double> p_in) {
    Hadronic_Current_R::p = p_in;
    this->p0 = p_in[0];
    double pp_aux[] = {p_in[0], -p_in[1], -p_in[2], -p_in[3]};
    this->pp.assign(pp_aux, pp_aux + sizeof(pp_aux) / sizeof(double) );
}

std::complex<double> Hadronic_Current_R::getR(int i, int j) {
    return 0;
}

void Hadronic_Current_R::setFF(double q2) {}

Hadronic_Current_R::~Hadronic_Current_R() {}




//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////


Hadronic_Current_R_Sum::Hadronic_Current_R_Sum(const std::vector<Hadronic_Current_R *> &vectorOfCurrents)
        : vector_of_currents(vectorOfCurrents) {}

void Hadronic_Current_R_Sum::setQ(const std::vector<double> &q_in) {
    for (int i = 0; i < static_cast<unsigned int>(this->vector_of_currents.size()); ++i) {
        vector_of_currents[i]->setQ(q_in);
    }
}

void Hadronic_Current_R_Sum::setKg(const std::vector<double> &kg_in) {
    for (int i = 0; i < static_cast<unsigned int>(this->vector_of_currents.size()); ++i) {
        vector_of_currents[i]->setKg(kg_in);
    }
}

void Hadronic_Current_R_Sum::setP(std::vector<double> p_in) {
    for (int i = 0; i < static_cast<unsigned int>(this->vector_of_currents.size()); ++i) {
        vector_of_currents[i]->setP(p_in);
    }
}

std::complex<double> Hadronic_Current_R_Sum::getR(int i, int j) {
    std::complex<double> R  (0.0,0.0);
    for (int k = 0; k < static_cast<unsigned int>(this->vector_of_currents.size()); ++k) {
        R += vector_of_currents[k]->getR(i,j);
    }
    return R;
}

void Hadronic_Current_R_Sum::setFF(double p_in) {
    for (int i = 0; i < static_cast<unsigned int>(this->vector_of_currents.size()); ++i) {
        vector_of_currents[i]->setFF(p_in);
    }
}





//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////



Hadronic_Current_R_Delta::Hadronic_Current_R_Delta(Nuclear_FF *nuclearFf) : nuclearFF(nuclearFf) {
    this->mn = NCgamma_Param::mp;
    this->mn2 = NCgamma_Param::mp2;
    this->mDelta = NCgamma_Param::mDelta;
    this->mDelta2 = NCgamma_Param::mDelta*NCgamma_Param::mDelta;
    this->mpi = NCgamma_Param::mpi;
    this->mpi2 = NCgamma_Param::mpi2;
    this->c_i = NCgamma_Param::c_i;
    this->ff_Delta = new Form_Factors_Delta();

    this->ff_Delta->setFF(0.0);
    this->C3v = this->ff_Delta->getC3v();
    this->C4v = this->ff_Delta->getC4v();
    this->C5v = this->ff_Delta->getC5v();

    this->propagator = new Delta_in_medium();
}

Hadronic_Current_R_Delta::~Hadronic_Current_R_Delta() {
    delete this->ff_Delta;
    delete propagator;
}


void Hadronic_Current_R_Delta::setFF(double q2) {
    double Q2 = -q2;
    this->ff_Delta->setFF(Q2);
}

void Hadronic_Current_R_Delta::setKg(const std::vector<double> &kg_in) {
    this->kg = kg_in;

    this->C3vNC = this->ff_Delta->getC3vNC();
    this->C4vNC = this->ff_Delta->getC4vNC();
    this->C5vNC = this->ff_Delta->getC5vNC();
    this->C3aNC = this->ff_Delta->getC3aNC();
    this->C4aNC = this->ff_Delta->getC4aNC();
    this->C5aNC = this->ff_Delta->getC5aNC();

    this->set_tr_dir(tr_p_dir);
    this->set_tr_crs(tr_p_crs);

    this->set_tr_dir(tr_n_dir);
    this->set_tr_crs(tr_n_crs);

    std::vector<double> pm = this->q;
    std::transform(pm.begin( ), pm.end( ), this->kg.begin( ), pm.begin( ),std::minus<double>( ));
    this->q_minus_kg = sqrt(pm[1]*pm[1]+pm[2]*pm[2]+pm[3]*pm[3]);

    this->nuclearFF->setFF(this->q_minus_kg);
    this->N_FF_p = this->nuclearFF->getFfP();
    this->N_FF_n = this->nuclearFF->getFfN();

    this->dir_or_crs = "dir";
    Propagator_Delta(this->p_dir);
    this->dir_or_crs = "crs";
    Propagator_Delta(this->p_crs);
}


void Hadronic_Current_R_Delta::set_tr_dir(Array4x4 &tr) {

    tr[0][0] = -(C3v*C3vNC*(mDelta*p0 + this->mn*(p0 + this->q[0]))*(4*mDelta2*this->q[1]*this->kg[1] - 4*pow(p0,2)*this->q[1]*this->kg[1] + pow(this->q[1],3)*this->kg[1] + this->q[1]*this->q[3] * this->q[3]*this->kg[1] +
                                                                     this->q[1] * this->q[1]*this->kg[1] * this->kg[1] - this->q[3] * this->q[3]*this->kg[1] * this->kg[1] - this->q[3] * this->q[3]*this->kg[2] * this->kg[2] +
                                                                     this->q[0] * this->q[0]*(this->q[1] * this->q[1] - 3*this->q[1]*this->kg[1] + this->q[3]*(2*this->q[3] - 3*this->kg[3])) +
                                                                     this->q[3]*(4*mDelta2 - 4*pow(p0,2) + this->q[1] * this->q[1] + this->q[3] * this->q[3] + 2*this->q[1]*this->kg[1])*this->kg[3] - 8*p0*this->q[0]*(this->q[1]*this->kg[1] + this->q[3]*this->kg[3])))/
               (3.*mDelta2*pow(this->mn,3));

    tr[0][1] = (this->c_i*C3v*C5aNC*mDelta*(mDelta*p0 + this->mn*(p0 + this->q[0]))*this->q[3]*this->kg[2]*
                (8*this->mn2 - 2*this->q[0]*(p0 + this->q[0]) + this->q[1]*(this->q[1] + this->kg[1]) + this->q[3]*(this->q[3] + this->kg[3])) -
                C3v*C3vNC*this->mn*(16*mDelta2*this->mn*this->q[0]*(p0 + this->q[0])*this->kg[1] + 4*pow(mDelta,3)*(4*p0*this->q[0]*this->kg[1] + this->q[3]*(-(this->q[3]*this->kg[1]) + this->q[1]*this->kg[3])) -
                                    4*this->mn*(p0 + this->q[0])*(4*pow(p0,2)*this->q[0]*this->kg[1] + 8*p0*this->q[0] * this->q[0]*this->kg[1] + 4*pow(this->q[0],3)*this->kg[1] - this->q[0]*this->q[1] * this->q[1]*this->kg[1] +
                                                                  2*p0*this->q[3]*(-(this->q[3]*this->kg[1]) + this->q[1]*this->kg[3]) - this->q[0]*this->q[1]*(2*this->kg[1] * this->kg[1] + this->kg[2] * this->kg[2] + this->kg[3]*(-this->q[3] + this->kg[3])) -
                                                                  this->q[0]*this->kg[1]*(this->q[0] * this->q[0] + this->q[3]*(2*this->q[3] + this->kg[3]))) +
                                    mDelta*(-16*pow(p0,3)*this->q[0]*this->kg[1] - this->q[3]*(this->q[3]*this->kg[1] - this->q[1]*this->kg[3])*(-3*this->q[0] * this->q[0] + this->q[1] * this->q[1] + 2*this->q[1]*this->kg[1] + this->q[3]*(this->q[3] + 2*this->kg[3])) -
                                            4*pow(p0,2)*(8*this->q[0] * this->q[0]*this->kg[1] + 3*this->q[3]*(-(this->q[3]*this->kg[1]) + this->q[1]*this->kg[3])) +
                                            4*p0*this->q[0]*(this->q[0] * this->q[0]*(this->q[1] - 3*this->kg[1]) + this->q[1] * this->q[1]*this->kg[1] + this->q[3]*this->kg[1]*(4*this->q[3] + this->kg[3]) + this->q[1]*(this->kg[1] * this->kg[1] - 3*this->q[3]*this->kg[3])))))/
               (12.*mDelta2*pow(this->mn,4));

    tr[0][2] = (C3v*(std::complex<double>(0,-1)*C5aNC*mDelta*(mDelta*p0 + this->mn*(p0 + this->q[0]))*(this->q[3]*this->kg[1] - this->q[1]*this->kg[3])*
                     (8*this->mn2 - 2*this->q[0]*(p0 + this->q[0]) + this->q[1]*(this->q[1] + this->kg[1]) + this->q[3]*(this->q[3] + this->kg[3])) +
                     C3vNC*this->mn*this->kg[2]*(-16*mDelta2*this->mn*this->q[0]*(p0 + this->q[0]) + 4*pow(mDelta,3)*(-4*p0*this->q[0] + this->q[1] * this->q[1] + this->q[3] * this->q[3]) +
                                                 4*this->mn*(p0 + this->q[0])*(4*pow(p0,2)*this->q[0] + 8*p0*this->q[0] * this->q[0] + 3*pow(this->q[0],3) - 2*p0*(this->q[1] * this->q[1] + this->q[3] * this->q[3]) -
                                                                               this->q[0]*this->q[1]*(2*this->q[1] + this->kg[1]) - this->q[0]*this->q[3]*(2*this->q[3] + this->kg[3])) +
                                                 mDelta*(16*pow(p0,3)*this->q[0] + 4*pow(p0,2)*(8*this->q[0] * this->q[0] - 3*(this->q[1] * this->q[1] + this->q[3] * this->q[3])) +
                                                         (this->q[1] * this->q[1] + this->q[3] * this->q[3])*(-3*this->q[0] * this->q[0] + this->q[1] * this->q[1] + this->q[3] * this->q[3] + 2*this->q[1]*this->kg[1] + 2*this->q[3]*this->kg[3]) +
                                                         4*p0*this->q[0]*(3*this->q[0] * this->q[0] - this->q[1]*(4*this->q[1] + this->kg[1]) - this->q[3]*(4*this->q[3] + this->kg[3]))))))/(12.*mDelta2*pow(this->mn,4));

    tr[0][3] = (C3v*(std::complex<double>(0,-1)*C5aNC*mDelta*(mDelta*p0 + this->mn*(p0 + this->q[0]))*this->q[1]*this->kg[2]*
                     (8*this->mn2 - 2*this->q[0]*(p0 + this->q[0]) + this->q[1]*(this->q[1] + this->kg[1]) + this->q[3]*(this->q[3] + this->kg[3])) +
                     C3vNC*this->mn*(-16*mDelta2*this->mn*this->q[0]*(p0 + this->q[0])*this->kg[3] + 4*pow(mDelta,3)*(-(this->q[1]*this->q[3]*this->kg[1]) - 4*p0*this->q[0]*this->kg[3] + this->q[1] * this->q[1]*this->kg[3]) +
                                     4*this->mn*(p0 + this->q[0])*(this->q[0]*this->q[3]*(this->kg[1]*(this->q[1] + this->kg[1]) + this->kg[2] * this->kg[2]) + 8*p0*this->q[0] * this->q[0]*this->kg[3] -
                                                                   this->q[0]*(-4*pow(p0,2) + this->q[3] * this->q[3] + this->q[1]*(2*this->q[1] + this->kg[1]))*this->kg[3] + pow(this->q[0],3)*(-2*this->q[3] + 3*this->kg[3]) + 2*p0*this->q[1]*(this->q[3]*this->kg[1] - this->q[1]*this->kg[3]))\
           + mDelta*(16*pow(p0,3)*this->q[0]*this->kg[3] + 4*pow(p0,2)*(3*this->q[1]*this->q[3]*this->kg[1] + 8*this->q[0] * this->q[0]*this->kg[3] - 3*this->q[1] * this->q[1]*this->kg[3]) +
                     this->q[1]*(-(this->q[3]*this->kg[1]) + this->q[1]*this->kg[3])*(-3*this->q[0] * this->q[0] + this->q[1] * this->q[1] + this->q[3] * this->q[3] + 2*this->q[1]*this->kg[1] + 2*this->q[3]*this->kg[3]) -
                     4*p0*this->q[0]*(-(this->q[3]*(3*this->q[1]*this->kg[1] + this->kg[1] * this->kg[1] + this->kg[2] * this->kg[2])) + this->q[0] * this->q[0]*(2*this->q[3] - 3*this->kg[3]) +
                                      (this->q[3] * this->q[3] + this->q[1]*(4*this->q[1] + this->kg[1]))*this->kg[3])))))/(12.*mDelta2*pow(this->mn,4));


    tr[1][0] = (C3v*(this->c_i*C5aNC*this->q[3]*this->kg[2]*(8*pow(mDelta,3)*this->mn*(p0 + this->q[0]) -
                                                             4*this->mn2*(p0 + this->q[0])*(2*p0*this->q[0] + this->q[0] * this->q[0] - this->q[1]*this->kg[1] - this->q[3]*this->kg[3]) -
                                                             2*mDelta*this->mn*(p0 + this->q[0])*(8*this->mn2 + 2*p0*this->q[0] + this->q[0] * this->q[0] - this->q[1]*this->kg[1] - this->q[3]*this->kg[3]) +
                                                             mDelta2*(8*pow(p0,3) + 12*pow(p0,2)*this->q[0] + 4*p0*this->q[0] * this->q[0] + 8*this->mn2*(-2*p0 + this->q[0]) +
                                                                      this->q[0]*(this->q[1] * this->q[1] - this->q[1]*this->kg[1] + this->q[3]*(this->q[3] - this->kg[3])) - 2*p0*(this->q[1] * this->q[1] + 2*this->q[1]*this->kg[1] + this->q[3]*(this->q[3] + 2*this->kg[3])))) -
                     2*C3vNC*this->mn*(16*mDelta2*this->mn*this->q[0]*(p0 + this->q[0])*this->q[1] +
                                       4*pow(mDelta,3)*(this->q[1]*(4*p0*this->q[0] - this->q[0] * this->q[0] + this->kg[1] * this->kg[1]) + this->q[3]*this->kg[1]*this->kg[3]) +
                                       mDelta*(-16*pow(p0,3)*this->q[0]*this->q[1] + 4*p0*this->q[0]*(pow(this->q[1],3) + 2*this->q[1] * this->q[1]*this->kg[1] - 3*this->q[1]*this->kg[1] * this->kg[1] + this->q[3]*this->kg[1]*(this->q[3] - 3*this->kg[3]) +
                                                                                                      this->q[1]*this->q[3]*(this->q[3] + this->kg[3])) + (3*this->q[0] * this->q[0] - this->q[1]*(this->q[1] + 2*this->kg[1]) - this->q[3]*(this->q[3] + 2*this->kg[3]))*(this->q[0] * this->q[0]*this->q[1] - this->kg[1]*(this->q[1]*this->kg[1] + this->q[3]*this->kg[3])) -
                                               4*pow(p0,2)*(5*this->q[0] * this->q[0]*this->q[1] + 3*this->kg[1]*(this->q[1]*this->kg[1] + this->q[3]*this->kg[3]))) -
                                       4*this->mn*(p0 + this->q[0])*(4*pow(p0,2)*this->q[0]*this->q[1] + 6*p0*this->q[0] * this->q[0]*this->q[1] + 2*p0*this->kg[1]*(this->q[1]*this->kg[1] + this->q[3]*this->kg[3]) +
                                                                     this->q[0]*(2*this->q[0] * this->q[0]*this->q[1] - pow(this->q[1],3) - 2*this->q[1] * this->q[1]*this->kg[1] + this->q[3]*this->kg[1]*(-this->q[3] + this->kg[3]) + this->q[1]*(this->kg[1] * this->kg[1] - this->q[3]*(this->q[3] + this->kg[3])))))))/
               (24.*mDelta2*pow(this->mn,4));

    tr[1][1] = (C3v*(this->c_i*C5aNC*this->q[3]*this->kg[2]*(4*pow(mDelta,3)*this->mn*(this->q[1] + this->kg[1]) + 2*this->mn2*(this->q[1] + this->kg[1])*(-(this->q[0]*(2*p0 + this->q[0])) + this->q[1]*this->kg[1] + this->q[3]*this->kg[3]) +
                                                             mDelta*this->mn*(this->q[1] + this->kg[1])*(-4*p0*this->q[0] - 3*this->q[0] * this->q[0] + this->q[1] * this->q[1] + this->q[3] * this->q[3] + 2*this->q[1]*this->kg[1] + 2*this->q[3]*this->kg[3]) +
                                                             mDelta2*(this->q[0] * this->q[0]*(this->q[1] - this->kg[1]) + 8*this->mn2*this->kg[1] - 2*p0*this->q[0]*this->kg[1] + this->q[3] * this->q[3]*this->kg[1] +
                                                                      4*pow(p0,2)*(this->q[1] + this->kg[1]) - this->q[1]*this->q[3]*this->kg[3])) -
                     2*C3vNC*this->mn*(16*mDelta2*this->mn*(p0 + this->q[0])*(this->q[0] * this->q[0] - this->q[3]*this->kg[3]) +
                                       4*pow(mDelta,3)*(-(this->q[0]*(-4*p0*this->q[0] + this->q[0] * this->q[0] + this->q[3] * this->q[3] - this->kg[1] * this->kg[1])) + 2*(-2*p0 + this->q[0])*this->q[3]*this->kg[3]) -
                                       4*this->mn*(p0 + this->q[0])*(4*pow(p0,2)*this->q[0] * this->q[0] + 2*pow(this->q[0],4) -
                                                                     this->q[0] * this->q[0]*(this->q[1] * this->q[1] + 2*this->q[1]*this->kg[1] - this->kg[1] * this->kg[1] + 3*this->q[3]*this->kg[3]) +
                                                                     this->q[3]*(-(this->q[3]*(2*this->kg[1] * this->kg[1] + this->kg[2] * this->kg[2])) + (this->q[1] * this->q[1] + this->q[3] * this->q[3] + 2*this->q[1]*this->kg[1])*this->kg[3]) +
                                                                     2*p0*this->q[0]*(3*this->q[0] * this->q[0] + this->kg[1] * this->kg[1] - this->q[3]*(this->q[3] + 2*this->kg[3]))) -
                                       mDelta*(16*pow(p0,3)*this->q[0] * this->q[0] + 4*pow(p0,2)*this->q[0]*(5*this->q[0] * this->q[0] - 3*this->q[3] * this->q[3] + 3*this->kg[1] * this->kg[1] - 2*this->q[3]*this->kg[3]) +
                                               4*p0*(-2*this->q[3] * this->q[3]*this->kg[1] * this->kg[1] - this->q[0] * this->q[0]*(2*this->q[3] * this->q[3] + (this->q[1] - this->kg[1])*(this->q[1] + 3*this->kg[1])) +
                                                     this->q[3]*(this->q[0] * this->q[0] + this->q[1] * this->q[1] + this->q[3] * this->q[3] + 2*this->q[1]*this->kg[1])*this->kg[3]) -
                                               this->q[0]*(this->q[0] * this->q[0] + this->q[3] * this->q[3] - this->kg[1] * this->kg[1] - 2*this->q[3]*this->kg[3])*(3*this->q[0] * this->q[0] - this->q[1]*(this->q[1] + 2*this->kg[1]) - this->q[3]*(this->q[3] + 2*this->kg[3]))))))/
               (24.*mDelta2*pow(this->mn,4));

    tr[1][2] = (C3v*(-2*C3vNC*this->mn*this->kg[2]*(16*mDelta2*this->mn*(p0 + this->q[0])*this->q[1] + 4*pow(mDelta,3)*(4*p0*this->q[1] + this->q[0]*(-2*this->q[1] + this->kg[1])) +
                                                    4*this->mn*(p0 + this->q[0])*(pow(this->q[1],3) + this->q[1]*this->q[3] * this->q[3] + 2*this->q[1] * this->q[1]*this->kg[1] + this->q[3] * this->q[3]*this->kg[1] - 2*p0*this->q[0]*(2*this->q[1] + this->kg[1]) -
                                                                                  this->q[0] * this->q[0]*(2*this->q[1] + this->kg[1]) + this->q[1]*this->q[3]*this->kg[3]) +
                                                    mDelta*(-4*pow(p0,2)*this->q[0]*(2*this->q[1] + 3*this->kg[1]) +
                                                            4*p0*(this->q[0] * this->q[0]*(2*this->q[1] - 3*this->kg[1]) + (this->q[1] * this->q[1] + this->q[3] * this->q[3])*(this->q[1] + 2*this->kg[1])) -
                                                            this->q[0]*(2*this->q[1] - this->kg[1])*(-3*this->q[0] * this->q[0] + this->q[1] * this->q[1] + this->q[3] * this->q[3] + 2*this->q[1]*this->kg[1] + 2*this->q[3]*this->kg[3]))) +
                     this->c_i*C5aNC*(2*this->mn2*this->q[3]*this->kg[2] * this->kg[2]*(-(this->q[0]*(2*p0 + this->q[0])) + this->q[1]*this->kg[1] + this->q[3]*this->kg[3]) +
                                      4*pow(mDelta,3)*this->mn*(-(this->q[3]*this->kg[1] * this->kg[1]) + this->q[0] * this->q[0]*(this->q[3] - 2*this->kg[3]) - 2*p0*this->q[0]*this->kg[3] +
                                                                (8*this->mn2 + this->q[3] * this->q[3] + this->q[1]*(this->q[1] + this->kg[1]))*this->kg[3]) +
                                      mDelta2*(-2*p0*this->q[3]*(2*p0*(-2*this->q[0] * this->q[0] + this->kg[1] * this->kg[1]) +
                                                                 this->q[0]*(-3*this->q[0] * this->q[0] + this->q[1] * this->q[1] + this->q[3] * this->q[3] + this->q[1]*this->kg[1] + this->kg[1] * this->kg[1])) +
                                               this->q[3]*(-(this->q[0]*(4*p0 + this->q[0])) + this->q[1] * this->q[1] + this->q[3] * this->q[3])*this->kg[2] * this->kg[2] +
                                               2*p0*(-4*pow(p0,2)*this->q[0] - 2*pow(this->q[0],3) + this->q[0]*this->q[1]*(this->q[1] + this->kg[1]) +
                                                     2*p0*(-3*this->q[0] * this->q[0] + this->q[3] * this->q[3] + this->q[1]*(this->q[1] + this->kg[1])))*this->kg[3] +
                                               8*this->mn2*(this->q[3]*(-2*p0*this->q[0] + this->kg[2] * this->kg[2]) + 2*p0*(2*p0 + this->q[0])*this->kg[3])) +
                                      mDelta*this->mn*(this->q[3]*(-pow(this->q[0],4) - this->kg[1] * this->kg[1]*(3*(this->q[1] * this->q[1] + this->q[3] * this->q[3]) + 4*this->q[1]*this->kg[1]) -
                                                                   2*(this->q[3] * this->q[3] + this->q[1]*(this->q[1] + this->kg[1]))*this->kg[2] * this->kg[2] +
                                                                   this->q[0] * this->q[0]*(this->q[1] * this->q[1] + this->q[3] * this->q[3] + 2*this->q[1]*this->kg[1] + 5*this->kg[1] * this->kg[1] + 2*this->kg[2] * this->kg[2])) +
                                                       (2*pow(this->q[0],4) + pow(this->q[1] * this->q[1] + this->q[3] * this->q[3],2) + 3*this->q[1]*(this->q[1] * this->q[1] + this->q[3] * this->q[3])*this->kg[1] +
                                                        2*(this->q[1] - this->q[3])*(this->q[1] + this->q[3])*this->kg[1] * this->kg[1] - this->q[0] * this->q[0]*(3*(this->q[1] * this->q[1] + this->q[3] * this->q[3]) + 5*this->q[1]*this->kg[1]))*this->kg[3] +
                                                       4*pow(p0,2)*this->q[0] * this->q[0]*(this->q[3] + this->kg[3]) -
                                                       8*this->mn2*(2*p0*this->q[0]*(this->q[3] + this->kg[3]) + this->q[0] * this->q[0]*(2*this->q[3] + this->kg[3]) -
                                                                    this->kg[3]*(this->q[1] * this->q[1] + this->q[3] * this->q[3] + 2*this->q[1]*this->kg[1] + 2*this->q[3]*this->kg[3])) +
                                                       2*p0*this->q[0]*(-(this->q[1] * this->q[1]*(this->q[3] + 2*this->kg[3])) + this->q[0] * this->q[0]*(this->q[3] + 3*this->kg[3]) - this->q[1]*this->kg[1]*(this->q[3] + 3*this->kg[3]) +
                                                                        this->q[3]*(3*this->kg[1] * this->kg[1] + this->kg[2] * this->kg[2] - this->q[3]*(this->q[3] + 3*this->kg[3])))))))/(24.*mDelta2*pow(this->mn,4));

    tr[1][3] = -(C3v*(this->c_i*C5aNC*this->kg[2]*(4*pow(mDelta,3)*this->mn*(8*this->mn2 - 2*this->q[0]*(p0 + this->q[0]) + this->q[1]*(this->q[1] + this->kg[1])) -
                                                   2*this->mn2*this->q[3]*(this->q[3] + this->kg[3])*(-(this->q[0]*(2*p0 + this->q[0])) + this->q[1]*this->kg[1] + this->q[3]*this->kg[3]) +
                                                   mDelta*this->mn*(4*pow(p0,2)*this->q[0] * this->q[0] + 6*p0*pow(this->q[0],3) + 2*pow(this->q[0],4) - 4*p0*this->q[0]*this->q[1] * this->q[1] -
                                                                    3*this->q[0] * this->q[0]*this->q[1] * this->q[1] + pow(this->q[1],4) + this->q[1] * this->q[1]*this->q[3] * this->q[3] - 6*p0*this->q[0]*this->q[1]*this->kg[1] - 5*this->q[0] * this->q[0]*this->q[1]*this->kg[1] +
                                                                    3*pow(this->q[1],3)*this->kg[1] + this->q[1]*this->q[3] * this->q[3]*this->kg[1] + 2*this->q[1] * this->q[1]*this->kg[1] * this->kg[1] + 2*this->q[3]*(-(this->q[0]*(p0 + this->q[0])) + this->q[1]*(this->q[1] + this->kg[1]))*this->kg[3] +
                                                                    8*this->mn2*(-(this->q[0]*(2*p0 + this->q[0])) + this->q[1] * this->q[1] + 2*this->q[1]*this->kg[1] + this->q[3]*(this->q[3] + 2*this->kg[3]))) +
                                                   mDelta2*(-8*pow(p0,3)*this->q[0] + 4*pow(p0,2)*(-3*this->q[0] * this->q[0] + this->q[1]*(this->q[1] + this->kg[1])) +
                                                            8*this->mn2*(2*p0*(2*p0 + this->q[0]) - this->q[3]*this->kg[3]) + 2*p0*this->q[0]*(-2*this->q[0] * this->q[0] + this->q[1]*(this->q[1] + this->kg[1]) + this->q[3]*(this->q[3] + 2*this->kg[3])) +
                                                            this->q[3]*(this->q[0] * this->q[0]*(-this->q[3] + this->kg[3]) + this->q[1]*(this->q[3]*this->kg[1] - this->q[1]*this->kg[3])))) +
                      2*C3vNC*this->mn*(16*mDelta2*this->mn*(p0 + this->q[0])*this->q[1]*this->kg[3] + 4*pow(mDelta,3)*(this->q[0]*this->q[1]*(this->q[3] - 2*this->kg[3]) + 4*p0*this->q[1]*this->kg[3] + this->q[0]*this->kg[1]*this->kg[3]) -
                                        4*this->mn*(p0 + this->q[0])*(this->q[1]*this->q[3]*(2*this->kg[1] * this->kg[1] + this->kg[2] * this->kg[2]) - this->q[1]*(this->q[1] * this->q[1] + this->q[3] * this->q[3] + 2*this->q[1]*this->kg[1])*this->kg[3] +
                                                                      this->q[0] * this->q[0]*(-(this->q[3]*(this->q[1] + this->kg[1])) + (2*this->q[1] + this->kg[1])*this->kg[3]) + 2*p0*this->q[0]*(this->kg[1]*this->kg[3] + this->q[1]*(this->q[3] + 2*this->kg[3]))) +
                                        mDelta*(-4*pow(p0,2)*this->q[0]*(3*this->q[1]*this->q[3] + 2*this->q[1]*this->kg[3] + 3*this->kg[1]*this->kg[3]) -
                                                this->q[0]*(this->q[1]*(this->q[3] - 2*this->kg[3]) + this->kg[1]*this->kg[3])*(3*this->q[0] * this->q[0] - this->q[1]*(this->q[1] + 2*this->kg[1]) - this->q[3]*(this->q[3] + 2*this->kg[3])) +
                                                4*p0*(-2*this->q[1]*this->q[3]*this->kg[1] * this->kg[1] + this->q[1]*(this->q[1] * this->q[1] + this->q[3] * this->q[3] + 2*this->q[1]*this->kg[1])*this->kg[3] +
                                                      this->q[0] * this->q[0]*(-(this->q[1]*this->q[3]) + this->q[3]*this->kg[1] + 2*this->q[1]*this->kg[3] - 3*this->kg[1]*this->kg[3]))))))/(24.*mDelta2*pow(this->mn,4));


    tr[2][0] = (C3v*(-2*C3vNC*this->mn*this->kg[2]*(4*pow(mDelta,3)*(this->q[1]*this->kg[1] + this->q[3]*this->kg[3]) -
                                                    4*this->mn*(p0 + this->q[0])*(-(this->q[0]*(this->q[1] * this->q[1] + this->q[3] * this->q[3] - this->q[1]*this->kg[1] - this->q[3]*this->kg[3])) + 2*p0*(this->q[1]*this->kg[1] + this->q[3]*this->kg[3])) +
                                                    mDelta*(4*p0*this->q[0]*(this->q[1] * this->q[1] + this->q[3] * this->q[3] - 3*this->q[1]*this->kg[1] - 3*this->q[3]*this->kg[3]) - 12*pow(p0,2)*(this->q[1]*this->kg[1] + this->q[3]*this->kg[3]) +
                                                            (this->q[1]*this->kg[1] + this->q[3]*this->kg[3])*(-3*this->q[0] * this->q[0] + this->q[1] * this->q[1] + this->q[3] * this->q[3] + 2*this->q[1]*this->kg[1] + 2*this->q[3]*this->kg[3]))) -
                     this->c_i*C5aNC*(this->q[3]*this->kg[1] - this->q[1]*this->kg[3])*(8*pow(mDelta,3)*this->mn*(p0 + this->q[0]) -
                                                                                        4*this->mn2*(p0 + this->q[0])*(2*p0*this->q[0] + this->q[0] * this->q[0] - this->q[1]*this->kg[1] - this->q[3]*this->kg[3]) -
                                                                                        2*mDelta*this->mn*(p0 + this->q[0])*(8*this->mn2 + 2*p0*this->q[0] + this->q[0] * this->q[0] - this->q[1]*this->kg[1] - this->q[3]*this->kg[3]) +
                                                                                        mDelta2*(8*pow(p0,3) + 12*pow(p0,2)*this->q[0] + 4*p0*this->q[0] * this->q[0] + 8*this->mn2*(-2*p0 + this->q[0]) +
                                                                                                 this->q[0]*(this->q[1] * this->q[1] - this->q[1]*this->kg[1] + this->q[3]*(this->q[3] - this->kg[3])) - 2*p0*(this->q[1] * this->q[1] + 2*this->q[1]*this->kg[1] + this->q[3]*(this->q[3] + 2*this->kg[3]))))))/
               (24.*mDelta2*pow(this->mn,4));

    tr[2][1] = (C3v*(-2*C3vNC*this->mn*this->kg[2]*(4*pow(mDelta,3)*this->q[0]*this->kg[1] - 4*this->mn*(p0 + this->q[0])*
                                                                                             (2*p0*this->q[0]*this->kg[1] + this->q[0] * this->q[0]*(-this->q[1] + this->kg[1]) + this->q[3]*(-(this->q[3]*this->kg[1]) + this->q[1]*this->kg[3])) +
                                                    mDelta*(4*p0*this->q[0] * this->q[0]*(this->q[1] - 3*this->kg[1]) - 12*pow(p0,2)*this->q[0]*this->kg[1] + 8*p0*this->q[3]*(this->q[3]*this->kg[1] - this->q[1]*this->kg[3]) +
                                                            this->q[0]*this->kg[1]*(-3*this->q[0] * this->q[0] + this->q[1] * this->q[1] + this->q[3] * this->q[3] + 2*this->q[1]*this->kg[1] + 2*this->q[3]*this->kg[3]))) -
                     this->c_i*C5aNC*(2*this->mn2*(this->q[1] + this->kg[1])*(-(this->q[3]*this->kg[1]) + this->q[1]*this->kg[3])*(2*p0*this->q[0] + this->q[0] * this->q[0] - this->q[1]*this->kg[1] - this->q[3]*this->kg[3]) +
                                      4*pow(mDelta,3)*this->mn*(this->q[3]*(this->q[0] * this->q[0] + this->q[1]*this->kg[1] - this->kg[2] * this->kg[2]) + (8*this->mn2 - 2*this->q[0]*(p0 + this->q[0]) + this->q[3] * this->q[3])*this->kg[3]) +
                                      mDelta2*(-8*pow(p0,3)*this->q[0]*this->kg[3] + 2*p0*this->q[0]*
                                                                                     (-(this->q[3]*(this->q[1] * this->q[1] + this->q[3] * this->q[3] + this->q[1]*this->kg[1] + 2*this->kg[1] * this->kg[1] + this->kg[2] * this->kg[2])) + this->q[0] * this->q[0]*(3*this->q[3] - 2*this->kg[3]) +
                                                                                      this->q[1]*(this->q[1] + 2*this->kg[1])*this->kg[3]) + 4*pow(p0,2)*(this->q[3]*this->kg[1]*(this->q[1] + this->kg[1]) + this->q[0] * this->q[0]*(this->q[3] - 3*this->kg[3]) + this->q[3]*this->kg[3]*(this->q[3] + this->kg[3])) +
                                               (this->q[3]*this->kg[1] - this->q[1]*this->kg[3])*(this->q[0] * this->q[0]*(this->q[1] - this->kg[1]) + this->q[3]*(this->q[3]*this->kg[1] - this->q[1]*this->kg[3])) +
                                               8*this->mn2*(4*pow(p0,2)*this->kg[3] + 2*p0*this->q[0]*(-this->q[3] + this->kg[3]) + this->kg[1]*(this->q[3]*this->kg[1] - this->q[1]*this->kg[3]))) +
                                      mDelta*this->mn*(this->q[3]*(-pow(this->q[0],4) + this->kg[1]*(pow(this->q[1],3) + this->q[1]*this->q[3] * this->q[3] + 2*this->q[1] * this->q[1]*this->kg[1] - 2*this->q[3] * this->q[3]*this->kg[1]) -
                                                                   (this->q[1] * this->q[1] + 3*this->q[3] * this->q[3] + 2*this->q[1]*this->kg[1])*this->kg[2] * this->kg[2] +
                                                                   this->q[0] * this->q[0]*(-this->q[1] * this->q[1] + this->q[3] * this->q[3] - 3*this->q[1]*this->kg[1] + 2*this->kg[1] * this->kg[1] + 5*this->kg[2] * this->kg[2])) +
                                                       (2*pow(this->q[0],4) - this->q[0] * this->q[0]*(3*this->q[3] * this->q[3] + 2*this->q[1]*this->kg[1]) +
                                                        this->q[3] * this->q[3]*(this->q[1] * this->q[1] + this->q[3] * this->q[3] + 4*this->q[1]*this->kg[1] - 2*this->kg[2] * this->kg[2]))*this->kg[3] + 4*pow(p0,2)*this->q[0] * this->q[0]*(this->q[3] + this->kg[3]) -
                                                       8*this->mn2*(2*p0*this->q[0]*(this->q[3] + this->kg[3]) + this->q[0] * this->q[0]*(2*this->q[3] + this->kg[3]) -
                                                                    this->kg[3]*(this->q[1] * this->q[1] + this->q[3] * this->q[3] + 2*this->q[1]*this->kg[1] + 2*this->q[3]*this->kg[3])) +
                                                       2*p0*this->q[0]*(-(this->q[1] * this->q[1]*this->q[3]) - this->q[1]*this->kg[1]*(3*this->q[3] + this->kg[3]) + this->q[0] * this->q[0]*(this->q[3] + 3*this->kg[3]) +
                                                                        this->q[3]*(this->kg[1] * this->kg[1] + 3*this->kg[2] * this->kg[2] - this->q[3]*(this->q[3] + 3*this->kg[3])))))))/(24.*mDelta2*pow(this->mn,4));

    tr[2][2] = (C3v*(std::complex<double>(0,-1)*C5aNC*this->kg[2]*(this->q[3]*this->kg[1] - this->q[1]*this->kg[3])*(4*pow(mDelta,3)*this->mn +
                                                                                                                     mDelta2*(8*this->mn2 + 4*pow(p0,2) - 2*p0*this->q[0] - this->q[0] * this->q[0] + this->q[1] * this->q[1] + this->q[3] * this->q[3]) +
                                                                                                                     2*this->mn2*(-(this->q[0]*(2*p0 + this->q[0])) + this->q[1]*this->kg[1] + this->q[3]*this->kg[3]) +
                                                                                                                     mDelta*this->mn*(-4*p0*this->q[0] - 3*this->q[0] * this->q[0] + this->q[1] * this->q[1] + this->q[3] * this->q[3] + 2*this->q[1]*this->kg[1] + 2*this->q[3]*this->kg[3])) +
                     2*C3vNC*this->mn*(-16*mDelta2*this->mn*(p0 + this->q[0])*(this->q[0] * this->q[0] - this->q[1]*this->kg[1] - this->q[3]*this->kg[3]) +
                                       4*pow(mDelta,3)*(this->q[0]*(this->q[0] * this->q[0] + this->q[1] * this->q[1] + this->q[3] * this->q[3] - 2*this->q[1]*this->kg[1] - this->kg[2] * this->kg[2] - 2*this->q[3]*this->kg[3]) +
                                                        4*p0*(-this->q[0] * this->q[0] + this->q[1]*this->kg[1] + this->q[3]*this->kg[3])) +
                                       4*this->mn*(p0 + this->q[0])*(4*pow(p0,2)*this->q[0] * this->q[0] + 2*pow(this->q[0],4) + pow(this->q[1],3)*this->kg[1] + this->q[1]*this->q[3] * this->q[3]*this->kg[1] +
                                                                     this->q[1] * this->q[1]*this->kg[1] * this->kg[1] - this->q[3] * this->q[3]*this->kg[1] * this->kg[1] - this->q[1] * this->q[1]*this->kg[2] * this->kg[2] - 2*this->q[3] * this->q[3]*this->kg[2] * this->kg[2] +
                                                                     this->q[3]*(this->q[1] * this->q[1] + this->q[3] * this->q[3] + 2*this->q[1]*this->kg[1])*this->kg[3] - this->q[0] * this->q[0]*(this->q[1] * this->q[1] + 3*this->q[1]*this->kg[1] - this->kg[2] * this->kg[2] + 3*this->q[3]*this->kg[3]) +
                                                                     2*p0*this->q[0]*(3*this->q[0] * this->q[0] - this->q[1]*(this->q[1] + 2*this->kg[1]) + this->kg[2] * this->kg[2] - this->q[3]*(this->q[3] + 2*this->kg[3]))) +
                                       mDelta*(16*pow(p0,3)*this->q[0] * this->q[0] + 4*pow(p0,2)*this->q[0]*
                                                                                      (5*this->q[0] * this->q[0] - 3*this->q[1] * this->q[1] - 3*this->q[3] * this->q[3] - 2*this->q[1]*this->kg[1] + 3*this->kg[2] * this->kg[2] - 2*this->q[3]*this->kg[3]) -
                                               this->q[0]*(this->q[0] * this->q[0] + this->q[1] * this->q[1] + this->q[3] * this->q[3] - 2*this->q[1]*this->kg[1] - this->kg[2] * this->kg[2] - 2*this->q[3]*this->kg[3])*
                                               (3*this->q[0] * this->q[0] - this->q[1]*(this->q[1] + 2*this->kg[1]) - this->q[3]*(this->q[3] + 2*this->kg[3])) +
                                               4*p0*(pow(this->q[1],3)*this->kg[1] + this->q[1]*this->q[3] * this->q[3]*this->kg[1] + this->q[3] * this->q[3]*(-2*this->kg[2] * this->kg[2] + this->q[3]*this->kg[3]) +
                                                     this->q[0] * this->q[0]*(this->q[1]*(-4*this->q[1] + this->kg[1]) + 3*this->kg[2] * this->kg[2] + this->q[3]*(-2*this->q[3] + this->kg[3])) + this->q[1] * this->q[1]*(2*this->kg[1] * this->kg[1] + this->kg[3]*(this->q[3] + 2*this->kg[3]))
                                               )))))/(24.*mDelta2*pow(this->mn,4));

    tr[2][3] = (C3v*(-2*C3vNC*this->mn*this->kg[2]*(4*pow(mDelta,3)*this->q[0]*this->kg[3] - 4*this->mn*(p0 + this->q[0])*(2*p0*this->q[0]*this->kg[3] + this->q[0] * this->q[0]*(-this->q[3] + this->kg[3]) + this->q[1]*(this->q[3]*this->kg[1] - this->q[1]*this->kg[3])) +
                                                    mDelta*(4*p0*this->q[0] * this->q[0]*(this->q[3] - 3*this->kg[3]) - 12*pow(p0,2)*this->q[0]*this->kg[3] + 8*p0*this->q[1]*(-(this->q[3]*this->kg[1]) + this->q[1]*this->kg[3]) +
                                                            this->q[0]*this->kg[3]*(-3*this->q[0] * this->q[0] + this->q[1] * this->q[1] + this->q[3] * this->q[3] + 2*this->q[1]*this->kg[1] + 2*this->q[3]*this->kg[3]))) +
                     this->c_i*C5aNC*(-2*this->mn2*(this->q[3] + this->kg[3])*(this->q[3]*this->kg[1] - this->q[1]*this->kg[3])*(-(this->q[0]*(2*p0 + this->q[0])) + this->q[1]*this->kg[1] + this->q[3]*this->kg[3]) +
                                      4*pow(mDelta,3)*this->mn*(this->q[0] * this->q[0]*(this->q[1] - 2*this->kg[1]) + 8*this->mn2*this->kg[1] - 2*p0*this->q[0]*this->kg[1] + this->q[1]*(this->q[1]*this->kg[1] - this->kg[2] * this->kg[2] + this->q[3]*this->kg[3])) +
                                      mDelta2*(-8*pow(p0,3)*this->q[0]*this->kg[1] + 4*pow(p0,2)*
                                                                                     (this->q[0] * this->q[0]*(this->q[1] - 3*this->kg[1]) + this->q[1]*this->kg[1]*(this->q[1] + this->kg[1]) + this->q[1]*this->kg[3]*(this->q[3] + this->kg[3])) +
                                               (this->q[3]*this->kg[1] - this->q[1]*this->kg[3])*(this->q[0] * this->q[0]*(-this->q[3] + this->kg[3]) + this->q[1]*(this->q[3]*this->kg[1] - this->q[1]*this->kg[3])) +
                                               8*this->mn2*(4*pow(p0,2)*this->kg[1] + 2*p0*this->q[0]*(-this->q[1] + this->kg[1]) + this->kg[3]*(-(this->q[3]*this->kg[1]) + this->q[1]*this->kg[3])) +
                                               2*p0*this->q[0]*(-pow(this->q[1],3) + this->q[0] * this->q[0]*(this->q[1] - 2*this->kg[1]) + this->q[3]*this->kg[1]*(this->q[3] + 2*this->kg[3]) +
                                                                this->q[1]*(2*this->kg[1] * this->kg[1] + this->kg[2] * this->kg[2] - this->q[3]*(this->q[3] + this->kg[3])))) +
                                      mDelta*this->mn*(pow(this->q[0],4)*this->q[1] - this->q[0] * this->q[0]*pow(this->q[1],3) + this->q[0] * this->q[0]*this->q[1]*this->q[3] * this->q[3] + 2*pow(this->q[0],4)*this->kg[1] -
                                                       3*this->q[0] * this->q[0]*this->q[1] * this->q[1]*this->kg[1] + pow(this->q[1],4)*this->kg[1] + this->q[1] * this->q[1]*this->q[3] * this->q[3]*this->kg[1] - 2*this->q[0] * this->q[0]*this->q[1]*this->kg[1] * this->kg[1] +
                                                       2*pow(this->q[1],3)*this->kg[1] * this->kg[1] - 2*this->q[1]*this->q[3] * this->q[3]*this->kg[1] * this->kg[1] + 4*pow(p0,2)*this->q[0] * this->q[0]*(this->q[1] + this->kg[1]) +
                                                       3*this->q[0] * this->q[0]*this->q[1]*this->kg[2] * this->kg[2] - pow(this->q[1],3)*this->kg[2] * this->kg[2] - 3*this->q[1]*this->q[3] * this->q[3]*this->kg[2] * this->kg[2] -
                                                       2*this->q[1] * this->q[1]*this->kg[1]*this->kg[2] * this->kg[2] - 3*this->q[0] * this->q[0]*this->q[1]*this->q[3]*this->kg[3] + pow(this->q[1],3)*this->q[3]*this->kg[3] + this->q[1]*pow(this->q[3],3)*this->kg[3] -
                                                       2*this->q[0] * this->q[0]*this->q[3]*this->kg[1]*this->kg[3] + 4*this->q[1] * this->q[1]*this->q[3]*this->kg[1]*this->kg[3] - 2*this->q[1]*this->q[3]*this->kg[2] * this->kg[2]*this->kg[3] -
                                                       8*this->mn2*(2*p0*this->q[0]*(this->q[1] + this->kg[1]) + this->q[0] * this->q[0]*(2*this->q[1] + this->kg[1]) -
                                                                    this->kg[1]*(this->q[1] * this->q[1] + this->q[3] * this->q[3] + 2*this->q[1]*this->kg[1] + 2*this->q[3]*this->kg[3])) +
                                                       2*p0*this->q[0]*(-pow(this->q[1],3) - 3*this->q[1] * this->q[1]*this->kg[1] + 2*this->q[0] * this->q[0]*(this->q[1] + 2*this->kg[1]) -
                                                                        this->q[1]*(this->q[3] * this->q[3] + this->kg[1] * this->kg[1] - 2*this->kg[2] * this->kg[2] + 3*this->q[3]*this->kg[3]) - this->kg[1]*(this->kg[1] * this->kg[1] + this->kg[2] * this->kg[2] + this->kg[3]*(this->q[3] + this->kg[3]))))))
               )/(24.*mDelta2*pow(this->mn,4));


    tr[3][0] = -(C3v*(this->c_i*C5aNC*this->q[1]*this->kg[2]*(8*pow(mDelta,3)*this->mn*(p0 + this->q[0]) -
                                                              4*this->mn2*(p0 + this->q[0])*(2*p0*this->q[0] + this->q[0] * this->q[0] - this->q[1]*this->kg[1] - this->q[3]*this->kg[3]) -
                                                              2*mDelta*this->mn*(p0 + this->q[0])*(8*this->mn2 + 2*p0*this->q[0] + this->q[0] * this->q[0] - this->q[1]*this->kg[1] - this->q[3]*this->kg[3]) +
                                                              mDelta2*(8*pow(p0,3) + 12*pow(p0,2)*this->q[0] + 4*p0*this->q[0] * this->q[0] + 8*this->mn2*(-2*p0 + this->q[0]) +
                                                                       this->q[0]*(this->q[1] * this->q[1] - this->q[1]*this->kg[1] + this->q[3]*(this->q[3] - this->kg[3])) - 2*p0*(this->q[1] * this->q[1] + 2*this->q[1]*this->kg[1] + this->q[3]*(this->q[3] + 2*this->kg[3])))) +
                      2*C3vNC*this->mn*(16*mDelta2*this->mn*this->q[0]*(p0 + this->q[0])*this->q[3] +
                                        4*pow(mDelta,3)*(4*p0*this->q[0]*this->q[3] - this->q[3]*(this->kg[1] * this->kg[1] + this->kg[2] * this->kg[2]) + this->q[1]*this->kg[1]*this->kg[3]) +
                                        4*this->mn*(p0 + this->q[0])*(-4*pow(p0,2)*this->q[0]*this->q[3] - 8*p0*this->q[0] * this->q[0]*this->q[3] - 3*pow(this->q[0],3)*this->q[3] +
                                                                      2*p0*this->q[3]*(this->kg[1] * this->kg[1] + this->kg[2] * this->kg[2]) + this->q[0]*this->q[3]*(this->q[1] * this->q[1] + this->q[3] * this->q[3] + this->q[1]*this->kg[1] + this->kg[1] * this->kg[1] + this->kg[2] * this->kg[2]) -
                                                                      2*p0*this->q[1]*this->kg[1]*this->kg[3] + this->q[0]*(this->q[1] * this->q[1] + 2*this->q[3] * this->q[3] - this->q[1]*this->kg[1])*this->kg[3]) +
                                        mDelta*(-16*pow(p0,3)*this->q[0]*this->q[3] - 4*pow(p0,2)*(8*this->q[0] * this->q[0]*this->q[3] - 3*this->q[3]*(this->kg[1] * this->kg[1] + this->kg[2] * this->kg[2]) + 3*this->q[1]*this->kg[1]*this->kg[3]) -
                                                (this->q[3]*(this->kg[1] * this->kg[1] + this->kg[2] * this->kg[2]) - this->q[1]*this->kg[1]*this->kg[3])*(-3*this->q[0] * this->q[0] + this->q[1] * this->q[1] + 2*this->q[1]*this->kg[1] + this->q[3]*(this->q[3] + 2*this->kg[3])) +
                                                4*p0*this->q[0]*(-3*this->q[0] * this->q[0]*this->q[3] + this->q[1]*this->kg[1]*(this->q[3] - 3*this->kg[3]) + this->q[1] * this->q[1]*(this->q[3] + this->kg[3]) +
                                                                 this->q[3]*(3*(this->kg[1] * this->kg[1] + this->kg[2] * this->kg[2]) + this->q[3]*(this->q[3] + 2*this->kg[3])))))))/(24.*mDelta2*pow(this->mn,4));

    tr[3][1] = (this->c_i*C3v*C5aNC*this->kg[2]*(-2*this->mn2*this->q[1]*(this->q[1] + this->kg[1])*(-(this->q[0]*(2*p0 + this->q[0])) + this->q[1]*this->kg[1] + this->q[3]*this->kg[3]) +
                                                 4*pow(mDelta,3)*this->mn*(8*this->mn2 - 2*this->q[0]*(p0 + this->q[0]) + this->q[3]*(this->q[3] + this->kg[3])) -
                                                 mDelta*this->mn*(-4*pow(p0,2)*this->q[0] * this->q[0] - 2*pow(this->q[0],4) + this->q[0] * this->q[0]*this->q[3] * this->q[3] - this->q[1] * this->q[1]*this->q[3] * this->q[3] - pow(this->q[3],4) +
                                                                  2*this->q[0] * this->q[0]*this->q[1]*this->kg[1] - 2*this->q[1]*this->q[3] * this->q[3]*this->kg[1] + 2*this->q[3] * this->q[3]*this->kg[1] * this->kg[1] + 2*this->q[3] * this->q[3]*this->kg[2] * this->kg[2] -
                                                                  this->q[3]*(-5*this->q[0] * this->q[0] + this->q[1] * this->q[1] + 3*this->q[3] * this->q[3] + 2*this->q[1]*this->kg[1])*this->kg[3] +
                                                                  2*p0*this->q[0]*(-3*this->q[0] * this->q[0] + 2*this->q[3] * this->q[3] + this->q[1]*this->kg[1] + 3*this->q[3]*this->kg[3]) +
                                                                  8*this->mn2*(2*p0*this->q[0] + this->q[0] * this->q[0] - this->q[1]*(this->q[1] + 2*this->kg[1]) - this->q[3]*(this->q[3] + 2*this->kg[3]))) +
                                                 mDelta2*(-8*pow(p0,3)*this->q[0] + 8*this->mn2*(2*p0*(2*p0 + this->q[0]) - this->q[1]*this->kg[1]) +
                                                          4*pow(p0,2)*(-3*this->q[0] * this->q[0] + this->q[3]*(this->q[3] + this->kg[3])) + 2*p0*this->q[0]*(-2*this->q[0] * this->q[0] + this->q[1]*(this->q[1] + 2*this->kg[1]) + this->q[3]*(this->q[3] + this->kg[3])) +
                                                          this->q[1]*(this->q[0] * this->q[0]*(-this->q[1] + this->kg[1]) + this->q[3]*(-(this->q[3]*this->kg[1]) + this->q[1]*this->kg[3])))) -
                2*C3v*C3vNC*this->mn*(16*mDelta2*this->mn*(p0 + this->q[0])*this->q[3]*this->kg[1] + 4*pow(mDelta,3)*(this->q[0]*this->q[1]*this->q[3] + 4*p0*this->q[3]*this->kg[1] + this->q[0]*this->kg[1]*(-2*this->q[3] + this->kg[3])) -
                                      4*this->mn*(p0 + this->q[0])*(2*p0*this->q[0]*(this->q[1]*this->q[3] + this->kg[1]*(2*this->q[3] + this->kg[3])) + this->q[0] * this->q[0]*(this->q[1]*(this->q[3] - this->kg[3]) + this->kg[1]*(2*this->q[3] + this->kg[3])) -
                                                                    this->q[3]*(this->q[1] * this->q[1]*this->kg[1] + this->q[1]*(2*this->kg[1] * this->kg[1] + this->kg[2] * this->kg[2]) + this->q[3]*this->kg[1]*(this->q[3] + 2*this->kg[3]))) +
                                      mDelta*(-4*pow(p0,2)*this->q[0]*(3*this->q[1]*this->q[3] + 2*this->q[3]*this->kg[1] + 3*this->kg[1]*this->kg[3]) +
                                              this->q[0]*(this->q[1]*this->q[3] + this->kg[1]*(-2*this->q[3] + this->kg[3]))*(-3*this->q[0] * this->q[0] + this->q[1] * this->q[1] + 2*this->q[1]*this->kg[1] + this->q[3]*(this->q[3] + 2*this->kg[3])) +
                                              4*p0*(this->q[0] * this->q[0]*(-3*this->q[1]*this->q[3] + 2*this->q[3]*this->kg[1] + this->q[1]*this->kg[3] - 3*this->kg[1]*this->kg[3]) +
                                                    this->q[3]*(this->q[1] * this->q[1]*this->kg[1] + 2*this->q[1]*(this->kg[1] * this->kg[1] + this->kg[2] * this->kg[2]) + this->q[3]*this->kg[1]*(this->q[3] + 2*this->kg[3]))))))/(24.*mDelta2*pow(this->mn,4));

    tr[3][2] = (C3v*(-2*C3vNC*this->mn*this->kg[2]*(16*mDelta2*this->mn*(p0 + this->q[0])*this->q[3] -
                                                    4*this->mn*(p0 + this->q[0])*(-(this->q[3]*(-2*this->q[0]*(2*p0 + this->q[0]) + this->q[3] * this->q[3] + this->q[1]*(this->q[1] + this->kg[1]))) +
                                                                                  (2*p0*this->q[0] + this->q[0] * this->q[0] - this->q[1] * this->q[1] - 2*this->q[3] * this->q[3])*this->kg[3]) + 4*pow(mDelta,3)*(4*p0*this->q[3] + this->q[0]*(-2*this->q[3] + this->kg[3])) +
                                                    mDelta*(-4*pow(p0,2)*this->q[0]*(2*this->q[3] + 3*this->kg[3]) -
                                                            this->q[0]*(2*this->q[3] - this->kg[3])*(-3*this->q[0] * this->q[0] + this->q[1] * this->q[1] + this->q[3] * this->q[3] + 2*this->q[1]*this->kg[1] + 2*this->q[3]*this->kg[3]) +
                                                            4*p0*(this->q[0] * this->q[0]*(2*this->q[3] - 3*this->kg[3]) + (this->q[1] * this->q[1] + this->q[3] * this->q[3])*(this->q[3] + 2*this->kg[3])))) -
                     this->c_i*C5aNC*(2*this->mn2*this->q[1]*this->kg[2] * this->kg[2]*(-(this->q[0]*(2*p0 + this->q[0])) + this->q[1]*this->kg[1] + this->q[3]*this->kg[3]) +
                                      4*pow(mDelta,3)*this->mn*(8*this->mn2*this->kg[1] - 2*p0*this->q[0]*this->kg[1] - 2*this->q[0] * this->q[0]*this->kg[1] + this->q[1] * this->q[1]*this->kg[1] + this->q[3] * this->q[3]*this->kg[1] +
                                                                this->q[1]*this->kg[1] * this->kg[1] + this->q[1]*this->kg[2] * this->kg[2] + this->q[3]*this->kg[1]*this->kg[3]) +
                                      mDelta2*(-8*pow(p0,3)*this->q[0]*this->kg[1] + this->q[1]*(-this->q[0] * this->q[0] + this->q[1] * this->q[1] + this->q[3] * this->q[3])*this->kg[2] * this->kg[2] +
                                               8*this->mn2*(4*pow(p0,2)*this->kg[1] + 2*p0*this->q[0]*(-this->q[1] + this->kg[1]) + this->q[1]*this->kg[2] * this->kg[2]) +
                                               4*pow(p0,2)*(this->q[0] * this->q[0]*(this->q[1] - 3*this->kg[1]) + this->q[1] * this->q[1]*this->kg[1] + this->q[1]*(this->kg[1] * this->kg[1] + this->kg[2] * this->kg[2]) + this->q[3]*this->kg[1]*(this->q[3] + this->kg[3])) +
                                               2*p0*this->q[0]*(-pow(this->q[1],3) + 2*this->q[0] * this->q[0]*(this->q[1] - this->kg[1]) + this->q[3]*this->kg[1]*(this->q[3] + this->kg[3]) -
                                                                this->q[1]*(-this->kg[1] * this->kg[1] + this->kg[2] * this->kg[2] + this->q[3]*(this->q[3] + this->kg[3])))) +
                                      mDelta*this->mn*(4*pow(this->q[0],4)*this->q[1] - 2*this->q[0] * this->q[0]*pow(this->q[1],3) - 2*this->q[0] * this->q[0]*this->q[1]*this->q[3] * this->q[3] + 2*pow(this->q[0],4)*this->kg[1] -
                                                       5*this->q[0] * this->q[0]*this->q[1] * this->q[1]*this->kg[1] + pow(this->q[1],4)*this->kg[1] - this->q[0] * this->q[0]*this->q[3] * this->q[3]*this->kg[1] + 2*this->q[1] * this->q[1]*this->q[3] * this->q[3]*this->kg[1] +
                                                       pow(this->q[3],4)*this->kg[1] - 5*this->q[0] * this->q[0]*this->q[1]*this->kg[1] * this->kg[1] + 3*pow(this->q[1],3)*this->kg[1] * this->kg[1] + 3*this->q[1]*this->q[3] * this->q[3]*this->kg[1] * this->kg[1] +
                                                       2*this->q[1] * this->q[1]*pow(this->kg[1],3) - 2*this->q[3] * this->q[3]*pow(this->kg[1],3) + 4*pow(p0,2)*this->q[0] * this->q[0]*(this->q[1] + this->kg[1]) -
                                                       3*this->q[0] * this->q[0]*this->q[1]*this->kg[2] * this->kg[2] + pow(this->q[1],3)*this->kg[2] * this->kg[2] + this->q[1]*this->q[3] * this->q[3]*this->kg[2] * this->kg[2] +
                                                       2*this->q[1] * this->q[1]*this->kg[1]*this->kg[2] * this->kg[2] - 2*this->q[3] * this->q[3]*this->kg[1]*this->kg[2] * this->kg[2] - 2*this->q[0] * this->q[0]*this->q[1]*this->q[3]*this->kg[3] - 5*this->q[0] * this->q[0]*this->q[3]*this->kg[1]*this->kg[3] +
                                                       3*this->q[1] * this->q[1]*this->q[3]*this->kg[1]*this->kg[3] + 3*pow(this->q[3],3)*this->kg[1]*this->kg[3] + 4*this->q[1]*this->q[3]*this->kg[1] * this->kg[1]*this->kg[3] + 2*this->q[1]*this->q[3]*this->kg[2] * this->kg[2]*this->kg[3] -
                                                       8*this->mn2*(2*p0*this->q[0]*(this->q[1] + this->kg[1]) + this->q[0] * this->q[0]*(2*this->q[1] + this->kg[1]) -
                                                                    this->kg[1]*(this->q[1] * this->q[1] + this->q[3] * this->q[3] + 2*this->q[1]*this->kg[1] + 2*this->q[3]*this->kg[3])) +
                                                       2*p0*this->q[0]*(-pow(this->q[1],3) - 3*this->q[1] * this->q[1]*this->kg[1] + this->q[0] * this->q[0]*(4*this->q[1] + 3*this->kg[1]) - this->q[3]*this->kg[1]*(2*this->q[3] + 3*this->kg[3]) -
                                                                        this->q[1]*(3*this->kg[1] * this->kg[1] + 2*this->kg[2] * this->kg[2] + this->q[3]*(this->q[3] + this->kg[3])))))))/(24.*mDelta2*pow(this->mn,4));

    tr[3][3] = -(C3v*(this->c_i*C5aNC*this->q[1]*this->kg[2]*(4*pow(mDelta,3)*this->mn*(this->q[3] + this->kg[3]) +
                                                              2*this->mn2*(this->q[3] + this->kg[3])*(-(this->q[0]*(2*p0 + this->q[0])) + this->q[1]*this->kg[1] + this->q[3]*this->kg[3]) +
                                                              mDelta*this->mn*(this->q[3] + this->kg[3])*(-4*p0*this->q[0] - 3*this->q[0] * this->q[0] + this->q[1] * this->q[1] + this->q[3] * this->q[3] + 2*this->q[1]*this->kg[1] + 2*this->q[3]*this->kg[3]) +
                                                              mDelta2*(-(this->q[1]*this->q[3]*this->kg[1]) + this->q[0] * this->q[0]*(this->q[3] - this->kg[3]) + 8*this->mn2*this->kg[3] - 2*p0*this->q[0]*this->kg[3] + this->q[1] * this->q[1]*this->kg[3] +
                                                                       4*pow(p0,2)*(this->q[3] + this->kg[3]))) + 2*C3vNC*this->mn*
                                                                                                                  (16*mDelta2*this->mn*(p0 + this->q[0])*(this->q[0] * this->q[0] - this->q[1]*this->kg[1]) +
                                                                                                                   4*pow(mDelta,3)*(4*p0*(this->q[0] * this->q[0] - this->q[1]*this->kg[1]) - this->q[0]*(pow(this->q[1] - this->kg[1],2) + this->kg[2] * this->kg[2])) -
                                                                                                                   4*this->mn*(p0 + this->q[0])*(4*pow(p0,2)*this->q[0] * this->q[0] + 3*pow(this->q[0],4) +
                                                                                                                                                 2*p0*this->q[0]*(4*this->q[0] * this->q[0] - pow(this->q[1] + this->kg[1],2) - this->kg[2] * this->kg[2]) -
                                                                                                                                                 this->q[0] * this->q[0]*(this->q[3] * this->q[3] + (this->q[1] + this->kg[1])*(2*this->q[1] + this->kg[1]) + this->kg[2] * this->kg[2] + 2*this->q[3]*this->kg[3]) +
                                                                                                                                                 this->q[1]*(this->q[1] * this->q[1]*this->kg[1] + this->q[1]*(2*this->kg[1] * this->kg[1] + this->kg[2] * this->kg[2]) + this->q[3]*this->kg[1]*(this->q[3] + 2*this->kg[3]))) -
                                                                                                                   mDelta*(16*pow(p0,3)*this->q[0] * this->q[0] + 4*pow(p0,2)*this->q[0]*
                                                                                                                                                                  (8*this->q[0] * this->q[0] - 3*this->q[1] * this->q[1] - 2*this->q[1]*this->kg[1] - 3*(this->kg[1] * this->kg[1] + this->kg[2] * this->kg[2])) +
                                                                                                                           this->q[0]*(pow(this->q[1] - this->kg[1],2) + this->kg[2] * this->kg[2])*(-3*this->q[0] * this->q[0] + this->q[1] * this->q[1] + this->q[3] * this->q[3] + 2*this->q[1]*this->kg[1] + 2*this->q[3]*this->kg[3]) +
                                                                                                                           4*p0*(3*pow(this->q[0],4) - this->q[0] * this->q[0]*(4*this->q[1] * this->q[1] - this->q[1]*this->kg[1] + 3*(this->kg[1] * this->kg[1] + this->kg[2] * this->kg[2]) + this->q[3]*(this->q[3] + 2*this->kg[3])) +
                                                                                                                                 this->q[1]*(this->q[1] * this->q[1]*this->kg[1] + 2*this->q[1]*(this->kg[1] * this->kg[1] + this->kg[2] * this->kg[2]) + this->q[3]*this->kg[1]*(this->q[3] + 2*this->kg[3])))))))/
               (24.*mDelta2*pow(this->mn,4));

}

void Hadronic_Current_R_Delta::set_tr_crs(Array4x4 &tr) {

    tr[0][0] = -(C3v*C3vNC*((mDelta + this->mn)*p0 - this->mn*this->q[0])*(4*mDelta2*this->q[1]*this->kg[1] - 4*pow(p0,2)*this->q[1]*this->kg[1] + pow(this->q[1],3)*this->kg[1] + this->q[1]*this->q[3] * this->q[3]*this->kg[1] +
                                                                           this->q[1] * this->q[1]*this->kg[1] * this->kg[1] - this->q[3] * this->q[3]*this->kg[1] * this->kg[1] - this->q[3] * this->q[3]*this->kg[2] * this->kg[2] +
                                                                           this->q[0] * this->q[0]*(this->q[1] * this->q[1] - 3*this->q[1]*this->kg[1] + this->q[3]*(2*this->q[3] - 3*this->kg[3])) +
                                                                           this->q[3]*(4*mDelta2 - 4*pow(p0,2) + this->q[1] * this->q[1] + this->q[3] * this->q[3] + 2*this->q[1]*this->kg[1])*this->kg[3] + 8*p0*this->q[0]*(this->q[1]*this->kg[1] + this->q[3]*this->kg[3])))/
               (3.*mDelta2*pow(this->mn,3));

    tr[0][1] = -(C3v*(this->c_i*C5aNC*mDelta*((mDelta + this->mn)*p0 - this->mn*this->q[0])*this->q[3]*this->kg[2]*
                      (8*this->mn2 + 2*(p0 - this->q[0])*this->q[0] + this->q[1]*(this->q[1] + this->kg[1]) + this->q[3]*(this->q[3] + this->kg[3])) +
                      C3vNC*this->mn*(16*mDelta2*this->mn*(p0 - this->q[0])*this->q[0]*this->kg[1] + 4*pow(mDelta,3)*(4*p0*this->q[0]*this->kg[1] + this->q[3]*(this->q[3]*this->kg[1] - this->q[1]*this->kg[3])) +
                                      4*this->mn*(p0 - this->q[0])*(pow(this->q[0],3)*(this->q[1] - 3*this->kg[1]) + 8*p0*this->q[0] * this->q[0]*this->kg[1] + 2*p0*this->q[3]*(-(this->q[3]*this->kg[1]) + this->q[1]*this->kg[3]) +
                                                                    this->q[0]*(-4*pow(p0,2)*this->kg[1] + 2*this->q[3] * this->q[3]*this->kg[1] + this->q[1]*this->kg[1]*(this->q[1] + this->kg[1]) + this->q[3]*(-this->q[1] + this->kg[1])*this->kg[3])) +
                                      mDelta*(-16*pow(p0,3)*this->q[0]*this->kg[1] + this->q[3]*(this->q[3]*this->kg[1] - this->q[1]*this->kg[3])*(-3*this->q[0] * this->q[0] + this->q[1] * this->q[1] + 2*this->q[1]*this->kg[1] + this->q[3]*(this->q[3] + 2*this->kg[3])) +
                                              4*pow(p0,2)*(8*this->q[0] * this->q[0]*this->kg[1] + 3*this->q[3]*(-(this->q[3]*this->kg[1]) + this->q[1]*this->kg[3])) +
                                              4*p0*this->q[0]*(this->q[0] * this->q[0]*(this->q[1] - 3*this->kg[1]) + this->q[1] * this->q[1]*this->kg[1] + this->q[3]*this->kg[1]*(4*this->q[3] + this->kg[3]) + this->q[1]*(this->kg[1] * this->kg[1] - 3*this->q[3]*this->kg[3]))))))/
               (12.*mDelta2*pow(this->mn,4));

    tr[0][2] = (this->c_i*C3v*C5aNC*mDelta*((mDelta + this->mn)*p0 - this->mn*this->q[0])*(this->q[3]*this->kg[1] - this->q[1]*this->kg[3])*
                (8*this->mn2 + 2*(p0 - this->q[0])*this->q[0] + this->q[1]*(this->q[1] + this->kg[1]) + this->q[3]*(this->q[3] + this->kg[3])) -
                C3v*C3vNC*this->mn*this->kg[2]*(16*mDelta2*this->mn*(p0 - this->q[0])*this->q[0] + 4*pow(mDelta,3)*(4*p0*this->q[0] + this->q[1] * this->q[1] + this->q[3] * this->q[3]) +
                                                4*this->mn*(p0 - this->q[0])*(-4*pow(p0,2)*this->q[0] + 8*p0*this->q[0] * this->q[0] - 3*pow(this->q[0],3) - 2*p0*(this->q[1] * this->q[1] + this->q[3] * this->q[3]) +
                                                                              this->q[0]*this->q[1]*(2*this->q[1] + this->kg[1]) + this->q[0]*this->q[3]*(2*this->q[3] + this->kg[3])) +
                                                mDelta*(-16*pow(p0,3)*this->q[0] + 4*pow(p0,2)*(8*this->q[0] * this->q[0] - 3*(this->q[1] * this->q[1] + this->q[3] * this->q[3])) +
                                                        (this->q[1] * this->q[1] + this->q[3] * this->q[3])*(-3*this->q[0] * this->q[0] + this->q[1] * this->q[1] + this->q[3] * this->q[3] + 2*this->q[1]*this->kg[1] + 2*this->q[3]*this->kg[3]) +
                                                        4*p0*this->q[0]*(-3*this->q[0] * this->q[0] + this->q[1]*(4*this->q[1] + this->kg[1]) + this->q[3]*(4*this->q[3] + this->kg[3])))))/(12.*mDelta2*pow(this->mn,4));

    tr[0][3] = (this->c_i*C3v*C5aNC*mDelta*((mDelta + this->mn)*p0 - this->mn*this->q[0])*this->q[1]*this->kg[2]*(8*this->mn2 + 2*(p0 - this->q[0])*this->q[0] + this->q[1]*(this->q[1] + this->kg[1]) + this->q[3]*(this->q[3] + this->kg[3])) -
                C3v*C3vNC*this->mn*(16*mDelta2*this->mn*(p0 - this->q[0])*this->q[0]*this->kg[3] + 4*pow(mDelta,3)*(-(this->q[1]*this->q[3]*this->kg[1]) + 4*p0*this->q[0]*this->kg[3] + this->q[1] * this->q[1]*this->kg[3]) +
                                    4*this->mn*(p0 - this->q[0])*(-(this->q[0]*this->q[3]*(this->kg[1]*(this->q[1] + this->kg[1]) + this->kg[2] * this->kg[2])) + pow(this->q[0],3)*(2*this->q[3] - 3*this->kg[3]) + 8*p0*this->q[0] * this->q[0]*this->kg[3] +
                                                                  this->q[0]*(-4*pow(p0,2) + this->q[3] * this->q[3] + this->q[1]*(2*this->q[1] + this->kg[1]))*this->kg[3] + 2*p0*this->q[1]*(this->q[3]*this->kg[1] - this->q[1]*this->kg[3])) +
                                    mDelta*(-16*pow(p0,3)*this->q[0]*this->kg[3] + 4*pow(p0,2)*(3*this->q[1]*this->q[3]*this->kg[1] + 8*this->q[0] * this->q[0]*this->kg[3] - 3*this->q[1] * this->q[1]*this->kg[3]) +
                                            this->q[1]*(-(this->q[3]*this->kg[1]) + this->q[1]*this->kg[3])*(-3*this->q[0] * this->q[0] + this->q[1] * this->q[1] + this->q[3] * this->q[3] + 2*this->q[1]*this->kg[1] + 2*this->q[3]*this->kg[3]) +
                                            4*p0*this->q[0]*(-(this->q[3]*(3*this->q[1]*this->kg[1] + this->kg[1] * this->kg[1] + this->kg[2] * this->kg[2])) + this->q[0] * this->q[0]*(2*this->q[3] - 3*this->kg[3]) +
                                                             (this->q[3] * this->q[3] + this->q[1]*(4*this->q[1] + this->kg[1]))*this->kg[3]))))/(12.*mDelta2*pow(this->mn,4));


    tr[1][0] = (this->c_i*C3v*C5aNC*this->q[3]*this->kg[2]*(8*pow(mDelta,3)*this->mn*(-p0 + this->q[0]) +
                                                            2*mDelta*this->mn*(p0 - this->q[0])*(8*this->mn2 - 2*p0*this->q[0] + this->q[0] * this->q[0] - this->q[1]*this->kg[1] - this->q[3]*this->kg[3]) -
                                                            4*this->mn2*(p0 - this->q[0])*(2*p0*this->q[0] - this->q[0] * this->q[0] + this->q[1]*this->kg[1] + this->q[3]*this->kg[3]) +
                                                            mDelta2*(-8*pow(p0,3) + 12*pow(p0,2)*this->q[0] + 8*this->mn2*(2*p0 + this->q[0]) +
                                                                     this->q[0]*(this->q[1] * this->q[1] + this->q[3] * this->q[3] - this->q[1]*this->kg[1] - this->q[3]*this->kg[3]) +
                                                                     2*p0*(-2*this->q[0] * this->q[0] + this->q[1] * this->q[1] + this->q[3] * this->q[3] + 2*this->q[1]*this->kg[1] + 2*this->q[3]*this->kg[3]))) -
                2*C3v*C3vNC*this->mn*(16*mDelta2*this->mn*(p0 - this->q[0])*this->q[0]*this->q[1] +
                                      4*pow(mDelta,3)*(4*p0*this->q[0]*this->q[1] + this->q[0] * this->q[0]*this->q[1] - this->kg[1]*(this->q[1]*this->kg[1] + this->q[3]*this->kg[3])) +
                                      mDelta*(-16*pow(p0,3)*this->q[0]*this->q[1] + 4*p0*this->q[0]*(pow(this->q[1],3) + 2*this->q[1] * this->q[1]*this->kg[1] - 3*this->q[1]*this->kg[1] * this->kg[1] + this->q[3]*this->kg[1]*(this->q[3] - 3*this->kg[3]) +
                                                                                                     this->q[1]*this->q[3]*(this->q[3] + this->kg[3])) - (3*this->q[0] * this->q[0] - this->q[1]*(this->q[1] + 2*this->kg[1]) - this->q[3]*(this->q[3] + 2*this->kg[3]))*(this->q[0] * this->q[0]*this->q[1] - this->kg[1]*(this->q[1]*this->kg[1] + this->q[3]*this->kg[3])) +
                                              4*pow(p0,2)*(5*this->q[0] * this->q[0]*this->q[1] + 3*this->kg[1]*(this->q[1]*this->kg[1] + this->q[3]*this->kg[3]))) -
                                      4*this->mn*(p0 - this->q[0])*(4*pow(p0,2)*this->q[0]*this->q[1] - 2*p0*(3*this->q[0] * this->q[0]*this->q[1] + this->kg[1]*(this->q[1]*this->kg[1] + this->q[3]*this->kg[3])) +
                                                                    this->q[0]*(2*this->q[0] * this->q[0]*this->q[1] - pow(this->q[1],3) - 2*this->q[1] * this->q[1]*this->kg[1] + this->q[3]*this->kg[1]*(-this->q[3] + this->kg[3]) + this->q[1]*(this->kg[1] * this->kg[1] - this->q[3]*(this->q[3] + this->kg[3]))))))/
               (24.*mDelta2*pow(this->mn,4));

    tr[1][1] = (this->c_i*C3v*C5aNC*this->q[3]*this->kg[2]*(4*pow(mDelta,3)*this->mn*(this->q[1] + this->kg[1]) +
                                                            2*this->mn2*(this->q[1] + this->kg[1])*(2*p0*this->q[0] - this->q[0] * this->q[0] + this->q[1]*this->kg[1] + this->q[3]*this->kg[3]) +
                                                            mDelta*this->mn*(this->q[1] + this->kg[1])*(4*p0*this->q[0] - 3*this->q[0] * this->q[0] + this->q[1] * this->q[1] + this->q[3] * this->q[3] + 2*this->q[1]*this->kg[1] + 2*this->q[3]*this->kg[3]) +
                                                            mDelta2*(this->q[0] * this->q[0]*(this->q[1] - this->kg[1]) + 8*this->mn2*this->kg[1] + 2*p0*this->q[0]*this->kg[1] + this->q[3] * this->q[3]*this->kg[1] + 4*pow(p0,2)*(this->q[1] + this->kg[1]) -
                                                                     this->q[1]*this->q[3]*this->kg[3])) - 2*C3v*C3vNC*this->mn*(16*mDelta2*this->mn*(p0 - this->q[0])*(this->q[0] * this->q[0] - this->q[3]*this->kg[3]) +
                                                                                                                                 4*pow(mDelta,3)*(this->q[0]*(4*p0*this->q[0] + this->q[0] * this->q[0] + this->q[3] * this->q[3] - this->kg[1] * this->kg[1]) - 2*(2*p0 + this->q[0])*this->q[3]*this->kg[3]) -
                                                                                                                                 4*this->mn*(p0 - this->q[0])*(4*pow(p0,2)*this->q[0] * this->q[0] + 2*pow(this->q[0],4) +
                                                                                                                                                               2*p0*this->q[0]*(-3*this->q[0] * this->q[0] + this->q[3] * this->q[3] - this->kg[1] * this->kg[1] + 2*this->q[3]*this->kg[3]) -
                                                                                                                                                               this->q[0] * this->q[0]*(this->q[1] * this->q[1] + 2*this->q[1]*this->kg[1] - this->kg[1] * this->kg[1] + 3*this->q[3]*this->kg[3]) +
                                                                                                                                                               this->q[3]*(-(this->q[3]*(2*this->kg[1] * this->kg[1] + this->kg[2] * this->kg[2])) + (this->q[1] * this->q[1] + this->q[3] * this->q[3] + 2*this->q[1]*this->kg[1])*this->kg[3])) +
                                                                                                                                 mDelta*(-16*pow(p0,3)*this->q[0] * this->q[0] + 4*pow(p0,2)*this->q[0]*(5*this->q[0] * this->q[0] - 3*this->q[3] * this->q[3] + 3*this->kg[1] * this->kg[1] - 2*this->q[3]*this->kg[3]) -
                                                                                                                                         this->q[0]*(this->q[0] * this->q[0] + this->q[3] * this->q[3] - this->kg[1] * this->kg[1] - 2*this->q[3]*this->kg[3])*(3*this->q[0] * this->q[0] - this->q[1]*(this->q[1] + 2*this->kg[1]) - this->q[3]*(this->q[3] + 2*this->kg[3])) +
                                                                                                                                         4*p0*(this->q[0] * this->q[0]*((this->q[1] - this->kg[1])*(this->q[1] + 3*this->kg[1]) + this->q[3]*(2*this->q[3] - this->kg[3])) +
                                                                                                                                               this->q[3]*(2*this->q[3]*this->kg[1] * this->kg[1] - (this->q[1] * this->q[1] + this->q[3] * this->q[3] + 2*this->q[1]*this->kg[1])*this->kg[3])))))/(24.*mDelta2*pow(this->mn,4));

    tr[1][2] = (C3v*(-2*C3vNC*this->mn*this->kg[2]*(16*mDelta2*this->mn*(p0 - this->q[0])*this->q[1] + 4*pow(mDelta,3)*(4*p0*this->q[1] + 2*this->q[0]*this->q[1] - this->q[0]*this->kg[1]) +
                                                    4*this->mn*(p0 - this->q[0])*(pow(this->q[1],3) + this->q[1]*this->q[3] * this->q[3] + 2*this->q[1] * this->q[1]*this->kg[1] + this->q[3] * this->q[3]*this->kg[1] + 2*p0*this->q[0]*(2*this->q[1] + this->kg[1]) -
                                                                                  this->q[0] * this->q[0]*(2*this->q[1] + this->kg[1]) + this->q[1]*this->q[3]*this->kg[3]) +
                                                    mDelta*(4*pow(p0,2)*this->q[0]*(2*this->q[1] + 3*this->kg[1]) + 4*p0*
                                                                                                                    (this->q[0] * this->q[0]*(2*this->q[1] - 3*this->kg[1]) + (this->q[1] * this->q[1] + this->q[3] * this->q[3])*(this->q[1] + 2*this->kg[1])) +
                                                            this->q[0]*(2*this->q[1] - this->kg[1])*(-3*this->q[0] * this->q[0] + this->q[1] * this->q[1] + this->q[3] * this->q[3] + 2*this->q[1]*this->kg[1] + 2*this->q[3]*this->kg[3]))) +
                     this->c_i*C5aNC*(2*this->mn2*this->q[3]*this->kg[2] * this->kg[2]*(2*p0*this->q[0] - this->q[0] * this->q[0] + this->q[1]*this->kg[1] + this->q[3]*this->kg[3]) +
                                      4*pow(mDelta,3)*this->mn*(-(this->q[3]*this->kg[1] * this->kg[1]) + this->q[0] * this->q[0]*(this->q[3] - 2*this->kg[3]) + 2*p0*this->q[0]*this->kg[3] +
                                                                (8*this->mn2 + this->q[3] * this->q[3] + this->q[1]*(this->q[1] + this->kg[1]))*this->kg[3]) +
                                      mDelta2*(2*p0*this->q[3]*(this->q[0]*(4*p0*this->q[0] - 3*this->q[0] * this->q[0] + this->q[1] * this->q[1] + this->q[3] * this->q[3]) + this->q[0]*this->q[1]*this->kg[1] +
                                                                (-2*p0 + this->q[0])*this->kg[1] * this->kg[1]) + this->q[3]*(4*p0*this->q[0] - this->q[0] * this->q[0] + this->q[1] * this->q[1] + this->q[3] * this->q[3])*this->kg[2] * this->kg[2] +
                                               2*p0*(4*pow(p0,2)*this->q[0] + 2*pow(this->q[0],3) - this->q[0]*this->q[1]*(this->q[1] + this->kg[1]) + 2*p0*(-3*this->q[0] * this->q[0] + this->q[3] * this->q[3] + this->q[1]*(this->q[1] + this->kg[1])))*
                                               this->kg[3] + 8*this->mn2*(this->q[3]*this->kg[2] * this->kg[2] + 2*p0*this->q[0]*(this->q[3] - this->kg[3]) + 4*pow(p0,2)*this->kg[3])) +
                                      mDelta*this->mn*(this->q[3]*(4*pow(p0,2)*this->q[0] * this->q[0] - pow(this->q[0],4) - this->kg[1] * this->kg[1]*(3*(this->q[1] * this->q[1] + this->q[3] * this->q[3]) + 4*this->q[1]*this->kg[1]) -
                                                                   2*(this->q[3] * this->q[3] + this->q[1]*(this->q[1] + this->kg[1]))*this->kg[2] * this->kg[2] +
                                                                   2*p0*this->q[0]*(-this->q[0] * this->q[0] + this->q[1] * this->q[1] + this->q[3] * this->q[3] + this->q[1]*this->kg[1] - 3*this->kg[1] * this->kg[1] - this->kg[2] * this->kg[2]) +
                                                                   this->q[0] * this->q[0]*(this->q[1] * this->q[1] + this->q[3] * this->q[3] + 2*this->q[1]*this->kg[1] + 5*this->kg[1] * this->kg[1] + 2*this->kg[2] * this->kg[2])) +
                                                       (4*pow(p0,2)*this->q[0] * this->q[0] + 2*pow(this->q[0],4) + pow(this->q[1] * this->q[1] + this->q[3] * this->q[3],2) +
                                                        3*this->q[1]*(this->q[1] * this->q[1] + this->q[3] * this->q[3])*this->kg[1] + 2*(this->q[1] - this->q[3])*(this->q[1] + this->q[3])*this->kg[1] * this->kg[1] +
                                                        2*p0*this->q[0]*(-3*this->q[0] * this->q[0] + 2*this->q[1] * this->q[1] + 3*this->q[3] * this->q[3] + 3*this->q[1]*this->kg[1]) -
                                                        this->q[0] * this->q[0]*(3*(this->q[1] * this->q[1] + this->q[3] * this->q[3]) + 5*this->q[1]*this->kg[1]))*this->kg[3] +
                                                       8*this->mn2*(2*p0*this->q[0]*(this->q[3] + this->kg[3]) - this->q[0] * this->q[0]*(2*this->q[3] + this->kg[3]) + this->kg[3]*(this->q[1] * this->q[1] + this->q[3] * this->q[3] + 2*this->q[1]*this->kg[1] + 2*this->q[3]*this->kg[3]))
                                      ))))/(24.*mDelta2*pow(this->mn,4));

    tr[1][3] = (C3v*(2*C3vNC*this->mn*(16*mDelta2*this->mn*(-p0 + this->q[0])*this->q[1]*this->kg[3] + 4*pow(mDelta,3)*(this->q[0]*this->q[1]*(this->q[3] - 2*this->kg[3]) - 4*p0*this->q[1]*this->kg[3] + this->q[0]*this->kg[1]*this->kg[3]) -
                                       4*this->mn*(p0 - this->q[0])*(-(this->q[1]*this->q[3]*(2*this->kg[1] * this->kg[1] + this->kg[2] * this->kg[2])) + this->q[1]*(this->q[1] * this->q[1] + this->q[3] * this->q[3] + 2*this->q[1]*this->kg[1])*this->kg[3] +
                                                                     this->q[0] * this->q[0]*(this->q[3]*(this->q[1] + this->kg[1]) - (2*this->q[1] + this->kg[1])*this->kg[3]) + 2*p0*this->q[0]*(this->kg[1]*this->kg[3] + this->q[1]*(this->q[3] + 2*this->kg[3]))) +
                                       mDelta*(4*p0*this->q[3]*(this->q[0] * this->q[0]*(this->q[1] - this->kg[1]) + 2*this->q[1]*this->kg[1] * this->kg[1]) -
                                               4*p0*(this->q[0] * this->q[0]*(2*this->q[1] - 3*this->kg[1]) + this->q[1]*(this->q[1] * this->q[1] + this->q[3] * this->q[3] + 2*this->q[1]*this->kg[1]))*this->kg[3] -
                                               4*pow(p0,2)*this->q[0]*(3*this->q[1]*this->q[3] + 2*this->q[1]*this->kg[3] + 3*this->kg[1]*this->kg[3]) +
                                               this->q[0]*(this->q[1]*(this->q[3] - 2*this->kg[3]) + this->kg[1]*this->kg[3])*(-3*this->q[0] * this->q[0] + this->q[1] * this->q[1] + 2*this->q[1]*this->kg[1] + this->q[3]*(this->q[3] + 2*this->kg[3])))) -
                     this->c_i*C5aNC*this->kg[2]*(4*pow(mDelta,3)*this->mn*(8*this->mn2 + 2*(p0 - this->q[0])*this->q[0] + this->q[1]*(this->q[1] + this->kg[1])) -
                                                  2*this->mn2*this->q[3]*(this->q[3] + this->kg[3])*(2*p0*this->q[0] - this->q[0] * this->q[0] + this->q[1]*this->kg[1] + this->q[3]*this->kg[3]) +
                                                  mDelta*this->mn*(4*pow(p0,2)*this->q[0] * this->q[0] - 6*p0*pow(this->q[0],3) + 2*pow(this->q[0],4) + 4*p0*this->q[0]*this->q[1] * this->q[1] -
                                                                   3*this->q[0] * this->q[0]*this->q[1] * this->q[1] + pow(this->q[1],4) + this->q[1] * this->q[1]*this->q[3] * this->q[3] + 6*p0*this->q[0]*this->q[1]*this->kg[1] - 5*this->q[0] * this->q[0]*this->q[1]*this->kg[1] +
                                                                   3*pow(this->q[1],3)*this->kg[1] + this->q[1]*this->q[3] * this->q[3]*this->kg[1] + 2*this->q[1] * this->q[1]*this->kg[1] * this->kg[1] + 2*this->q[3]*((p0 - this->q[0])*this->q[0] + this->q[1]*(this->q[1] + this->kg[1]))*this->kg[3] +
                                                                   8*this->mn2*(2*p0*this->q[0] - this->q[0] * this->q[0] + this->q[1] * this->q[1] + this->q[3] * this->q[3] + 2*this->q[1]*this->kg[1] + 2*this->q[3]*this->kg[3])) +
                                                  mDelta2*(8*pow(p0,3)*this->q[0] + 4*pow(p0,2)*(-3*this->q[0] * this->q[0] + this->q[1]*(this->q[1] + this->kg[1])) +
                                                           8*this->mn2*(4*pow(p0,2) - 2*p0*this->q[0] - this->q[3]*this->kg[3]) + 2*p0*this->q[0]*(2*this->q[0] * this->q[0] - this->q[1]*(this->q[1] + this->kg[1]) - this->q[3]*(this->q[3] + 2*this->kg[3])) +
                                                           this->q[3]*(this->q[0] * this->q[0]*(-this->q[3] + this->kg[3]) + this->q[1]*(this->q[3]*this->kg[1] - this->q[1]*this->kg[3]))))))/(24.*mDelta2*pow(this->mn,4));


    tr[2][0] = (C3v*(2*C3vNC*this->mn*this->kg[2]*(4*pow(mDelta,3)*(this->q[1]*this->kg[1] + this->q[3]*this->kg[3]) -
                                                   4*this->mn*(p0 - this->q[0])*(this->q[0]*(this->q[1] * this->q[1] + this->q[3] * this->q[3] - this->q[1]*this->kg[1] - this->q[3]*this->kg[3]) + 2*p0*(this->q[1]*this->kg[1] + this->q[3]*this->kg[3])) +
                                                   mDelta*(-4*p0*this->q[0]*(this->q[1] * this->q[1] + this->q[3] * this->q[3] - 3*this->q[1]*this->kg[1] - 3*this->q[3]*this->kg[3]) - 12*pow(p0,2)*(this->q[1]*this->kg[1] + this->q[3]*this->kg[3]) +
                                                           (this->q[1]*this->kg[1] + this->q[3]*this->kg[3])*(-3*this->q[0] * this->q[0] + this->q[1] * this->q[1] + this->q[3] * this->q[3] + 2*this->q[1]*this->kg[1] + 2*this->q[3]*this->kg[3]))) -
                     this->c_i*C5aNC*(this->q[3]*this->kg[1] - this->q[1]*this->kg[3])*(8*pow(mDelta,3)*this->mn*(-p0 + this->q[0]) +
                                                                                        2*mDelta*this->mn*(p0 - this->q[0])*(8*this->mn2 - 2*p0*this->q[0] + this->q[0] * this->q[0] - this->q[1]*this->kg[1] - this->q[3]*this->kg[3]) -
                                                                                        4*this->mn2*(p0 - this->q[0])*(2*p0*this->q[0] - this->q[0] * this->q[0] + this->q[1]*this->kg[1] + this->q[3]*this->kg[3]) +
                                                                                        mDelta2*(-8*pow(p0,3) + 12*pow(p0,2)*this->q[0] + 8*this->mn2*(2*p0 + this->q[0]) +
                                                                                                 this->q[0]*(this->q[1] * this->q[1] + this->q[3] * this->q[3] - this->q[1]*this->kg[1] - this->q[3]*this->kg[3]) +
                                                                                                 2*p0*(-2*this->q[0] * this->q[0] + this->q[1] * this->q[1] + this->q[3] * this->q[3] + 2*this->q[1]*this->kg[1] + 2*this->q[3]*this->kg[3])))))/(24.*mDelta2*pow(this->mn,4));

    tr[2][1] = (C3v*(2*C3vNC*this->mn*this->kg[2]*(4*pow(mDelta,3)*this->q[0]*this->kg[1] - 4*this->mn*(p0 - this->q[0])*(this->q[0] * this->q[0]*(this->q[1] - this->kg[1]) + 2*p0*this->q[0]*this->kg[1] + this->q[3]*(this->q[3]*this->kg[1] - this->q[1]*this->kg[3])) +
                                                   mDelta*(-12*pow(p0,2)*this->q[0]*this->kg[1] + this->q[0]*this->kg[1]*(-3*this->q[0] * this->q[0] + this->q[1] * this->q[1] + this->q[3] * this->q[3] + 2*this->q[1]*this->kg[1] + 2*this->q[3]*this->kg[3]) -
                                                           4*p0*(this->q[0] * this->q[0]*(this->q[1] - 3*this->kg[1]) + 2*this->q[3]*(this->q[3]*this->kg[1] - this->q[1]*this->kg[3])))) -
                     this->c_i*C5aNC*(-2*this->mn2*(this->q[1] + this->kg[1])*(-(this->q[3]*this->kg[1]) + this->q[1]*this->kg[3])*(2*p0*this->q[0] - this->q[0] * this->q[0] + this->q[1]*this->kg[1] + this->q[3]*this->kg[3]) +
                                      4*pow(mDelta,3)*this->mn*(this->q[3]*(this->q[0] * this->q[0] + this->q[1]*this->kg[1] - this->kg[2] * this->kg[2]) + (8*this->mn2 + 2*(p0 - this->q[0])*this->q[0] + this->q[3] * this->q[3])*this->kg[3]) +
                                      mDelta2*(8*pow(p0,3)*this->q[0]*this->kg[3] + 4*pow(p0,2)*
                                                                                    (this->q[3]*this->kg[1]*(this->q[1] + this->kg[1]) + this->q[0] * this->q[0]*(this->q[3] - 3*this->kg[3]) + this->q[3]*this->kg[3]*(this->q[3] + this->kg[3])) +
                                               2*p0*this->q[0]*(this->q[3]*(this->q[1] * this->q[1] + this->q[3] * this->q[3] + this->q[1]*this->kg[1] + 2*this->kg[1] * this->kg[1] + this->kg[2] * this->kg[2]) - this->q[1]*(this->q[1] + 2*this->kg[1])*this->kg[3] +
                                                                this->q[0] * this->q[0]*(-3*this->q[3] + 2*this->kg[3])) + (this->q[3]*this->kg[1] - this->q[1]*this->kg[3])*(this->q[0] * this->q[0]*(this->q[1] - this->kg[1]) + this->q[3]*(this->q[3]*this->kg[1] - this->q[1]*this->kg[3])) +
                                               8*this->mn2*(2*p0*this->q[0]*(this->q[3] - this->kg[3]) + 4*pow(p0,2)*this->kg[3] + this->kg[1]*(this->q[3]*this->kg[1] - this->q[1]*this->kg[3]))) +
                                      mDelta*this->mn*(this->q[3]*(-pow(this->q[0],4) + this->kg[1]*(pow(this->q[1],3) + this->q[1]*this->q[3] * this->q[3] + 2*this->q[1] * this->q[1]*this->kg[1] - 2*this->q[3] * this->q[3]*this->kg[1]) -
                                                                   (this->q[1] * this->q[1] + 3*this->q[3] * this->q[3] + 2*this->q[1]*this->kg[1])*this->kg[2] * this->kg[2] +
                                                                   this->q[0] * this->q[0]*(-this->q[1] * this->q[1] + this->q[3] * this->q[3] - 3*this->q[1]*this->kg[1] + 2*this->kg[1] * this->kg[1] + 5*this->kg[2] * this->kg[2])) +
                                                       (2*pow(this->q[0],4) - this->q[0] * this->q[0]*(3*this->q[3] * this->q[3] + 2*this->q[1]*this->kg[1]) +
                                                        this->q[3] * this->q[3]*(this->q[1] * this->q[1] + this->q[3] * this->q[3] + 4*this->q[1]*this->kg[1] - 2*this->kg[2] * this->kg[2]))*this->kg[3] + 4*pow(p0,2)*this->q[0] * this->q[0]*(this->q[3] + this->kg[3]) +
                                                       8*this->mn2*(2*p0*this->q[0]*(this->q[3] + this->kg[3]) - this->q[0] * this->q[0]*(2*this->q[3] + this->kg[3]) +
                                                                    this->kg[3]*(this->q[1] * this->q[1] + this->q[3] * this->q[3] + 2*this->q[1]*this->kg[1] + 2*this->q[3]*this->kg[3])) +
                                                       2*p0*this->q[0]*(this->q[1] * this->q[1]*this->q[3] + this->q[1]*this->kg[1]*(3*this->q[3] + this->kg[3]) - this->q[0] * this->q[0]*(this->q[3] + 3*this->kg[3]) +
                                                                        this->q[3]*(this->q[3] * this->q[3] - this->kg[1] * this->kg[1] - 3*this->kg[2] * this->kg[2] + 3*this->q[3]*this->kg[3]))))))/(24.*mDelta2*pow(this->mn,4));

    tr[2][2] = -(C3v*(this->c_i*C5aNC*this->kg[2]*(this->q[3]*this->kg[1] - this->q[1]*this->kg[3])*(4*pow(mDelta,3)*this->mn +
                                                                                                     mDelta2*(8*this->mn2 + 4*pow(p0,2) + 2*p0*this->q[0] - this->q[0] * this->q[0] + this->q[1] * this->q[1] + this->q[3] * this->q[3]) +
                                                                                                     2*this->mn2*(2*p0*this->q[0] - this->q[0] * this->q[0] + this->q[1]*this->kg[1] + this->q[3]*this->kg[3]) +
                                                                                                     mDelta*this->mn*(4*p0*this->q[0] - 3*this->q[0] * this->q[0] + this->q[1] * this->q[1] + this->q[3] * this->q[3] + 2*this->q[1]*this->kg[1] + 2*this->q[3]*this->kg[3])) +
                      2*C3vNC*this->mn*(16*mDelta2*this->mn*(p0 - this->q[0])*(this->q[0] * this->q[0] - this->q[1]*this->kg[1] - this->q[3]*this->kg[3]) +
                                        4*pow(mDelta,3)*(this->q[0]*(this->q[0] * this->q[0] + this->q[1] * this->q[1] + this->q[3] * this->q[3] - 2*this->q[1]*this->kg[1] - this->kg[2] * this->kg[2] - 2*this->q[3]*this->kg[3]) +
                                                         4*p0*(this->q[0] * this->q[0] - this->q[1]*this->kg[1] - this->q[3]*this->kg[3])) -
                                        4*this->mn*(p0 - this->q[0])*(4*pow(p0,2)*this->q[0] * this->q[0] + 2*pow(this->q[0],4) + pow(this->q[1],3)*this->kg[1] + this->q[1]*this->q[3] * this->q[3]*this->kg[1] +
                                                                      this->q[1] * this->q[1]*this->kg[1] * this->kg[1] - this->q[3] * this->q[3]*this->kg[1] * this->kg[1] - this->q[1] * this->q[1]*this->kg[2] * this->kg[2] - 2*this->q[3] * this->q[3]*this->kg[2] * this->kg[2] +
                                                                      this->q[3]*(this->q[1] * this->q[1] + this->q[3] * this->q[3] + 2*this->q[1]*this->kg[1])*this->kg[3] +
                                                                      2*p0*this->q[0]*(-3*this->q[0] * this->q[0] + this->q[1] * this->q[1] + this->q[3] * this->q[3] + 2*this->q[1]*this->kg[1] - this->kg[2] * this->kg[2] + 2*this->q[3]*this->kg[3]) -
                                                                      this->q[0] * this->q[0]*(this->q[1] * this->q[1] + 3*this->q[1]*this->kg[1] - this->kg[2] * this->kg[2] + 3*this->q[3]*this->kg[3])) -
                                        mDelta*(16*pow(p0,3)*this->q[0] * this->q[0] + 4*pow(p0,2)*this->q[0]*
                                                                                       (-5*this->q[0] * this->q[0] + 3*this->q[1] * this->q[1] + 3*this->q[3] * this->q[3] + 2*this->q[1]*this->kg[1] - 3*this->kg[2] * this->kg[2] + 2*this->q[3]*this->kg[3]) +
                                                this->q[0]*(this->q[0] * this->q[0] + this->q[1] * this->q[1] + this->q[3] * this->q[3] - 2*this->q[1]*this->kg[1] - this->kg[2] * this->kg[2] - 2*this->q[3]*this->kg[3])*
                                                (3*this->q[0] * this->q[0] - this->q[1]*(this->q[1] + 2*this->kg[1]) - this->q[3]*(this->q[3] + 2*this->kg[3])) +
                                                4*p0*(pow(this->q[1],3)*this->kg[1] + this->q[1]*this->q[3] * this->q[3]*this->kg[1] + this->q[3] * this->q[3]*(-2*this->kg[2] * this->kg[2] + this->q[3]*this->kg[3]) +
                                                      this->q[0] * this->q[0]*(this->q[1]*(-4*this->q[1] + this->kg[1]) + 3*this->kg[2] * this->kg[2] + this->q[3]*(-2*this->q[3] + this->kg[3])) +
                                                      this->q[1] * this->q[1]*(2*this->kg[1] * this->kg[1] + this->kg[3]*(this->q[3] + 2*this->kg[3])))))))/(24.*mDelta2*pow(this->mn,4));

    tr[2][3] = (C3v*(2*C3vNC*this->mn*this->kg[2]*(4*pow(mDelta,3)*this->q[0]*this->kg[3] - 4*this->mn*(p0 - this->q[0])*
                                                                                            (this->q[0] * this->q[0]*(this->q[3] - this->kg[3]) + 2*p0*this->q[0]*this->kg[3] + this->q[1]*(-(this->q[3]*this->kg[1]) + this->q[1]*this->kg[3])) +
                                                   mDelta*(-12*pow(p0,2)*this->q[0]*this->kg[3] + this->q[0]*this->kg[3]*(-3*this->q[0] * this->q[0] + this->q[1] * this->q[1] + this->q[3] * this->q[3] + 2*this->q[1]*this->kg[1] + 2*this->q[3]*this->kg[3]) -
                                                           4*p0*(this->q[0] * this->q[0]*(this->q[3] - 3*this->kg[3]) + 2*this->q[1]*(-(this->q[3]*this->kg[1]) + this->q[1]*this->kg[3])))) +
                     this->c_i*C5aNC*(-2*this->mn2*(this->q[3] + this->kg[3])*(this->q[3]*this->kg[1] - this->q[1]*this->kg[3])*(2*p0*this->q[0] - this->q[0] * this->q[0] + this->q[1]*this->kg[1] + this->q[3]*this->kg[3]) +
                                      4*pow(mDelta,3)*this->mn*(this->q[0] * this->q[0]*(this->q[1] - 2*this->kg[1]) + 8*this->mn2*this->kg[1] + 2*p0*this->q[0]*this->kg[1] + this->q[1]*(this->q[1]*this->kg[1] - this->kg[2] * this->kg[2] + this->q[3]*this->kg[3])) +
                                      mDelta2*(8*pow(p0,3)*this->q[0]*this->kg[1] + 4*pow(p0,2)*
                                                                                    (this->q[0] * this->q[0]*(this->q[1] - 3*this->kg[1]) + this->q[1]*this->kg[1]*(this->q[1] + this->kg[1]) + this->q[1]*this->kg[3]*(this->q[3] + this->kg[3])) +
                                               (this->q[3]*this->kg[1] - this->q[1]*this->kg[3])*(this->q[0] * this->q[0]*(-this->q[3] + this->kg[3]) + this->q[1]*(this->q[3]*this->kg[1] - this->q[1]*this->kg[3])) +
                                               8*this->mn2*(2*p0*this->q[0]*(this->q[1] - this->kg[1]) + 4*pow(p0,2)*this->kg[1] + this->kg[3]*(-(this->q[3]*this->kg[1]) + this->q[1]*this->kg[3])) -
                                               2*p0*this->q[0]*(-pow(this->q[1],3) + this->q[0] * this->q[0]*(this->q[1] - 2*this->kg[1]) + this->q[3]*this->kg[1]*(this->q[3] + 2*this->kg[3]) +
                                                                this->q[1]*(2*this->kg[1] * this->kg[1] + this->kg[2] * this->kg[2] - this->q[3]*(this->q[3] + this->kg[3])))) +
                                      mDelta*this->mn*(pow(this->q[0],4)*this->q[1] - this->q[0] * this->q[0]*pow(this->q[1],3) + this->q[0] * this->q[0]*this->q[1]*this->q[3] * this->q[3] + 2*pow(this->q[0],4)*this->kg[1] -
                                                       3*this->q[0] * this->q[0]*this->q[1] * this->q[1]*this->kg[1] + pow(this->q[1],4)*this->kg[1] + this->q[1] * this->q[1]*this->q[3] * this->q[3]*this->kg[1] - 2*this->q[0] * this->q[0]*this->q[1]*this->kg[1] * this->kg[1] +
                                                       2*pow(this->q[1],3)*this->kg[1] * this->kg[1] - 2*this->q[1]*this->q[3] * this->q[3]*this->kg[1] * this->kg[1] + 4*pow(p0,2)*this->q[0] * this->q[0]*(this->q[1] + this->kg[1]) +
                                                       3*this->q[0] * this->q[0]*this->q[1]*this->kg[2] * this->kg[2] - pow(this->q[1],3)*this->kg[2] * this->kg[2] - 3*this->q[1]*this->q[3] * this->q[3]*this->kg[2] * this->kg[2] -
                                                       2*this->q[1] * this->q[1]*this->kg[1]*this->kg[2] * this->kg[2] - 3*this->q[0] * this->q[0]*this->q[1]*this->q[3]*this->kg[3] + pow(this->q[1],3)*this->q[3]*this->kg[3] + this->q[1]*pow(this->q[3],3)*this->kg[3] -
                                                       2*this->q[0] * this->q[0]*this->q[3]*this->kg[1]*this->kg[3] + 4*this->q[1] * this->q[1]*this->q[3]*this->kg[1]*this->kg[3] - 2*this->q[1]*this->q[3]*this->kg[2] * this->kg[2]*this->kg[3] +
                                                       8*this->mn2*(2*p0*this->q[0]*(this->q[1] + this->kg[1]) - this->q[0] * this->q[0]*(2*this->q[1] + this->kg[1]) +
                                                                    this->kg[1]*(this->q[1] * this->q[1] + this->q[3] * this->q[3] + 2*this->q[1]*this->kg[1] + 2*this->q[3]*this->kg[3])) +
                                                       2*p0*this->q[0]*(pow(this->q[1],3) + 3*this->q[1] * this->q[1]*this->kg[1] - 2*this->q[0] * this->q[0]*(this->q[1] + 2*this->kg[1]) +
                                                                        this->q[1]*(this->q[3] * this->q[3] + this->kg[1] * this->kg[1] - 2*this->kg[2] * this->kg[2] + 3*this->q[3]*this->kg[3]) + this->kg[1]*(this->kg[1] * this->kg[1] + this->kg[2] * this->kg[2] + this->kg[3]*(this->q[3] + this->kg[3]))))))
               )/(24.*mDelta2*pow(this->mn,4));


    tr[3][0] = -(C3v*(this->c_i*C5aNC*this->q[1]*this->kg[2]*(8*pow(mDelta,3)*this->mn*(-p0 + this->q[0]) +
                                                              2*mDelta*this->mn*(p0 - this->q[0])*(8*this->mn2 - 2*p0*this->q[0] + this->q[0] * this->q[0] - this->q[1]*this->kg[1] - this->q[3]*this->kg[3]) -
                                                              4*this->mn2*(p0 - this->q[0])*(2*p0*this->q[0] - this->q[0] * this->q[0] + this->q[1]*this->kg[1] + this->q[3]*this->kg[3]) +
                                                              mDelta2*(-8*pow(p0,3) + 12*pow(p0,2)*this->q[0] + 8*this->mn2*(2*p0 + this->q[0]) +
                                                                       this->q[0]*(this->q[1] * this->q[1] + this->q[3] * this->q[3] - this->q[1]*this->kg[1] - this->q[3]*this->kg[3]) +
                                                                       2*p0*(-2*this->q[0] * this->q[0] + this->q[1] * this->q[1] + this->q[3] * this->q[3] + 2*this->q[1]*this->kg[1] + 2*this->q[3]*this->kg[3]))) +
                      2*C3vNC*this->mn*(16*mDelta2*this->mn*(p0 - this->q[0])*this->q[0]*this->q[3] +
                                        4*pow(mDelta,3)*(this->q[3]*(4*p0*this->q[0] + this->kg[1] * this->kg[1] + this->kg[2] * this->kg[2]) - this->q[1]*this->kg[1]*this->kg[3]) +
                                        4*this->mn*(p0 - this->q[0])*(-4*pow(p0,2)*this->q[0]*this->q[3] + 8*p0*this->q[0] * this->q[0]*this->q[3] - 3*pow(this->q[0],3)*this->q[3] -
                                                                      2*p0*this->q[3]*(this->kg[1] * this->kg[1] + this->kg[2] * this->kg[2]) + this->q[0]*this->q[3]*(this->q[1] * this->q[1] + this->q[3] * this->q[3] + this->q[1]*this->kg[1] + this->kg[1] * this->kg[1] + this->kg[2] * this->kg[2]) +
                                                                      2*p0*this->q[1]*this->kg[1]*this->kg[3] + this->q[0]*(this->q[1] * this->q[1] + 2*this->q[3] * this->q[3] - this->q[1]*this->kg[1])*this->kg[3]) +
                                        mDelta*(-16*pow(p0,3)*this->q[0]*this->q[3] + 4*pow(p0,2)*(8*this->q[0] * this->q[0]*this->q[3] - 3*this->q[3]*(this->kg[1] * this->kg[1] + this->kg[2] * this->kg[2]) + 3*this->q[1]*this->kg[1]*this->kg[3]) +
                                                (this->q[3]*(this->kg[1] * this->kg[1] + this->kg[2] * this->kg[2]) - this->q[1]*this->kg[1]*this->kg[3])*(-3*this->q[0] * this->q[0] + this->q[1] * this->q[1] + 2*this->q[1]*this->kg[1] + this->q[3]*(this->q[3] + 2*this->kg[3])) +
                                                4*p0*this->q[0]*(-3*this->q[0] * this->q[0]*this->q[3] + this->q[1]*this->kg[1]*(this->q[3] - 3*this->kg[3]) + this->q[1] * this->q[1]*(this->q[3] + this->kg[3]) +
                                                                 this->q[3]*(3*(this->kg[1] * this->kg[1] + this->kg[2] * this->kg[2]) + this->q[3]*(this->q[3] + 2*this->kg[3])))))))/(24.*mDelta2*pow(this->mn,4));

    tr[3][1] = (C3v*(this->c_i*C5aNC*this->kg[2]*(-2*this->mn2*this->q[1]*(this->q[1] + this->kg[1])*(2*p0*this->q[0] - this->q[0] * this->q[0] + this->q[1]*this->kg[1] + this->q[3]*this->kg[3]) +
                                                  4*pow(mDelta,3)*this->mn*(8*this->mn2 + 2*(p0 - this->q[0])*this->q[0] + this->q[3]*(this->q[3] + this->kg[3])) +
                                                  mDelta*this->mn*(4*pow(p0,2)*this->q[0] * this->q[0] - 6*p0*pow(this->q[0],3) + 2*pow(this->q[0],4) + 4*p0*this->q[0]*this->q[3] * this->q[3] -
                                                                   this->q[0] * this->q[0]*this->q[3] * this->q[3] + this->q[1] * this->q[1]*this->q[3] * this->q[3] + pow(this->q[3],4) + 2*p0*this->q[0]*this->q[1]*this->kg[1] - 2*this->q[0] * this->q[0]*this->q[1]*this->kg[1] +
                                                                   2*this->q[1]*this->q[3] * this->q[3]*this->kg[1] - 2*this->q[3] * this->q[3]*this->kg[1] * this->kg[1] - 2*this->q[3] * this->q[3]*this->kg[2] * this->kg[2] +
                                                                   this->q[3]*(6*p0*this->q[0] - 5*this->q[0] * this->q[0] + this->q[1] * this->q[1] + 3*this->q[3] * this->q[3] + 2*this->q[1]*this->kg[1])*this->kg[3] +
                                                                   8*this->mn2*(2*p0*this->q[0] - this->q[0] * this->q[0] + this->q[1] * this->q[1] + this->q[3] * this->q[3] + 2*this->q[1]*this->kg[1] + 2*this->q[3]*this->kg[3])) +
                                                  mDelta2*(8*pow(p0,3)*this->q[0] + 8*this->mn2*(4*pow(p0,2) - 2*p0*this->q[0] - this->q[1]*this->kg[1]) +
                                                           2*p0*this->q[0]*(2*this->q[0] * this->q[0] - this->q[1]*(this->q[1] + 2*this->kg[1]) - this->q[3]*(this->q[3] + this->kg[3])) + 4*pow(p0,2)*(-3*this->q[0] * this->q[0] + this->q[3]*(this->q[3] + this->kg[3])) +
                                                           this->q[1]*(this->q[0] * this->q[0]*(-this->q[1] + this->kg[1]) + this->q[3]*(-(this->q[3]*this->kg[1]) + this->q[1]*this->kg[3])))) +
                     2*C3vNC*this->mn*(16*mDelta2*this->mn*(-p0 + this->q[0])*this->q[3]*this->kg[1] + 4*pow(mDelta,3)*(this->q[0]*this->q[1]*this->q[3] - 4*p0*this->q[3]*this->kg[1] + this->q[0]*this->kg[1]*(-2*this->q[3] + this->kg[3])) -
                                       4*this->mn*(p0 - this->q[0])*(2*p0*this->q[0]*(this->q[1]*this->q[3] + this->kg[1]*(2*this->q[3] + this->kg[3])) - this->q[0] * this->q[0]*(this->q[1]*(this->q[3] - this->kg[3]) + this->kg[1]*(2*this->q[3] + this->kg[3])) +
                                                                     this->q[3]*(this->q[1] * this->q[1]*this->kg[1] + this->q[1]*(2*this->kg[1] * this->kg[1] + this->kg[2] * this->kg[2]) + this->q[3]*this->kg[1]*(this->q[3] + 2*this->kg[3]))) +
                                       mDelta*(-4*pow(p0,2)*this->q[0]*(3*this->q[1]*this->q[3] + 2*this->q[3]*this->kg[1] + 3*this->kg[1]*this->kg[3]) +
                                               this->q[0]*(-3*this->q[0] * this->q[0] + this->q[1] * this->q[1] + this->q[3] * this->q[3] + 2*this->q[1]*this->kg[1] + 2*this->q[3]*this->kg[3])*(this->q[1]*this->q[3] + this->kg[1]*(-2*this->q[3] + this->kg[3])) -
                                               4*p0*(this->q[0] * this->q[0]*(-3*this->q[1]*this->q[3] + 2*this->q[3]*this->kg[1] + this->q[1]*this->kg[3] - 3*this->kg[1]*this->kg[3]) +
                                                     this->q[3]*(this->q[1] * this->q[1]*this->kg[1] + 2*this->q[1]*(this->kg[1] * this->kg[1] + this->kg[2] * this->kg[2]) + this->q[3]*this->kg[1]*(this->q[3] + 2*this->kg[3])))))))/(24.*mDelta2*pow(this->mn,4));

    tr[3][2] = (C3v*(-2*C3vNC*this->mn*this->kg[2]*(16*mDelta2*this->mn*(p0 - this->q[0])*this->q[3] + 4*pow(mDelta,3)*(4*p0*this->q[3] + 2*this->q[0]*this->q[3] - this->q[0]*this->kg[3]) +
                                                    4*this->mn*(p0 - this->q[0])*(this->q[3]*(4*p0*this->q[0] - 2*this->q[0] * this->q[0] + this->q[3] * this->q[3] + this->q[1]*(this->q[1] + this->kg[1])) +
                                                                                  (2*p0*this->q[0] - this->q[0] * this->q[0] + this->q[1] * this->q[1] + 2*this->q[3] * this->q[3])*this->kg[3]) +
                                                    mDelta*(4*pow(p0,2)*this->q[0]*(2*this->q[3] + 3*this->kg[3]) + this->q[0]*(2*this->q[3] - this->kg[3])*
                                                                                                                    (-3*this->q[0] * this->q[0] + this->q[1] * this->q[1] + this->q[3] * this->q[3] + 2*this->q[1]*this->kg[1] + 2*this->q[3]*this->kg[3]) +
                                                            4*p0*(this->q[0] * this->q[0]*(2*this->q[3] - 3*this->kg[3]) + (this->q[1] * this->q[1] + this->q[3] * this->q[3])*(this->q[3] + 2*this->kg[3])))) -
                     this->c_i*C5aNC*(2*this->mn2*this->q[1]*this->kg[2] * this->kg[2]*(2*p0*this->q[0] - this->q[0] * this->q[0] + this->q[1]*this->kg[1] + this->q[3]*this->kg[3]) +
                                      4*pow(mDelta,3)*this->mn*(8*this->mn2*this->kg[1] + 2*p0*this->q[0]*this->kg[1] - 2*this->q[0] * this->q[0]*this->kg[1] + this->q[1] * this->q[1]*this->kg[1] + this->q[3] * this->q[3]*this->kg[1] +
                                                                this->q[1]*this->kg[1] * this->kg[1] + this->q[1]*this->kg[2] * this->kg[2] + this->q[3]*this->kg[1]*this->kg[3]) +
                                      mDelta2*(8*pow(p0,3)*this->q[0]*this->kg[1] + this->q[1]*(-this->q[0] * this->q[0] + this->q[1] * this->q[1] + this->q[3] * this->q[3])*this->kg[2] * this->kg[2] +
                                               8*this->mn2*(2*p0*this->q[0]*(this->q[1] - this->kg[1]) + 4*pow(p0,2)*this->kg[1] + this->q[1]*this->kg[2] * this->kg[2]) +
                                               4*pow(p0,2)*(this->q[0] * this->q[0]*(this->q[1] - 3*this->kg[1]) + this->q[1] * this->q[1]*this->kg[1] + this->q[1]*(this->kg[1] * this->kg[1] + this->kg[2] * this->kg[2]) + this->q[3]*this->kg[1]*(this->q[3] + this->kg[3])) -
                                               2*p0*this->q[0]*(-pow(this->q[1],3) + 2*this->q[0] * this->q[0]*(this->q[1] - this->kg[1]) + this->q[3]*this->kg[1]*(this->q[3] + this->kg[3]) -
                                                                this->q[1]*(-this->kg[1] * this->kg[1] + this->kg[2] * this->kg[2] + this->q[3]*(this->q[3] + this->kg[3])))) +
                                      mDelta*this->mn*(4*pow(this->q[0],4)*this->q[1] - 2*this->q[0] * this->q[0]*pow(this->q[1],3) - 2*this->q[0] * this->q[0]*this->q[1]*this->q[3] * this->q[3] + 2*pow(this->q[0],4)*this->kg[1] -
                                                       5*this->q[0] * this->q[0]*this->q[1] * this->q[1]*this->kg[1] + pow(this->q[1],4)*this->kg[1] - this->q[0] * this->q[0]*this->q[3] * this->q[3]*this->kg[1] + 2*this->q[1] * this->q[1]*this->q[3] * this->q[3]*this->kg[1] +
                                                       pow(this->q[3],4)*this->kg[1] - 5*this->q[0] * this->q[0]*this->q[1]*this->kg[1] * this->kg[1] + 3*pow(this->q[1],3)*this->kg[1] * this->kg[1] + 3*this->q[1]*this->q[3] * this->q[3]*this->kg[1] * this->kg[1] +
                                                       2*this->q[1] * this->q[1]*pow(this->kg[1],3) - 2*this->q[3] * this->q[3]*pow(this->kg[1],3) + 4*pow(p0,2)*this->q[0] * this->q[0]*(this->q[1] + this->kg[1]) -
                                                       3*this->q[0] * this->q[0]*this->q[1]*this->kg[2] * this->kg[2] + pow(this->q[1],3)*this->kg[2] * this->kg[2] + this->q[1]*this->q[3] * this->q[3]*this->kg[2] * this->kg[2] +
                                                       2*this->q[1] * this->q[1]*this->kg[1]*this->kg[2] * this->kg[2] - 2*this->q[3] * this->q[3]*this->kg[1]*this->kg[2] * this->kg[2] +
                                                       this->q[3]*(-(this->q[0] * this->q[0]*(2*this->q[1] + 5*this->kg[1])) + this->kg[1]*(3*(this->q[1] * this->q[1] + this->q[3] * this->q[3]) + 4*this->q[1]*this->kg[1]) + 2*this->q[1]*this->kg[2] * this->kg[2])*this->kg[3] +
                                                       8*this->mn2*(2*p0*this->q[0]*(this->q[1] + this->kg[1]) - this->q[0] * this->q[0]*(2*this->q[1] + this->kg[1]) +
                                                                    this->kg[1]*(this->q[1] * this->q[1] + this->q[3] * this->q[3] + 2*this->q[1]*this->kg[1] + 2*this->q[3]*this->kg[3])) +
                                                       2*p0*this->q[0]*(pow(this->q[1],3) + 3*this->q[1] * this->q[1]*this->kg[1] - this->q[0] * this->q[0]*(4*this->q[1] + 3*this->kg[1]) + this->q[3]*this->kg[1]*(2*this->q[3] + 3*this->kg[3]) +
                                                                        this->q[1]*(3*this->kg[1] * this->kg[1] + 2*this->kg[2] * this->kg[2] + this->q[3]*(this->q[3] + this->kg[3])))))))/(24.*mDelta2*pow(this->mn,4));

    tr[3][3] = -(C3v*(this->c_i*C5aNC*this->q[1]*this->kg[2]*(4*pow(mDelta,3)*this->mn*(this->q[3] + this->kg[3]) +
                                                              2*this->mn2*(this->q[3] + this->kg[3])*(2*p0*this->q[0] - this->q[0] * this->q[0] + this->q[1]*this->kg[1] + this->q[3]*this->kg[3]) +
                                                              mDelta*this->mn*(this->q[3] + this->kg[3])*(4*p0*this->q[0] - 3*this->q[0] * this->q[0] + this->q[1] * this->q[1] + this->q[3] * this->q[3] + 2*this->q[1]*this->kg[1] + 2*this->q[3]*this->kg[3]) +
                                                              mDelta2*(-(this->q[1]*this->q[3]*this->kg[1]) + this->q[0] * this->q[0]*(this->q[3] - this->kg[3]) + 8*this->mn2*this->kg[3] + 2*p0*this->q[0]*this->kg[3] + this->q[1] * this->q[1]*this->kg[3] +
                                                                       4*pow(p0,2)*(this->q[3] + this->kg[3]))) + 2*C3vNC*this->mn*
                                                                                                                  (16*mDelta2*this->mn*(p0 - this->q[0])*(this->q[0] * this->q[0] - this->q[1]*this->kg[1]) +
                                                                                                                   4*pow(mDelta,3)*(4*p0*(this->q[0] * this->q[0] - this->q[1]*this->kg[1]) + this->q[0]*(pow(this->q[1] - this->kg[1],2) + this->kg[2] * this->kg[2])) -
                                                                                                                   4*this->mn*(p0 - this->q[0])*(4*pow(p0,2)*this->q[0] * this->q[0] + 3*pow(this->q[0],4) +
                                                                                                                                                 2*p0*this->q[0]*(-4*this->q[0] * this->q[0] + pow(this->q[1] + this->kg[1],2) + this->kg[2] * this->kg[2]) -
                                                                                                                                                 this->q[0] * this->q[0]*(this->q[3] * this->q[3] + (this->q[1] + this->kg[1])*(2*this->q[1] + this->kg[1]) + this->kg[2] * this->kg[2] + 2*this->q[3]*this->kg[3]) +
                                                                                                                                                 this->q[1]*(this->q[1] * this->q[1]*this->kg[1] + this->q[1]*(2*this->kg[1] * this->kg[1] + this->kg[2] * this->kg[2]) + this->q[3]*this->kg[1]*(this->q[3] + 2*this->kg[3]))) +
                                                                                                                   mDelta*(-16*pow(p0,3)*this->q[0] * this->q[0] + 4*pow(p0,2)*this->q[0]*
                                                                                                                                                                   (8*this->q[0] * this->q[0] - 3*this->q[1] * this->q[1] - 2*this->q[1]*this->kg[1] - 3*(this->kg[1] * this->kg[1] + this->kg[2] * this->kg[2])) +
                                                                                                                           this->q[0]*(pow(this->q[1] - this->kg[1],2) + this->kg[2] * this->kg[2])*(-3*this->q[0] * this->q[0] + this->q[1] * this->q[1] + this->q[3] * this->q[3] + 2*this->q[1]*this->kg[1] + 2*this->q[3]*this->kg[3]) -
                                                                                                                           4*p0*(3*pow(this->q[0],4) - this->q[0] * this->q[0]*(4*this->q[1] * this->q[1] - this->q[1]*this->kg[1] + 3*(this->kg[1] * this->kg[1] + this->kg[2] * this->kg[2]) + this->q[3]*(this->q[3] + 2*this->kg[3])) +
                                                                                                                                 this->q[1]*(this->q[1] * this->q[1]*this->kg[1] + 2*this->q[1]*(this->kg[1] * this->kg[1] + this->kg[2] * this->kg[2]) + this->q[3]*this->kg[1]*(this->q[3] + 2*this->kg[3])))))))/
               (24.*mDelta2*pow(this->mn,4));

}

double Hadronic_Current_R_Delta::qcm(double s) {
//    Returns the 3-momentum of the pion formed after the decay of a
//    resonance (R->N pi) of inv. mass s in the rest frame of the resonance
    return sqrt((s-this->mpi2-this->mn2)*(s-this->mpi2-this->mn2) - 4.0 * this->mpi2 * this->mn2)/2.0/sqrt(s);
}

std::complex<double> Hadronic_Current_R_Delta::getR(int i, int j) {
    std::complex<double> tr_p = this->tr_p_dir[i][j]*this->propagator_dir_p
                                + this->tr_p_crs[i][j]*this->propagator_crs_p;
    std::complex<double> tr_n = this->tr_n_dir[i][j]*this->propagator_dir_n
                                + this->tr_n_crs[i][j]*this->propagator_crs_n;

    return (tr_p * this->N_FF_p + tr_n * this->N_FF_n) * this->mn / this->p0 / 2.0; //Factor 1/2 from the trace, see notes
}

void Hadronic_Current_R_Delta::setP(std::vector<double> p_in) {
    this->p = p_in;
    this->p0 = p_in[0];

    double aux_pp[] = {p_in[0], -p_in[1], -p_in[2], -p_in[3]};
    this->pp.assign(aux_pp, aux_pp + sizeof(aux_pp) / sizeof(double) );

    this->p_dir = this->p;
    std::transform(this->p_dir.begin( ), this->p_dir.end( ), this->q.begin( ), this->p_dir.begin( ),std::plus<double>( ));

    this->p_crs = this->pp;
    std::transform(this->p_crs.begin( ), this->p_crs.end( ), this->q.begin( ), this->p_crs.begin( ),std::minus<double>( ));
}

std::complex<double> Hadronic_Current_R_Delta::Propagator_Delta(std::vector<double> p_in) {
    double p2 = p_in[0] * p_in[0] - p_in[1] * p_in[1] - p_in[2] * p_in[2] - p_in[3] * p_in[3];

    std::vector<std::complex<double> > param(2);
    param[0] = p2;
    param[1] = this->q_minus_kg;
    this->propagator->change_other_parameters(param);

    if (this->dir_or_crs == "dir") {
        this->propagator->setNucleon("p");
        this->propagator_dir_p = this->propagator->Delta_propagator_avg_dir();
        this->propagator->setNucleon("n");
        this->propagator_dir_n = this->propagator->Delta_propagator_avg_dir();
    } else if (this->dir_or_crs == "crs") {
        this->propagator->setNucleon("p");
        this->propagator_crs_p = this->propagator->Delta_propagator_crs();
        this->propagator->setNucleon("n");
        this->propagator_crs_n = this->propagator->Delta_propagator_crs();
    }

    return 0;
}




//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////




