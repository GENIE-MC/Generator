//
// Created by edusaul on 23/04/19.
//

#include <utility>
#include "NCgamma_Diff_Cross_Section.h"
#include "NCgamma_Form_Factors.h"
#include "NCgamma_Parameters_GeV.h"
#include <algorithm>
#include <iostream>
#include <math.h>

using namespace NC_gamma;

Diff_Cross_Section::Diff_Cross_Section(const std::string &m, const std::string &n) : 
  mode(m), nucleus(n) {

  this->c_i = NCgamma_Param::c_i;
  double e2 = 4.0 * NCgamma_Param::pi * NCgamma_Param::alpha;
  
  this->constant_factors = NCgamma_Param::Gf2 * e2 / (2.0 * 8.0 * pow(2.0*NCgamma_Param::pi,4) );// * NCgamma_Param::hccm2;
    
  this->nuclearFF = new Nucleus_FF_DeVries(nucleus);
  this->vector_of_currents.push_back(new Hadronic_Current_R_Delta(nuclearFF));
  this->current_R = new Hadronic_Current_R_Sum(this->vector_of_currents);

}


Diff_Cross_Section::~Diff_Cross_Section() {
    delete this->nuclearFF;
    for (int i = 0; i <  static_cast<unsigned int>(this->vector_of_currents.size()); ++i) {
        delete this->vector_of_currents[i];
    }
    delete this->current_R;
}

double Diff_Cross_Section::getDiffCrossSection(double Enu, 
					       double Eg, 
					       double theta_l, 
					       double theta_g, double phi_g) {

  double params_aux[]= {Enu, Eg, phi_g, theta_g};
  std::vector<double> param (params_aux, params_aux + sizeof(params_aux) / sizeof(double) );
  
  this->change_other_parameters(param);
  
  double cs =this->integrand(theta_l);
  
  if(std::isnan(cs)) return 0.0;
  else return cs;

}



double Diff_Cross_Section::integrand(double th) {

    double kp_aux[] = {kp0, kp0 * sin(th), 0., kp0 * cos(th)};
    this->kp.assign(kp_aux, kp_aux + sizeof(kp_aux) / sizeof(double) );

    this->q = this->k;
    std::transform(q.begin( ), q.end( ), kp.begin( ), q.begin( ),std::minus<double>( ));

    this->current_R->setQ(this->q);

    double kg_aux[] = {this->k0g, this->k0g * sin(this->thg) * cos(this->phig),
                       this->k0g * sin(this->thg) * sin(this->phig), this->k0g * cos(this->thg)};
    this->kg.assign(kg_aux, kg_aux + sizeof(kg_aux) / sizeof(double) );

    std::vector<double> p(this->kg.size());
    for (int j = 0; j < static_cast<unsigned int>(this->kg.size()); ++j) {
        p[j] = (this->kg[j] - this->q[j])/2.0;
    }

    double pv2 = p[1]*p[1] + p[2]*p[2] + p[3]*p[3];
    p[0] = sqrt(NCgamma_Param::mp2 + pv2);

    this->current_R->setP(p);
    this->current_R->setKg(this->kg);

    double lh;
    if(this->mode == "nu") lh = this->LH_contraction_neutrino();
    else lh = this->LH_contraction_antineutrino();

    return lh*this->constant_factors*this->factors;
}

void Diff_Cross_Section::change_other_parameters(std::vector<double> k0_k0g_phig_thg) {
    this->k0 = k0_k0g_phig_thg[0];
    this->k0g = k0_k0g_phig_thg[1];
    this->phig = k0_k0g_phig_thg[2];
    this->thg = k0_k0g_phig_thg[3];

    //q0 = k0g approx
    this->kp0 = k0 - this->k0g;

    double k_aux[] = {k0, 0., 0., k0};
    this->k.assign(k_aux, k_aux + sizeof(k_aux) / sizeof(double) );

    this->factors = this->k0g * (this->k0-this->k0g) /this->k0;

}

double Diff_Cross_Section::LH_contraction_neutrino() {
    std::complex<double> lh = -8*k[0]*((2*k[0] - this->q[0] - this->q[3])*H(0,0,0,0) + this->q[1]*H(0,0,0,1) - this->c_i*this->q[1]*H(0,0,0,2) +
                       (-2*k[0] + this->q[0] + this->q[3])*H(0,0,0,3) + this->q[1]*H(0,1,0,0) + (-this->q[0] + this->q[3])*H(0,1,0,1) +
                       this->c_i*(this->q[0] - this->q[3])*H(0,1,0,2) - this->q[1]*H(0,1,0,3) + this->c_i*this->q[1]*H(0,2,0,0) -
                       this->c_i*(this->q[0] - this->q[3])*H(0,2,0,1) + (-this->q[0] + this->q[3])*H(0,2,0,2) - this->c_i*this->q[1]*H(0,2,0,3) +
                       (-2*k[0] + this->q[0] + this->q[3])*H(0,3,0,0) - this->q[1]*H(0,3,0,1) + this->c_i*this->q[1]*H(0,3,0,2) +
                       (2*k[0] - this->q[0] - this->q[3])*H(0,3,0,3) + (-2*k[0] + this->q[0] + this->q[3])*H(1,0,1,0) - this->q[1]*H(1,0,1,1) +
                       this->c_i*this->q[1]*H(1,0,1,2) + (2*k[0] - this->q[0] - this->q[3])*H(1,0,1,3) - this->q[1]*H(1,1,1,0) +
                       (this->q[0] - this->q[3])*H(1,1,1,1) - this->c_i*(this->q[0] - this->q[3])*H(1,1,1,2) + this->q[1]*H(1,1,1,3) -
                       this->c_i*this->q[1]*H(1,2,1,0) + this->c_i*(this->q[0] - this->q[3])*H(1,2,1,1) + (this->q[0] - this->q[3])*H(1,2,1,2) +
                       this->c_i*this->q[1]*H(1,2,1,3) + (2*k[0] - this->q[0] - this->q[3])*H(1,3,1,0) + this->q[1]*H(1,3,1,1) -
                       this->c_i*this->q[1]*H(1,3,1,2) + (-2*k[0] + this->q[0] + this->q[3])*H(1,3,1,3) + (-2*k[0] + this->q[0] + this->q[3])*H(2,0,2,0) -
                       this->q[1]*H(2,0,2,1) + this->c_i*this->q[1]*H(2,0,2,2) + (2*k[0] - this->q[0] - this->q[3])*H(2,0,2,3) - this->q[1]*H(2,1,2,0) +
                       (this->q[0] - this->q[3])*H(2,1,2,1) - this->c_i*(this->q[0] - this->q[3])*H(2,1,2,2) + this->q[1]*H(2,1,2,3) -
                       this->c_i*this->q[1]*H(2,2,2,0) + this->c_i*(this->q[0] - this->q[3])*H(2,2,2,1) + (this->q[0] - this->q[3])*H(2,2,2,2) +
                       this->c_i*this->q[1]*H(2,2,2,3) + (2*k[0] - this->q[0] - this->q[3])*H(2,3,2,0) + this->q[1]*H(2,3,2,1) -
                       this->c_i*this->q[1]*H(2,3,2,2) + (-2*k[0] + this->q[0] + this->q[3])*H(2,3,2,3) + (-2*k[0] + this->q[0] + this->q[3])*H(3,0,3,0) -
                       this->q[1]*H(3,0,3,1) + this->c_i*this->q[1]*H(3,0,3,2) + (2*k[0] - this->q[0] - this->q[3])*H(3,0,3,3) - this->q[1]*H(3,1,3,0) +
                       (this->q[0] - this->q[3])*H(3,1,3,1) - this->c_i*(this->q[0] - this->q[3])*H(3,1,3,2) + this->q[1]*H(3,1,3,3) -
                       this->c_i*this->q[1]*H(3,2,3,0) + this->c_i*(this->q[0] - this->q[3])*H(3,2,3,1) + (this->q[0] - this->q[3])*H(3,2,3,2) +
                       this->c_i*this->q[1]*H(3,2,3,3) + (2*k[0] - this->q[0] - this->q[3])*H(3,3,3,0) + this->q[1]*H(3,3,3,1) -
                       this->c_i*this->q[1]*H(3,3,3,2) + (-2*k[0] + this->q[0] + this->q[3])*H(3,3,3,3));

    return lh.real();
}

double Diff_Cross_Section::LH_contraction_antineutrino() {
    std::complex<double> lh = -8*this->k[0]*((2*this->k[0] - this->q[0] - this->q[3])*H(0,0,0,0) + this->q[1]*H(0,0,0,1) + this->c_i*this->q[1]*H(0,0,0,2) +
                             (-2*this->k[0] + this->q[0] + this->q[3])*H(0,0,0,3) + this->q[1]*H(0,1,0,0) + (-this->q[0] + this->q[3])*H(0,1,0,1) -
                             this->c_i*(this->q[0] - this->q[3])*H(0,1,0,2) - this->q[1]*H(0,1,0,3) - this->c_i*this->q[1]*H(0,2,0,0) +
                             this->c_i*(this->q[0] - this->q[3])*H(0,2,0,1) + (-this->q[0] + this->q[3])*H(0,2,0,2) + this->c_i*this->q[1]*H(0,2,0,3) +
                             (-2*this->k[0] + this->q[0] + this->q[3])*H(0,3,0,0) - this->q[1]*H(0,3,0,1) - this->c_i*this->q[1]*H(0,3,0,2) +
                             (2*this->k[0] - this->q[0] - this->q[3])*H(0,3,0,3) + (-2*this->k[0] + this->q[0] + this->q[3])*H(1,0,1,0) - this->q[1]*H(1,0,1,1) -
                             this->c_i*this->q[1]*H(1,0,1,2) + (2*this->k[0] - this->q[0] - this->q[3])*H(1,0,1,3) - this->q[1]*H(1,1,1,0) +
                             (this->q[0] - this->q[3])*H(1,1,1,1) + this->c_i*(this->q[0] - this->q[3])*H(1,1,1,2) + this->q[1]*H(1,1,1,3) +
                             this->c_i*this->q[1]*H(1,2,1,0) - this->c_i*(this->q[0] - this->q[3])*H(1,2,1,1) + (this->q[0] - this->q[3])*H(1,2,1,2) -
                             this->c_i*this->q[1]*H(1,2,1,3) + (2*this->k[0] - this->q[0] - this->q[3])*H(1,3,1,0) + this->q[1]*H(1,3,1,1) +
                             this->c_i*this->q[1]*H(1,3,1,2) + (-2*this->k[0] + this->q[0] + this->q[3])*H(1,3,1,3) + (-2*this->k[0] + this->q[0] + this->q[3])*H(2,0,2,0) -
                             this->q[1]*H(2,0,2,1) - this->c_i*this->q[1]*H(2,0,2,2) + (2*this->k[0] - this->q[0] - this->q[3])*H(2,0,2,3) - this->q[1]*H(2,1,2,0) +
                             (this->q[0] - this->q[3])*H(2,1,2,1) + this->c_i*(this->q[0] - this->q[3])*H(2,1,2,2) + this->q[1]*H(2,1,2,3) +
                             this->c_i*this->q[1]*H(2,2,2,0) - this->c_i*(this->q[0] - this->q[3])*H(2,2,2,1) + (this->q[0] - this->q[3])*H(2,2,2,2) -
                             this->c_i*this->q[1]*H(2,2,2,3) + (2*this->k[0] - this->q[0] - this->q[3])*H(2,3,2,0) + this->q[1]*H(2,3,2,1) +
                             this->c_i*this->q[1]*H(2,3,2,2) + (-2*this->k[0] + this->q[0] + this->q[3])*H(2,3,2,3) + (-2*this->k[0] + this->q[0] + this->q[3])*H(3,0,3,0) -
                             this->q[1]*H(3,0,3,1) - this->c_i*this->q[1]*H(3,0,3,2) + (2*this->k[0] - this->q[0] - this->q[3])*H(3,0,3,3) - this->q[1]*H(3,1,3,0) +
                             (this->q[0] - this->q[3])*H(3,1,3,1) + this->c_i*(this->q[0] - this->q[3])*H(3,1,3,2) + this->q[1]*H(3,1,3,3) +
                             this->c_i*this->q[1]*H(3,2,3,0) - this->c_i*(this->q[0] - this->q[3])*H(3,2,3,1) + (this->q[0] - this->q[3])*H(3,2,3,2) -
                             this->c_i*this->q[1]*H(3,2,3,3) + (2*this->k[0] - this->q[0] - this->q[3])*H(3,3,3,0) + this->q[1]*H(3,3,3,1) +
                             this->c_i*this->q[1]*H(3,3,3,2) + (-2*this->k[0] + this->q[0] + this->q[3])*H(3,3,3,3));

    return lh.real();
}


std::complex<double> Diff_Cross_Section::H(int l, int m, int n, int o) {
    std::complex<double> r1c = std::conj(this->current_R->getR(l,m));
    std::complex<double> r1 = this->current_R->getR(n,o);
    std::complex<double> r = r1 * r1c;
    return r;
}

void Diff_Cross_Section::setMode(const std::string &mode_input) {
    Diff_Cross_Section::mode = mode_input;
}

void Diff_Cross_Section::setNucleus(const std::string &nucleus_input) {
    Diff_Cross_Section::nucleus = nucleus_input;
}















