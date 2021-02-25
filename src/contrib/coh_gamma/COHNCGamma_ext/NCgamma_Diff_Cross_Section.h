//
// Created by edusaul on 23/04/19.
//

#ifndef GENIE_XSECTION_NCGAMMA_DIFF_CROSS_SECTION_H
#define GENIE_XSECTION_NCGAMMA_DIFF_CROSS_SECTION_H

#include <complex>
#include <vector>
#include "NCgamma_Hadronic_Current_R.h"

namespace NC_gamma {

    class Diff_Cross_Section {
    protected:
        std::string mode;
        std::string nucleus;
        std::vector<double> k;
        std::vector<double> kp;
        std::vector<double> q;
        std::vector<double> kg;

        std::complex<double> c_i;

        double k0;
        double k0g;
        double thg;
        double phig;
        double kp0;
        double constant_factors;
        double factors;

        Nucleus_FF_DeVries *nuclearFF;
        std::vector<Hadronic_Current_R *> vector_of_currents;
        Hadronic_Current_R *current_R;


    public:
        Diff_Cross_Section(const std::string &m, const std::string &n);

        virtual ~Diff_Cross_Section();

        double getDiffCrossSection(double Enu, double Enu_final, double theta_l, double theta_g, double phi_g);

        double integrand(double);

        void change_other_parameters(std::vector<double>);

        double LH_contraction_neutrino();

        double LH_contraction_antineutrino();

        std::complex<double> H(int, int, int, int);

        void setMode(const std::string &mode_input);

        void setNucleus(const std::string &nucleus_input);
    };

}

#endif //GENIE_XSECTION_NCGAMMA_DIFF_CROSS_SECTION_H
