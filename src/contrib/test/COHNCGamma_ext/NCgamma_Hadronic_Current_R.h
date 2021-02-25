//
// Created by edusaul on 11/03/19.
//

#ifndef COHERENT_PHOTON_HADRONIC_CURRENT_R_H
#define COHERENT_PHOTON_HADRONIC_CURRENT_R_H

#include <vector>
#include <complex>
#include "NCgamma_Form_Factors.h"



namespace NC_gamma {

    class Delta_in_medium {
    protected:
        std::string nucleon;
        std::string dir_or_crs;
        double p;
        double p2;
        double q_minus_kg;
        double qcm;
        double mDelta;
        double mDelta2;
        double mn;
        double mn2;
        double mpi;
        double V0;
        double rho0;
        double coupling_DeltaNpi;
        std::complex<double> c_i;
        double hc3;

        double gamma_vac;
        double rho_r;

        double rho_avg; // average density in units of rho0
        std::complex<double> propagator_avg;

        std::complex<double> propagator_crs;

    public:

        Delta_in_medium();

        virtual ~Delta_in_medium();

        void setNucleon(const std::string &n);

        virtual std::complex<double> Delta_propagator_dir(double);

        virtual std::complex<double> Delta_propagator_crs();

        void change_other_parameters(std::vector<std::complex<double> >);

        std::complex<double> Sigma();

        double Gamma_tilde_Delta(double);

        double Gamma_vacuum(double);

        double I_series(double);

        double lambda_func(double, double, double);

        double kf(double);

        void setDirOrCrs(const std::string &dirOrCrs);

        std::complex<double> Delta_propagator_avg_dir();

    };


//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////



    class Hadronic_Current_R {
    protected:
        std::complex<double> c_i;

        std::vector<double> q;
        std::vector<double> kg;
        std::vector<double> p;
        std::vector<double> pp;
        double p0;

        typedef std::complex<double> Array4x4[4][4];

        Array4x4 tr_p_dir;
        Array4x4 tr_p_crs;
        Array4x4 tr_n_dir;
        Array4x4 tr_n_crs;

        virtual void set_tr_dir(Array4x4 &);

        virtual void set_tr_crs(Array4x4 &);

    public:

        explicit Hadronic_Current_R();

        virtual void setQ(const std::vector<double> &q_in);

        virtual void setKg(const std::vector<double> &kg_in);

        virtual void setP(std::vector<double>);

        virtual std::complex<double> getR(int, int);

        virtual void setFF(double);

        virtual ~Hadronic_Current_R();
    };


//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////


    class Hadronic_Current_R_Sum : public Hadronic_Current_R {
    protected:
        std::vector<Hadronic_Current_R *> vector_of_currents;

    public:
        explicit Hadronic_Current_R_Sum(const std::vector<Hadronic_Current_R *> &vectorOfCurrents);

        virtual void setQ(const std::vector<double> &q);

        virtual void setKg(const std::vector<double> &kg);

        virtual void setP(std::vector<double>);

        virtual std::complex<double> getR(int, int);

        virtual void setFF(double);

    };



//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////


    class Hadronic_Current_R_Delta : public Hadronic_Current_R {
    protected:
        Form_Factors_Delta *ff_Delta;

        double mn;
        double mn2;
        double mDelta;
        double mDelta2;
        double mpi;
        double mpi2;
        std::complex<double> c_i;
        double q_minus_kg;

        double C3v;
        double C4v;
        double C5v;
        double C3vNC;
        double C4vNC;
        double C5vNC;
        double C3aNC;
        double C4aNC;
        double C5aNC;
        double C6aNC;

        std::complex<double> N_FF_p;
        std::complex<double> N_FF_n;

        std::complex<double> propagator_dir_p;
        std::complex<double> propagator_crs_p;
        std::complex<double> propagator_dir_n;
        std::complex<double> propagator_crs_n;

        std::vector<double> p_dir;
        std::vector<double> p_crs;

        Delta_in_medium *propagator;
        std::string dir_or_crs;

        Nuclear_FF *nuclearFF;


        virtual void set_tr_dir(Array4x4 &);

        virtual void set_tr_crs(Array4x4 &);

    public:

        explicit Hadronic_Current_R_Delta(Nuclear_FF *nuclearFf);

        virtual ~Hadronic_Current_R_Delta();

        virtual void setFF(double);

        virtual void setKg(const std::vector<double> &kg_in);

        std::complex<double> Propagator_Delta(std::vector<double>);

        double qcm(double);

        virtual std::complex<double> getR(int, int);

        virtual void setP(std::vector<double>);

    };



//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

}

#endif //COHERENT_PHOTON_HADRONIC_CURRENT_R_H
