//
// Created by edusaul on 8/04/19.
//

#ifndef COHERENT_PHOTON_NUCLEAR_FF_H
#define COHERENT_PHOTON_NUCLEAR_FF_H
#include <complex>
#include <string>

namespace NC_gamma {

    class Nuclear_FF {
    protected:
        std::complex<double> FF_p;
        std::complex<double> FF_n;
    public:

        const std::complex<double> &getFfP() const;

        const std::complex<double> &getFfN() const;

        virtual void setFF(double) = 0;

    };


/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////


    class Nucleus_FF_DeVries : public Nuclear_FF {
    protected:
        std::string nucleus;

        double a[16];
        double R;
        double pi;
        double hc;

    public:
        explicit Nucleus_FF_DeVries(const std::string &n);

        virtual void setFF(double);

        void set_12C();

        void set_40Ar();

        virtual ~Nucleus_FF_DeVries();

    };



/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////



    class Form_Factors_Delta {
    protected:
        double mn;
        double mn2;
        double mDelta;
        double mDelta2;
        double N_Delta_MA;
        double aNC;

        double kgcm0;
        double mpw2;
        double mmw2;

        double C3v;
        double C4v;
        double C5v;

        double C3a;
        double C4a;
        double C5a;
        double C6a;

        double C3vNC;
        double C4vNC;
        double C5vNC;

        double C3aNC;
        double C4aNC;
        double C5aNC;
        double C6aNC;

    public:
        Form_Factors_Delta();

        void setFF(double);

        double getC3v() const;

        double getC4v() const;

        double getC5v() const;

        double getC3a() const;

        double getC4a() const;

        double getC5a() const;

        double getC6a() const;

        double getC3vNC() const;

        double getC4vNC() const;

        double getC5vNC() const;

        double getC3aNC() const;

        double getC4aNC() const;

        double getC5aNC() const;

        double getC6aNC() const;

    };

}



#endif //COHERENT_PHOTON_NUCLEAR_FF_H
