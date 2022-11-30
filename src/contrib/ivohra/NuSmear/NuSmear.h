//____________________________________________________________________________
/*
\name     NuSmear smearing system
\brief    Provides preliminary simulation of energy smearing and angular
          smearing for neutrino-nucleon interactions via parameterized
          model-based presets in the DUNE-CDR and Default models.
          
          To see the full NuSmear paper, visit:
          https://inspirehep.net/literature/2150455
          
\author   Ishaan Vohra <ivohra@exeter.edu / ishaanklv@gmail.com>
          Phillips Exeter Academy
\created  August 16, 2022
*/
//____________________________________________________________________________

#include <iostream>
#include <string>
#include <map>
#include <random>
#include "GHepParticle.h"
#include "PDGCodes.h"
#include "TVector3.h"
#include "PDGUtils.h"

#ifndef _NUSMEAR_H
#define _NUSMEAR_H

using namespace genie;

// DUNE-CDR resolution functions

double PiP_Res_duneCdr(double myE, double myKE, double myPmag)
{
    if (myKE >= 0.1)
    {
        return 0.15;
    }
    else
    {
        return -1;
    }
}

double PiM_Res_duneCdr(double myE, double myKE, double myPmag)
{
    if (myKE >= 0.1)
    {
        return 0.15;
    }
    else
    {
        return -1;
    }
}

double Pi0_Res_duneCdr(double myE, double myKE, double myPmag)
{
    if (myKE >= 0.05)
    {
        return pow((pow(0.05, 2) + pow((0.3 / (pow(myE, 0.5))), 2)), 0.5);
    }
    else
    {
        return -1;
    }
}

double KP_Res_duneCdr(double myE, double myKE, double myPmag)
{
    if (myKE >= 0.05)
    {
        return pow((pow(0.05, 2) + pow((0.3 / (pow(myE, 0.5))), 2)), 0.5);
    }
    else
    {
        return -1;
    }
}

double KM_Res_duneCdr(double myE, double myKE, double myPmag)
{
    if (myKE >= 0.05)
    {
        return pow(pow(0.05, 2) + pow((0.3 / (pow(myE, 0.5))), 2), 0.5);
    }
    else
    {
        return -1;
    }
}

double Gamma_Res_duneCdr(double myE, double myKE, double myPmag)
{
    if (myKE >= 0.03)
    {
        return pow(pow(0.02, 2) + pow(0.15 / (pow(myE, 0.5)), 2), 0.5);
    }
    else
    {
        return -1;
    }
}

double Proton_Res_duneCdr(double myE, double myKE, double myPmag)
{
    if (myKE >= 0.05)
    {
        if (myPmag < 0.4)
        {
            return 0.1;
        }
        else
        {
            return pow(pow(0.05, 2) + pow(0.3 / (pow(myKE, 0.5)), 2), 0.5);
        }
    }
    else
    {
        return -1;
    }
}

double Electron_Res_duneCdr(double myE, double myKE, double myPmag)
{
    if (myKE >= 0.03)
    {
        return pow(pow(0.02, 2) + pow(0.15 / (pow(myE, 0.5)), 2), 0.5);
    }
    else
    {
        return -1;
    }
}

double Positron_Res_duneCdr(double myE, double myKE, double myPmag)
{
    if (myKE >= 0.03)
    {
        return pow(pow(0.02, 2) + pow(0.15 / (pow(myE, 0.5)), 2), 0.5);
    }
    else
    {
        return -1;
    }
}

double Muon_Res_duneCdr(double myE, double myKE, double myPmag)
{
    if (myKE >= 0.03)
    {
        return 0.15;
    }
    else
    {
        return -1;
    }
}

double AntiMuon_Res_duneCdr(double myE, double myKE, double myPmag)
{
    if (myKE >= 0.03)
    {
        return 0.15;
    }
    else
    {
        return -1;
    }
}

double Neutron_Res_duneCdr(double myE, double myKE, double myPmag)
{
    if (myKE >= 0.05)
    {
        return 0.4 / (pow(myKE, 0.5));
    }
    else
    {
        return -1;
    }
}

double K0_Res_duneCdr(double myE, double myKE, double myPmag)
{
    if (myKE >= 0.05)
    {
        return pow((pow(0.05, 2) + pow((0.3 / (pow(myE, 0.5))), 2)), 0.5);
    }
    else
    {
        return -1;
    }
}

double AntiK0_Res_duneCdr(double myE, double myKE, double myPmag)
{
    if (myKE >= 0.05)
    {
        return pow((pow(0.05, 2) + pow((0.3 / (pow(myE, 0.5))), 2)), 0.5);
    }
    else
    {
        return -1;
    }
}

double Lambda_Res_duneCdr(double myE, double myKE, double myPmag)
{
    if (myKE >= 0.05)
    {
        return pow((pow(0.05, 2) + pow((0.3 / (pow(myE, 0.5))), 2)), 0.5);
    }
    else
    {
        return -1;
    }
}

double SigmaP_Res_duneCdr(double myE, double myKE, double myPmag)
{
    if (myKE >= 0.05)
    {
        return pow((pow(0.05, 2) + pow((0.3 / (pow(myE, 0.5))), 2)), 0.5);
    }
    else
    {
        return -1;
    }
}

double Sigma0_Res_duneCdr(double myE, double myKE, double myPmag)
{
    if (myKE >= 0.05)
    {
        return pow((pow(0.05, 2) + pow((0.3 / (pow(myE, 0.5))), 2)), 0.5);
    }
    else
    {
        return -1;
    }
}

double SigmaM_Res_duneCdr(double myE, double myKE, double myPmag)
{
    if (myKE >= 0.05)
    {
        return pow((pow(0.05, 2) + pow((0.3 / (pow(myE, 0.5))), 2)), 0.5);
    }
    else
    {
        return -1;
    }
}

// Mersenne twister pseudo-random number generation

std::random_device rd;
std::mt19937 gen(rd());

// Energy smearing function

double smearE(int myPdg, double myE, double myKE, double myPx, double myPy, double myPz, std::string model)
{
    std::map<int, int> myMap;

    myMap[kPdgPiP] = 0;
    myMap[kPdgPiM] = 1;
    myMap[kPdgPi0] = 2;
    myMap[kPdgKP] = 3;
    myMap[kPdgKM] = 4;
    myMap[kPdgGamma] = 5;
    myMap[kPdgProton] = 6;
    myMap[kPdgElectron] = 7;
    myMap[kPdgPositron] = 8;
    myMap[kPdgMuon] = 9;
    myMap[kPdgAntiMuon] = 10;
    myMap[kPdgNeutron] = 11;
    myMap[kPdgK0] = 12;
    myMap[kPdgAntiK0] = 13;
    myMap[kPdgLambda] = 14;
    myMap[kPdgSigmaP] = 15;
    myMap[kPdgSigma0] = 16;
    myMap[kPdgSigmaM] = 17;

    TVector3 P(myPx, myPy, myPz);
    TVector3 Z(0, 0, 1);
    double myPmag = P.Mag();

    int in = myMap.find(myPdg)->second;

    if (model == "duneCdr")
    {

        double (*resolution_ptr_duneCdr[18])(double, double, double) = {
            PiP_Res_duneCdr,
            PiM_Res_duneCdr,
            Pi0_Res_duneCdr,
            KP_Res_duneCdr,
            KM_Res_duneCdr,
            Gamma_Res_duneCdr,
            Proton_Res_duneCdr,
            Electron_Res_duneCdr,
            Positron_Res_duneCdr,
            Muon_Res_duneCdr,
            AntiMuon_Res_duneCdr,
            Neutron_Res_duneCdr,
            K0_Res_duneCdr,
            AntiK0_Res_duneCdr,
            Lambda_Res_duneCdr,
            SigmaP_Res_duneCdr,
            Sigma0_Res_duneCdr,
            SigmaM_Res_duneCdr};

        double resolution = (*resolution_ptr_duneCdr[in])(myE, myKE, myPmag);

        if (resolution == -1)
        {
            return 0;
        }
        else
        {
            double var = 0;
            double eSq = 0;

            if (pdg::IsNeutronOrProton(myPdg))
            {
                var += pow((resolution * myKE), 2);
                eSq += pow(myKE, 2);
            }
            else
            {
                var += pow((resolution * myE), 2);
                eSq += pow(myE, 2);
            }

            double m = log(eSq / (pow(var + eSq, 0.5)));
            double s = pow(log(1 + (var / eSq)), 0.5);

            std::lognormal_distribution<double> distLognorm(m, s);

            if (myPdg == kPdgNeutron)
            {
                if (myPmag < 1)
                {
                    std::uniform_real_distribution<double> distUni(0, 1);

                    if (distUni(gen) < 0.1)
                    {
                        return 0;
                    }
                    else
                    {
                        return 0.6 * (distLognorm(gen));
                    }
                }
                else
                {
                    return 0.6 * (distLognorm(gen));
                }
            }
            else
            {
                return distLognorm(gen);
            }
        }
    }

    // Default model resolutions and particle detection dependencies

    else if (model == "default")
    {
        double info[18][2] = {

            {0.15, 1},
            {0.15, 1},
            {0.15, 1},
            {0.2, 1},
            {0.2, 1},
            {0.3, 0.5},
            {0.4, 1},
            {0.4, 1},
            {0.4, 1},
            {0.15, 1},
            {0.15, 1},
            {0.5, 0.5},
            {0.2, 1},
            {0.2, 1},
            {0.3, 1},
            {0.3, 1},
            {0.3, 1},
            {0.3, 1}};

        double resolution = info[in][0];
        double chanceToSee = info[in][1];
        double var = 0;
        double eSq = 0;

        if (pdg::IsNeutronOrProton(myPdg))
        {
            var += pow((resolution * myKE), 2);
            eSq += pow(myKE, 2);
        }
        else
        {
            var += pow((resolution * myE), 2);
            eSq += pow(myE, 2);
        }

        double m = log(eSq / (pow(var + eSq, 0.5)));
        double s = pow(log(1 + (var / eSq)), 0.5);

        std::lognormal_distribution<double> distLognorm(m, s);
        std::uniform_real_distribution<double> distUni(0, 1);

        if (distUni(gen) < chanceToSee)
        {
            if (myKE > 0.05)
            {
                return distLognorm(gen);
            }
            else
            {
                return 0;
            }
        }
        else
        {
            return 0;
        }
    }
    else
    {
        std::cout << "Error: Resolution model not found./n";
    }
}

// Angular smearing function

double smearA(int myPdg, double myPx, double myPy, double myPz, std::string model)
{
    std::map<int, int> myMap;

    myMap[kPdgPiP] = 0;
    myMap[kPdgPiM] = 1;
    myMap[kPdgPi0] = 2;
    myMap[kPdgKP] = 3;
    myMap[kPdgKM] = 4;
    myMap[kPdgGamma] = 5;
    myMap[kPdgProton] = 6;
    myMap[kPdgElectron] = 7;
    myMap[kPdgPositron] = 8;
    myMap[kPdgMuon] = 9;
    myMap[kPdgAntiMuon] = 10;
    myMap[kPdgNeutron] = 11;
    myMap[kPdgK0] = 12;
    myMap[kPdgAntiK0] = 13;
    myMap[kPdgLambda] = 14;
    myMap[kPdgSigmaP] = 15;
    myMap[kPdgSigma0] = 16;
    myMap[kPdgSigmaM] = 17;

    TVector3 P(myPx, myPy, myPz);
    TVector3 Z(0, 0, 1);
    double theta = (P.Angle(Z)) / M_PI * 180;

    int in = myMap.find(myPdg)->second;

    // DUNE-CDR and Default angular smearing values

    double angularResDeg_DuneCdr[18] = {
        1,
        1,
        5,
        5,
        5,
        1,
        5,
        1,
        1,
        1,
        1,
        5,
        5,
        5,
        5,
        5,
        5,
        5};

    double angularResDeg_Default[18] = {
        2,
        2,
        8,
        2,
        2,
        3,
        8,
        2,
        2,
        2,
        2,
        10,
        8,
        8,
        8,
        8,
        8,
        8};

    double resolution = 0;

    if (model == "duneCdr")
    {
        resolution += (angularResDeg_DuneCdr[in]);
    }
    else if (model == "default")
    {
        resolution += (angularResDeg_Default[in]);
    }

    std::normal_distribution<double> distNorm(theta, resolution);
    return distNorm(gen);
}

#endif
