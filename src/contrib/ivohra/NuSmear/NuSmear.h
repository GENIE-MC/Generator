//____________________________________________________________________________
/*
\name     NuSmear Smearing System
\brief    Provides generic parameterized smearing presets for energy smearing 
          and angular smearing with two smearing models - the Default model 
          and the DUNE CDR model.
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



//duneCdr resolution functions

double PiP_Res_duneCdr(double myE, double myKE, double myPmag)
{
    if (myKE >= 0.1){ //100 MeV Threshold for PiP/M
            return 0.15;

    } else {
        return -1;
    }

}

double PiM_Res_duneCdr(double myE, double myKE, double myPmag)
{
    if (myKE >= 0.1){ //100 MeV Threshold for PiP/M
            return 0.15;

    } else {
        return -1;
    }
}

double Pi0_Res_duneCdr(double myE, double myKE, double myPmag)
{
    if (myKE >= 0.05){ //50 MeV Threshold for other

        return pow((pow(0.05,2) + pow((0.3/(pow(myE,0.5))),2)),0.5); // 5% (sum in quad) 30%/√(E)
    } else {
        return -1;
    }
}

double KP_Res_duneCdr(double myE, double myKE, double myPmag)
{
    if (myKE >= 0.05){ //50 MeV Threshold for other

        return pow((pow(0.05,2) + pow((0.3/(pow(myE,0.5))),2)),0.5); // 5% (sum in quad) 30%/√(E)
    } else {
        return -1;
    }
}

double KM_Res_duneCdr(double myE, double myKE, double myPmag)
{
    if (myKE >= 0.05){ //50 MeV Threshold for other

        return pow(pow(0.05,2) + pow((0.3/(pow(myE,0.5))),2),0.5); // 5% (sum in quad) 30%/√(E)
        
    } else {
        return -1;
    }
}

double Gamma_Res_duneCdr(double myE, double myKE, double myPmag)
{
    if (myKE >= 0.03){ //30 MeV Threshold for gamma

        return pow(pow(0.02,2) + pow(0.15/(pow(myE,0.5)),2), 0.5); // 2% (sum in quad) 15%/√(E)

    } else {

        return -1;
    }
}

double Proton_Res_duneCdr(double myE, double myKE, double myPmag)//Uses KE instead of E
{
    if (myKE >= 0.05){ //50 MeV Threshold for proton

        if (myPmag < 0.4){ //if momentum is < 400 MeV
            return 0.1;

        } else { //if momentum is ≥ 400 MeV
            return pow(pow(0.05,2) + pow(0.3/(pow(myKE,0.5)),2),0.5); // 5% (sum in quad) 30%/√(E)
        }

    } else {
        return -1;
    }
}


double Electron_Res_duneCdr(double myE, double myKE, double myPmag)
{
    if (myKE >= 0.03){ //30 MeV Threshold for electron

        return pow(pow(0.02,2) + pow(0.15/(pow(myE,0.5)),2), 0.5); // 2% (sum in quad) 15%/√(E)

    } else {

        return -1;
    }
}

double Positron_Res_duneCdr(double myE, double myKE, double myPmag)
{
        if (myKE >= 0.03){ //30 MeV Threshold for positron

        return pow(pow(0.02,2) + pow(0.15/(pow(myE,0.5)),2), 0.5); // 2% (sum in quad) 15%/√(E)

    } else {

        return -1;
    }
}

double Muon_Res_duneCdr(double myE, double myKE, double myPmag)
{
        if (myKE >= 0.03){ //30 MeV Threshold for Muon

            return 0.15;

    } else {
        return -1;
    }
}

double AntiMuon_Res_duneCdr(double myE, double myKE, double myPmag)
{
        if (myKE >= 0.03){ //30 MeV Threshold for Muon

            return 0.15;

    } else {
        return -1;
    }
}

double Neutron_Res_duneCdr(double myE, double myKE, double myPmag) //Uses KE instead of E
{
        if (myKE >= 0.05){ //50 MeV Threshold for n

        return 0.4/(pow(myKE,0.5)); // 40%/√(E)
    } else {
        return -1;
    }
}

double K0_Res_duneCdr(double myE, double myKE, double myPmag)
{
if (myKE >= 0.05){ //50 MeV Threshold for other

        return pow((pow(0.05,2) + pow((0.3/(pow(myE,0.5))),2)),0.5); // 5% (sum in quad) 30%/√(E)
    } else {
        return -1;
    }
}

double AntiK0_Res_duneCdr(double myE, double myKE, double myPmag)
{
if (myKE >= 0.05){ //50 MeV Threshold for other

        return pow((pow(0.05,2) + pow((0.3/(pow(myE,0.5))),2)),0.5); // 5% (sum in quad) 30%/√(E)
    } else {
        return -1;
    }
}

double Lambda_Res_duneCdr(double myE, double myKE, double myPmag)
{
if (myKE >= 0.05){ //50 MeV Threshold for other

        return pow((pow(0.05,2) + pow((0.3/(pow(myE,0.5))),2)),0.5); // 5% (sum in quad) 30%/√(E)
    } else {
        return -1;
    }
}

double SigmaP_Res_duneCdr(double myE, double myKE, double myPmag)
{
if (myKE >= 0.05){ //50 MeV Threshold for other

        return pow((pow(0.05,2) + pow((0.3/(pow(myE,0.5))),2)),0.5); // 5% (sum in quad) 30%/√(E)
    } else {
        return -1;
    }
}

double Sigma0_Res_duneCdr(double myE, double myKE, double myPmag)
{
if (myKE >= 0.05){ //50 MeV Threshold for other

        return pow((pow(0.05,2) + pow((0.3/(pow(myE,0.5))),2)),0.5); // 5% (sum in quad) 30%/√(E)
    } else {
        return -1;
    }
}

double SigmaM_Res_duneCdr(double myE, double myKE, double myPmag)
{
if (myKE >= 0.05){ //50 MeV Threshold for other

        return pow((pow(0.05,2) + pow((0.3/(pow(myE,0.5))),2)),0.5); // 5% (sum in quad) 30%/√(E)
    } else {
        return -1;
    }
}


//RNG

    std::random_device rd;
    std::mt19937 gen(rd()); //generate random numbers according to mersenne twister


//SmearE energy smearing function

double smearE(int myPdg, double myE, double myKE, double myPx, double myPy, double myPz, std::string model)
{

//map from PDG to myMap index number

std::map<int,int> myMap;

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

TVector3 P (myPx, myPy, myPz); //particle momentum vector

TVector3 Z (0,0,1); //vector pointing downstream along Z axis

double myPmag = P.Mag(); //particle momentum magnitude



int in = myMap.find(myPdg)->second; //in is the index number for the given particle

if (model == "duneCdr"){

//duneCdr resolution function pointers

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
SigmaM_Res_duneCdr
};


double resolution = (*resolution_ptr_duneCdr[in])(myE, myKE, myPmag); //calls duneCdr resolution function corresponding to particle with index number 'in', and stores resulting resolution value

if (resolution == -1){ //if the particle is below the threshold energy
    return 0; //i.e. particle was not detected at all

} else {

        double var = 0;
        double eSq = 0;

        if (pdg::IsNeutronOrProton(myPdg)){ //for nucleons, use KE instead of E

        var += pow((resolution*myKE),2); //define Var[X]
        eSq += pow(myKE,2); //define E[X]^2

        } else {

        var += pow((resolution*myE),2); //define Var[X]
        eSq += pow(myE,2); //define E[X]^2

        }

    double m = log(eSq/(pow(var+eSq,0.5)));
    double s = pow(log(1+(var/eSq)),0.5);

     std::lognormal_distribution<double> distLognorm(m,s); //takes the parameters: m and s

    if (myPdg == kPdgNeutron){

        if (myPmag < 1){ //if neutron momentum is less than 1 GeV/c

             std::uniform_real_distribution<double> distUni(0,1);
    
          if (distUni(gen) < 0.1){ //10% chance of escaping detection
            return 0;
          } else {

            return 0.6*(distLognorm(gen)); // neutron 60% reconstructed energy

          }
        } else {

            return 0.6*(distLognorm(gen)); // neutron 60% reconstructed energy
        }

    } else {

        return distLognorm(gen);

    }
    
}


} else if (model == "default"){

double info[18][2] = {

{0.15, 1}, //first element is resolution, second element is chance of being observed (chanceToSee)
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
{0.3, 1}
};

double resolution = info[in][0];
double chanceToSee = info[in][1];

double var = 0;
double eSq = 0;

    if (pdg::IsNeutronOrProton(myPdg)){ //for nucleons, use KE instead of E
          var += pow((resolution*myKE),2); //define Var[X]
          eSq += pow(myKE,2); //define E[X]^2

    } else {

        var += pow((resolution*myE),2); //define Var[X]
        eSq += pow(myE,2); //define E[X]^2
    }

    double m = log(eSq/(pow(var+eSq,0.5)));
    double s = pow(log(1+(var/eSq)),0.5);

            std::lognormal_distribution<double> distLognorm(m,s); //takes the parameters: m and the s
           
            std::uniform_real_distribution<double> distUni(0,1); //generates a random double between 0 and 1 (i.e. in the interval [0,1) ) according to uniform distribution
           
          if (distUni(gen) < chanceToSee){ //i.e. particle is observed


                     if(myKE > 0.05){ //Threshold of 50 MeV 
                        return distLognorm(gen);
                     } else {
                        return 0;
                     }

          } else {
            return 0;
          }


} else {
    std::cout << "Error: Resolution model not found./n";
}

}




double smearA(int myPdg, double myPx, double myPy, double myPz, std::string model)
{

//map from PDG to myMap index number

std::map<int,int> myMap;

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

TVector3 P (myPx, myPy, myPz); //particle momentum vector

TVector3 Z (0,0,1); //vector pointing downstream along Z axis

double theta = (P.Angle(Z))/M_PI*180; //angle away from Z vector in degrees

int in = myMap.find(myPdg)->second; //in is the index number for the given particle

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
5
};

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
8
};

double resolution = 0;

if (model == "duneCdr"){

resolution += (angularResDeg_DuneCdr[in]); //DUNE-CDR angular resolution


} else if (model == "default") {

resolution += (angularResDeg_Default[in]); //Default angular resolution 

}
     
std::normal_distribution<double> distNorm(theta,resolution); //gaussian distribution centered at true theta with standard deviation = resolution
return distNorm(gen);

}

#endif