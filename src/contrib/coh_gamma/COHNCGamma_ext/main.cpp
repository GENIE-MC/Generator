#include <iostream>
#include <fstream>
#include "NCgamma_Diff_Cross_Section.h"
#include "NCgamma_Parameters_GeV.h"

using namespace NC_gamma;

int main() {

    std::ofstream file;
    file.open("/home/edusaul/Desktop/d5_Delta_avg.dat");

    std::string nucleus = "12C";
    std::string mode = "nu";

    Diff_Cross_Section cs(mode, nucleus);

    double k0 =        1.0; //GeV
    double kg0 =       0.0;//0.3; //GeV
    double phig =      10.0 * 3.14159/180.0;//* NCgamma_Param::pi/180;
    double thetag =    10.0 * 3.14159/180.0;//* NCgamma_Param::pi/180;
    double th =        10.0 * 3.14159/180.0;//* NCgamma_Param::pi/180;
//    double Enu_final = 0.00525802;//k0 - kg0;

    int npts = 50;
    double Enu_final_min = 0.0;
    double Enu_final_max = k0;
    double Enu_final = Enu_final_max;

    double diff_cs;

//    std::cout<<"k0 = "<<k0<<" ,  Enu_final = "<<Enu_final<<" ,  theta = "<<th
//             <<" , theta_g  = "<<thetag<<" , phi_g = "<<phig
//             <<" ,  ds5 = "<<diff_cs<<"   "<<diff_cs/NCgamma_Param::hccm2<<std::endl;
//
//    for (int i = 0; i < 10; ++i) {
//        for (int j = 0; j < 10; ++j) {
//            double diff_cs =  cs.getDiffCrossSection(k0,Enu_final,th,thetag,phig);;
//
//            std::cout<<"k0 = "<<k0<<" ,  kg0 = "<<kg0<<" ,  theta = "<<th
//                     <<" , theta_g  = "<<thetag<<" , phi_g = "<<phig
//                     <<" ,  ds5 = "<<diff_cs<<"   "<<diff_cs/NCgamma_Param::hccm2<<std::endl;
//            thetag += 1.5/9;
//        }
//        th += 0.5/9;
//    }

    for (int i = 0; i < npts; ++i) {

        kg0 = k0 - Enu_final;

        diff_cs = cs.getDiffCrossSection(k0,Enu_final,th,thetag,phig);
//        diff_cs *= 197.33e-16 * 197.33e-16;

        std::cout<<"k0 = "<<k0<<" ,  El = "<<Enu_final<<" ,  kg0 = "<<kg0<<" ,  theta = "<<th
                     <<" , theta_g  = "<<thetag<<" , phi_g = "<<phig
                     <<" ,  ds5 = "<<diff_cs<<"   ";//<<diff_cs<<std::endl;

        diff_cs *= 197.33e-16 * 197.33e-16;
        std::cout<<"   "<<diff_cs<<std::endl;

        file<<kg0<<"  "<<diff_cs<<std::endl;

        Enu_final -= (Enu_final_max-Enu_final_min)/double(npts-1);

    }

    file.close();



    return 0;
}