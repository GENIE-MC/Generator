cmake_minimum_required(VERSION 3.14)
project(Genie_XSection)

set(CMAKE_CXX_STANDARD)

add_library(Genie_XSection SHARED NCgamma_Diff_Cross_Section.cpp NCgamma_Diff_Cross_Section.h NCgamma_Form_Factors.h NCgamma_Form_Factors.cpp NCgamma_Hadronic_Current_R.h NCgamma_Hadronic_Current_R.cpp)

