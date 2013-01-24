#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ namespace genie;

#pragma link C++ class genie::RandomGen;
#pragma link C++ class genie::Spline;
#pragma link C++ class genie::BLI2DGrid;
#pragma link C++ class genie::BLI2DUnifGrid;
#pragma link C++ class genie::BLI2DNonUnifGrid;

//
// to be replaced with GSL/MathMore equivalents
//
#pragma link C++ class genie::GridSpacing;
#pragma link C++ class genie::GridDimension;
#pragma link C++ class genie::UnifGridDimension;
#pragma link C++ class genie::UnifGrid;
#pragma link C++ class genie::GBFunc;
#pragma link C++ class genie::GSFunc;
#pragma link C++ class genie::GVFunc;
#pragma link C++ class genie::FunctionMap;
#pragma link C++ class genie::IntegratorI;
#pragma link C++ class genie::Simpson1D;
#pragma link C++ class genie::Simpson2D;
#pragma link C++ class genie::Simpson2DFixN;

#endif
