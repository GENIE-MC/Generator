#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ namespace genie;
#pragma link C++ namespace genie::flux;

#pragma link C++ class genie::flux::GJPARCNuFlux;
#pragma link C++ class genie::flux::GJPARCNuFluxPassThroughInfo;

#pragma link C++ class genie::flux::GNuMIFlux;
#pragma link C++ class genie::flux::GNuMIFluxPassThroughInfo;
#pragma link C++ class genie::flux::GNuMIFlux::StdFluxWindow_t;

#pragma link C++ class genie::flux::GCylindTH1Flux;
#pragma link C++ class genie::flux::GMonoEnergeticFlux;

#pragma link C++ class genie::flux::GAtmoFlux;
#pragma link C++ class genie::flux::GFlukaAtmo3DFlux;
#pragma link C++ class genie::flux::GBartolAtmoFlux;

#pragma link C++ class genie::flux::GAstroFlux;
#pragma link C++ class genie::flux::GPointSourceAstroFlux;
#pragma link C++ class genie::flux::GDiffuseAstroFlux;

#pragma link C++ class genie::flux::GSimpleNtpEntry+;
#pragma link C++ class genie::flux::GSimpleNtpNuMI+;
#pragma link C++ class genie::flux::GSimpleNtpAux+;
#pragma link C++ class genie::flux::GSimpleNtpMeta+;

#pragma link C++ class genie::flux::GSimpleNtpFlux;

#pragma link C++ class genie::flux::GFluxBlender;

#pragma link C++ class genie::flux::GFlavorMixerI;
#pragma link C++ class genie::flux::GFlavorMixerFactory;
#pragma link C++ class genie::flux::GFlavorMap;

#pragma link C++ class genie::flux::GFluxDriverFactory;

#pragma link C++ enum          genie::flux::EExposure;
#pragma link C++ nestedtypedef genie::flux::Exposure_t;
#pragma link C++ class genie::flux::GFluxExposureI;
#pragma link C++ class genie::flux::GFluxFileConfigI;

#endif
