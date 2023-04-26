#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ namespace genie;
#pragma link C++ namespace genie::utils::hnl;

#pragma link C++ namespace genie::hnl;
#pragma link C++ namespace genie::hnl::enums;
#pragma link C++ namespace genie::hnl::selector;

#pragma link C++ class genie::hnl::ChannelCalculatorI;
#pragma link C++ class genie::hnl::GeomRecordVisitorI;
#pragma link C++ class genie::hnl::FluxRecordVisitorI;
#pragma link C++ class genie::hnl::DecayRecordVisitorI;

#pragma link C++ class genie::hnl::SimpleHNL;
#pragma link C++ class genie::hnl::BRCalculator;
#pragma link C++ class genie::hnl::FluxContainer;
#pragma link C++ class genie::hnl::FluxCreator;
#pragma link C++ class genie::hnl::Decayer;
#pragma link C++ class genie::hnl::VertexGenerator;

#pragma link C++ class genie::DummyHNLInteractionListGenerator;
#pragma link C++ class genie::DummyHNLPXSec;

#endif
