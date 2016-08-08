#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ namespace genie;

#pragma link C++ class genie::InitialState;
#pragma link C++ class genie::Interaction;
#pragma link C++ class genie::Target;
#pragma link C++ class genie::ProcessInfo;
#pragma link C++ class genie::Kinematics+;
#pragma link C++ class genie::XclsTag;
#pragma link C++ class genie::KPhaseSpace;

#pragma link C++ class std::map<genie::KineVar_t,double>+; // in Kinematics object
#pragma link C++ class std::pair<genie::KineVar_t,double>+; // in Kinematics object

#pragma link C++ ioctortype TRootIOCtor;

#endif
