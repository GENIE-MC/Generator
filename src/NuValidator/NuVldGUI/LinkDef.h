#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ nestedclass;
#pragma link C++ nestedtypedef;

#pragma link C++ namespace genie;
#pragma link C++ namespace genie::nuvld;

#pragma link C++ class genie::nuvld::MsgBox;
#pragma link C++ class genie::nuvld::HelpBox;
#pragma link C++ class genie::nuvld::MultiLineMsgBox;
#pragma link C++ class genie::nuvld::YNQuestionBox;
#pragma link C++ class genie::nuvld::DBConnectionDialog;
#pragma link C++ class genie::nuvld::DBConnection;
#pragma link C++ class genie::nuvld::TextEntryDialog;
#pragma link C++ class genie::nuvld::NeuGenConfigDialog;
#pragma link C++ class genie::nuvld::NeuGenInputDialog;
#pragma link C++ class genie::nuvld::NeuGenFitParamsDialog;
#pragma link C++ class genie::nuvld::SysLogSingleton;
#pragma link C++ class genie::nuvld::BrowserSingleton;
#pragma link C++ class genie::nuvld::DataSelectionDialog;
#pragma link C++ class genie::nuvld::vDataSelectionDialog;
#pragma link C++ class genie::nuvld::vMeasurementListDialog;
#pragma link C++ class genie::nuvld::vDataSelectionTab;
#pragma link C++ class genie::nuvld::eDataSelectionTab;
#pragma link C++ class genie::nuvld::SFDataSelectionTab;
#pragma link C++ class genie::nuvld::NuVldMainFrame;
#pragma link C++ class genie::nuvld::NuVldUserData;
#pragma link C++ class genie::nuvld::NeuGenCards;
#pragma link C++ class genie::nuvld::GuiStackHandler;
#pragma link C++ class genie::nuvld::GuiHelpHandler;
#pragma link C++ class genie::nuvld::GuiDBHandler;
#pragma link C++ class genie::nuvld::GuiFitKernel;
#pragma link C++ class genie::nuvld::GuiXmlFileHandler;
#pragma link C++ class genie::nuvld::GuiTablePrinter;
#pragma link C++ class genie::nuvld::GuiTableRenderer;

#pragma link C++ function xsec_fitfunc;
#pragma link C++ function xsec_fitfunc_e;
#pragma link C++ function exclusive_xsec_fitfunc; 
#pragma link C++ function exclusive_xsec_fitfunc_e; 

#endif
