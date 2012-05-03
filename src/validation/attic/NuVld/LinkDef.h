#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ nestedclass;
#pragma link C++ nestedtypedef;

#pragma link C++ namespace genie::nuvld;

#pragma link C++ class genie::nuvld::SqlQueryBuilder;
#pragma link C++ class genie::nuvld::MultiGraph;

#pragma link C++ class genie::nuvld::DBI;
#pragma link C++ class genie::nuvld::DBQueryString;
#pragma link C++ class genie::nuvld::DBXmlUploader;
#pragma link C++ class genie::nuvld::DBTableRow;
#pragma link C++ class genie::nuvld::DBTableFields;
#pragma link C++ class genie::nuvld::DBNuXSecTableFields;
#pragma link C++ class genie::nuvld::DBNuXSecTableRow;
#pragma link C++ class genie::nuvld::DBElDiffXSecTableFields;
#pragma link C++ class genie::nuvld::DBElDiffXSecTableRow;
#pragma link C++ class genie::nuvld::DBSFTableFields;
#pragma link C++ class genie::nuvld::DBSFTableRow;
#pragma link C++ class genie::nuvld::DBTable<genie::nuvld::DBNuXSecTableRow>+;
#pragma link C++ class genie::nuvld::DBTable<genie::nuvld::DBElDiffXSecTableRow>+;
#pragma link C++ class genie::nuvld::DBTable<genie::nuvld::DBSFTableRow>+;
#pragma link C++ class genie::nuvld::DBTableStack<genie::nuvld::DBNuXSecTableRow>+;
#pragma link C++ class genie::nuvld::DBTableStack<genie::nuvld::DBElDiffXSecTableRow>+;
#pragma link C++ class genie::nuvld::DBTableStack<genie::nuvld::DBSFTableRow>+;
#pragma link C++ class genie::nuvld::DBMeasurementId;
#pragma link C++ class genie::nuvld::DBMeasurementIdList;

#pragma link C++ class genie::nuvld::GuiMsgBox;
#pragma link C++ class genie::nuvld::GuiHelpBox;
#pragma link C++ class genie::nuvld::GuiMultiLineMsgBox;
#pragma link C++ class genie::nuvld::GuiYNQuestionBox;
#pragma link C++ class genie::nuvld::DBConnectionDialog;
#pragma link C++ class genie::nuvld::DBConnection;
#pragma link C++ class genie::nuvld::GuiTextEntryDialog;
#pragma link C++ class genie::nuvld::GuiSysLogSingleton;
#pragma link C++ class genie::nuvld::GuiBrowserSingleton;
#pragma link C++ class genie::nuvld::GuiDataSelectionDialog;
#pragma link C++ class genie::nuvld::GuiNuDataSelectionDialog;
#pragma link C++ class genie::nuvld::GuiNuMeasurementListDialog;
#pragma link C++ class genie::nuvld::GuiNuDataSelectionTab;
#pragma link C++ class genie::nuvld::GuiElDataSelectionTab;
#pragma link C++ class genie::nuvld::GuiSFDataSelectionTab;
#pragma link C++ class genie::nuvld::NuVldConfig;
#pragma link C++ class genie::nuvld::NuVldMainFrame;
#pragma link C++ class genie::nuvld::NuVldUserData;
#pragma link C++ class genie::nuvld::GuiStackHandler;
#pragma link C++ class genie::nuvld::GuiHelpHandler;
#pragma link C++ class genie::nuvld::GuiDBHandler;
#pragma link C++ class genie::nuvld::GuiXmlFileHandler;
#pragma link C++ class genie::nuvld::GuiTablePrinter;
#pragma link C++ class genie::nuvld::GuiTableRenderer;

#endif
