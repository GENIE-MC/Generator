#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ nestedclass;
#pragma link C++ nestedtypedef;

#pragma link C++ namespace genie::nuvld;

#pragma link C++ class genie::nuvld::DBI;
#pragma link C++ class genie::nuvld::DBQueryString;
#pragma link C++ class genie::nuvld::DBXmlUploader;
#pragma link C++ class genie::nuvld::DBTableRow;
#pragma link C++ class genie::nuvld::DBTableFields;
#pragma link C++ class genie::nuvld::SqlQueryBuilder;

#pragma link C++ class genie::nuvld::vXSecTableFields;
#pragma link C++ class genie::nuvld::vXSecTableRow;
#pragma link C++ class genie::nuvld::eDiffXSecTableFields;
#pragma link C++ class genie::nuvld::eDiffXSecTableRow;
#pragma link C++ class genie::nuvld::SFTableFields;
#pragma link C++ class genie::nuvld::SFTableRow;

#pragma link C++ class genie::nuvld::DBTable<genie::nuvld::vXSecTableRow>+;
#pragma link C++ class genie::nuvld::DBTable<genie::nuvld::eDiffXSecTableRow>+;
#pragma link C++ class genie::nuvld::DBTable<genie::nuvld::SFTableRow>+;
#pragma link C++ class genie::nuvld::DBTableStack<genie::nuvld::vXSecTableRow>+;
#pragma link C++ class genie::nuvld::DBTableStack<genie::nuvld::eDiffXSecTableRow>+;
#pragma link C++ class genie::nuvld::DBTableStack<genie::nuvld::SFTableRow>+;

#pragma link C++ class genie::nuvld::MultiGraph;
#pragma link C++ class genie::nuvld::MeasurementId;
#pragma link C++ class genie::nuvld::MeasurementIdList;

//--------rm
#pragma link C++ class genie::nuvld::SqlTag;


#endif
