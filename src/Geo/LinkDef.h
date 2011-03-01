#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ namespace genie;
#pragma link C++ namespace genie::geometry;

#pragma link C++ class genie::geometry::ROOTGeomAnalyzer;
#pragma link C++ class genie::geometry::PointGeomAnalyzer;

#pragma link C++ namespace genie::utils::geometry;

#pragma link C++ class genie::geometry::PathSegment;
#pragma link C++ class genie::geometry::PathSegmentList;
#pragma link C++ class genie::geometry::GeomVolSelectorI;
#pragma link C++ class genie::geometry::GeomVolSelectorBasic;

#pragma link C++ class genie::geometry::RayIntercept;
#pragma link C++ class genie::geometry::PlaneParam;
#pragma link C++ class genie::geometry::FidShape;
#pragma link C++ class genie::geometry::FidSphere;
#pragma link C++ class genie::geometry::FidCylinder;
#pragma link C++ class genie::geometry::FidPolyhedron;
#pragma link C++ class genie::geometry::GeomVolSelectorFiducial;

#pragma link C++ function genie::geometry::operator<<(ostream&, const genie::geometry::RayIntercept&);
#pragma link C++ function genie::geometry::operator<<(ostream&, const genie::geometry::PlaneParam&);
#pragma link C++ function genie::geometry::operator<<(ostream&, const genie::geometry::FidShape&);

#endif
