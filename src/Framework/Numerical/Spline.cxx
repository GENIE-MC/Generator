//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab
*/
//____________________________________________________________________________

#include <cassert>
#include <iomanip>
#include <cfloat>

#include "libxml/parser.h"
#include "libxml/xmlmemory.h"

#include <TFile.h>
#include <TNtupleD.h>
#include <TTree.h>
#include <TSQLServer.h>
#include <TGraph.h>
#include <TMath.h>
#include <TH2F.h>
#include <TROOT.h>

#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/Spline.h"
#include "Framework/Numerical/MathUtils.h"
#include "Framework/Utils/XmlParserUtils.h"

using std::endl;
using std::setw;
using std::setprecision;
using std::setfill;
using std::ios;

using namespace genie;

ClassImp(Spline)

//___________________________________________________________________________
namespace genie
{
  ostream & operator << (ostream & stream, const Spline & spl)
  {
     spl.Print(stream);
     return stream;
  }
}
//___________________________________________________________________________
Spline::Spline()
{
  this->InitSpline();
}
//___________________________________________________________________________
Spline::Spline(string filename, string xtag, string ytag, bool is_xml) :
TObject()
{
  string fmt = (is_xml) ? "XML" : "ASCII";

  LOG("Spline", pDEBUG)
      << "Constructing spline from data in " << fmt << " file: " << filename;

  this->InitSpline();

  if(is_xml)
     this->LoadFromXmlFile(filename, xtag, ytag);
  else
     this->LoadFromAsciiFile(filename);
}
//___________________________________________________________________________
Spline::Spline(TNtupleD * ntuple, string var, string cut) :
TObject()
{
  LOG("Spline", pDEBUG) << "Constructing spline from data in a TNtuple";

  this->InitSpline();
  this->LoadFromNtuple(ntuple, var, cut);
}
//___________________________________________________________________________
Spline::Spline(TTree * tree, string var, string cut) :
TObject()
{
  LOG("Spline", pDEBUG) << "Constructing spline from data in a TTree";

  this->InitSpline();
  this->LoadFromTree(tree, var, cut);
}
//___________________________________________________________________________
Spline::Spline(TSQLServer * db, string query) :
TObject()
{
  LOG("Spline", pDEBUG) << "Constructing spline from data in a MySQL server";

  this->InitSpline();
  this->LoadFromDBase(db, query);
}
//___________________________________________________________________________
Spline::Spline(int nentries, double x[], double y[]) :
TObject()
{
  LOG("Spline", pDEBUG)
                 << "Constructing spline from the arrays passed to the ctor";

  this->InitSpline();
  this->BuildSpline(nentries, x, y);
}
//___________________________________________________________________________
Spline::Spline(int nentries, float x[], float y[]) :
TObject()
{
  LOG("Spline", pDEBUG)
                 << "Constructing spline from the arrays passed to the ctor";

  double * dblx = new double[nentries];
  double * dbly = new double[nentries];

  for(int i = 0; i < nentries; i++) {
     dblx[i] = double( x[i] );
     dbly[i] = double( y[i] );
  }

  this->InitSpline();
  this->BuildSpline(nentries, dblx, dbly);

  delete [] dblx;
  delete [] dbly;
}
//___________________________________________________________________________
Spline::Spline(const Spline & spline) :
  TObject(), fInterpolator(0)
{
  LOG("Spline", pDEBUG) << "Spline copy constructor";

  this->LoadFromTSpline3( *spline.GetAsTSpline(), spline.NKnots() );
}
//___________________________________________________________________________
Spline::Spline(const TSpline3 & spline, int nknots) :
  TObject(), fInterpolator(0)
{
  LOG("Spline", pDEBUG)
                    << "Constructing spline from the input TSpline3 object";

  this->LoadFromTSpline3( spline, nknots );
}
//___________________________________________________________________________
Spline::~Spline()
{
  if(fInterpolator) delete fInterpolator;
}
//___________________________________________________________________________
bool Spline::LoadFromXmlFile(string filename, string xtag, string ytag)
{
  LOG("Spline", pDEBUG) << "Retrieving data from file: " << filename;

  xmlDocPtr xml_doc = xmlParseFile(filename.c_str());

  if(xml_doc==NULL) {
    LOG("Spline", pERROR)
           << "XML file could not be parsed! [filename: " << filename << "]";
    return false;
  }

  xmlNodePtr xmlCur = xmlDocGetRootElement(xml_doc);

  if(xmlCur==NULL) {
    LOG("Spline", pERROR)
         << "XML doc. has null root element! [filename: " << filename << "]";
    return false;
  }
  if( xmlStrcmp(xmlCur->name, (const xmlChar *) "spline") ) {
    LOG("Spline", pERROR)
      << "XML doc. has invalid root element! [filename: " << filename << "]";
    return false;
  }

  string name = utils::str::TrimSpaces(
                    utils::xml::GetAttribute(xmlCur, "name"));
  string snkn = utils::str::TrimSpaces(
                    utils::xml::GetAttribute(xmlCur, "nknots"));
  int nknots = atoi( snkn.c_str() );

  LOG("Spline", pINFO)
             << "Parsing XML spline: " << name << ", nknots = " << nknots;

  int      iknot = 0;
  double * vx    = new double[nknots];
  double * vy    = new double[nknots];

  xmlNodePtr xmlSplChild = xmlCur->xmlChildrenNode; // <knots>'s

  // loop over all xml tree nodes that are children of the <spline> node
  while (xmlSplChild != NULL) {
     LOG("Spline", pDEBUG)
                 << "Got <spline> children node: " << xmlSplChild->name;

     // enter everytime you find a <knot> tag
     if( (!xmlStrcmp(xmlSplChild->name, (const xmlChar *) "knot")) ) {

        xmlNodePtr xmlKnotChild = xmlSplChild->xmlChildrenNode;

        // loop over all xml tree nodes that are children of this <knot>
        while (xmlKnotChild != NULL) {
           LOG("Spline", pDEBUG)
                 << "Got <knot> children node: "  << xmlKnotChild->name;

           // enter everytime you find a <E> or a <xsec> tag
           const xmlChar * tag = xmlKnotChild->name;
           bool is_xtag = ! xmlStrcmp(tag,(const xmlChar *) xtag.c_str());
           bool is_ytag = ! xmlStrcmp(tag,(const xmlChar *) ytag.c_str());
           if (is_xtag || is_ytag) {
              xmlNodePtr xmlValTagChild = xmlKnotChild->xmlChildrenNode;
              string val = utils::xml::TrimSpaces(
                         xmlNodeListGetString(xml_doc, xmlValTagChild, 1));

              if (is_xtag) vx[iknot] = atof(val.c_str());
              if (is_ytag) vy[iknot] = atof(val.c_str());

              xmlFree(xmlValTagChild);
              LOG("Spline", pDEBUG) << "tag: " << tag << ", value: " << val;
           }//if current tag is <E>,<xsec>

           xmlKnotChild = xmlKnotChild->next;
        }//[end of] loop over tags within <knot>...</knot> tags
        xmlFree(xmlKnotChild);
        iknot++;
     } // if <knot>
     xmlSplChild = xmlSplChild->next;
  } //[end of] loop over tags within <spline>...</spline> tags
  xmlFree(xmlSplChild);

  this->BuildSpline(nknots, vx, vy);
  delete [] vx;
  delete [] vy;

  return true;
}
//___________________________________________________________________________
bool Spline::LoadFromAsciiFile(string filename)
{
  LOG("Spline", pDEBUG) << "Retrieving data from file: " << filename;

  TNtupleD nt("ntuple","","x:y");
  nt.ReadFile(filename.c_str());

  this->LoadFromNtuple(&nt, "x:y");
  return true;
}
//___________________________________________________________________________
bool Spline::LoadFromNtuple(TNtupleD * nt, string var, string cut)
{
  if(!nt) return false;

  TTree * tree = dynamic_cast<TTree *> (nt);
  return this->LoadFromTree(tree,var,cut);
}
//___________________________________________________________________________
bool Spline::LoadFromTree(TTree * tree, string var, string cut)
{
  LOG("Spline", pDEBUG) << "Retrieving data from tree: " << tree->GetName();

  if(!cut.size()) tree->Draw(var.c_str(), "",          "GOFF");
  else            tree->Draw(var.c_str(), cut.c_str(), "GOFF");

  TH2F * hst = (TH2F*)gROOT->FindObject("htemp");
  if(hst) { hst->SetDirectory(0); delete hst; }

  // Now, take into account that the data retrieved from the ntuple would
  // not be sorted in x and the resulting spline will be bogus...
  // Sort the x array, use the x-sorting to re-arrange the y array and only
  // then build the spline..

  int n = tree->GetSelectedRows();

  int *    idx = new int[n];
  double * x   = new double[n];
  double * y   = new double[n];

  TMath::Sort(n,tree->GetV1(),idx,false);

  for(int i=0; i<n; i++) {
     int ii = idx[i];
     x[i]   = (tree->GetV1())[ii];
     y[i]   = (tree->GetV2())[ii];
  }

  this->BuildSpline(n,x,y);

  delete [] idx;
  delete [] x;
  delete [] y;

  return true;
}
//___________________________________________________________________________
bool Spline::LoadFromDBase(TSQLServer * /*db*/,  string /*query*/)
{
  LOG("Spline", pDEBUG) << "Retrieving data from data-base: ";
  return false;
}
//___________________________________________________________________________
bool Spline::LoadFromTSpline3(const TSpline3 & spline, int nknots)
{
  int nentries = nknots;

  double * x = new double[nentries];
  double * y = new double[nentries];
  double cx = 0, cy = 0;

  for(int i = 0; i < nentries; i++) {
     spline.GetKnot(i, cx, cy);
     x[i] = cx;
     y[i] = cy;
  }
  this->BuildSpline(nentries, x, y);

  delete [] x;
  delete [] y;

  return true;
}
//___________________________________________________________________________
void Spline::GetKnot(int iknot, double & x, double & y) const
{
  if(!fInterpolator) {
     LOG("Spline", pWARN) << "Spline has not been built yet!";
     return;
  }
  fInterpolator->GetKnot(iknot,x,y);
}
//___________________________________________________________________________
double Spline::GetKnotX(int iknot) const
{
  if(!fInterpolator) {
     LOG("Spline", pWARN) << "Spline has not been built yet!";
     return 0;
  }
  double x,y;
  fInterpolator->GetKnot(iknot,x,y);
  return x;
}
//___________________________________________________________________________
double Spline::GetKnotY(int iknot) const
{
  if(!fInterpolator) {
     LOG("Spline", pWARN) << "Spline has not been built yet!";
     return 0;
  }
  double x,y;
  fInterpolator->GetKnot(iknot,x,y);
  return y;
}
//___________________________________________________________________________
bool Spline::IsWithinValidRange(double x) const
{
  bool is_in_range = (fXMin <= x && x <= fXMax);
  return is_in_range;
}
//___________________________________________________________________________
double Spline::Evaluate(double x) const
{
  LOG("Spline", pDEBUG) << "Evaluating spline at point x = " << x;
  assert(!TMath::IsNaN(x));

  double y = 0;
  if( this->IsWithinValidRange(x) ) {

    // we can interpolate within the range of spline knots - be careful with
    // strange cubic spline behaviour when close to knots with y=0
    bool is0p = this->ClosestKnotValueIsZero(x, "+");
    bool is0n = this->ClosestKnotValueIsZero(x, "-");

    if(!is0p && !is0n) {
      // both knots (on the left and right are non-zero) - just interpolate
      LOG("Spline", pDEBUG) << "Point is between non-zero knots";
      y = fInterpolator->Eval(x);
    } else {
      // at least one of the neighboring knots has y=0
      if(is0p && is0n) {
        // both neighboring knots have y=0
        LOG("Spline", pDEBUG) << "Point is between zero knots";
        y=0;
      } else {
        // just 1 neighboring knot has y=0 - do a linear interpolation
        LOG("Spline", pDEBUG) 
          << "Point has zero" << (is0n ? " left " : " right ") << "knot";
        double xpknot=0, ypknot=0, xnknot=0, ynknot=0;
        this->FindClosestKnot(x, xnknot, ynknot, "-");
        this->FindClosestKnot(x, xpknot, ypknot, "+");
        if(is0n) y = ypknot * (x-xnknot)/(xpknot-xnknot);
        else     y = ynknot * (x-xnknot)/(xpknot-xnknot);
      }
    }

  } else {
    LOG("Spline", pDEBUG) << "x = " << x
     << " is not within spline range [" << fXMin << ", " << fXMax << "]";
  }

  if(y<0 && !fYCanBeNegative) {
    LOG("Spline", pINFO) << "Negative y (" << y << ")";
    LOG("Spline", pINFO) << "x = " << x;
    LOG("Spline", pINFO) << "spline range [" << fXMin << ", " << fXMax << "]";
  }

  LOG("Spline", pDEBUG) << "Spline(x = " << x << ") = " << y;

  return y;
}
//___________________________________________________________________________
void Spline::SaveAsXml(
                string filename, string xtag, string ytag, string name) const
{
  ofstream outxml(filename.c_str());
  if(!outxml.is_open()) {
    LOG("Spline", pERROR) << "Couldn't create file = " << filename;
    return;
  }
  string spline_name = (name.size()>0 ? name : fName);

  this->SaveAsXml(outxml, xtag, ytag, spline_name);

  outxml << endl;
  outxml.close();
}
//___________________________________________________________________________
void Spline::SaveAsXml(
    ofstream & ofs, string xtag, string ytag, string name) const
{
  string spline_name = (name.size()>0 ? name : fName);

  // create a spline tag with the number of knots as an attribute
  int nknots = this->NKnots();
  string padding = "    ";
  ofs << padding << "<spline name=\"" << spline_name
      << "\" nknots=\"" << nknots << "\">" << endl;

  // start printing the knots
  double x=0, y=0;
  for(int iknot = 0; iknot < nknots; iknot++)
  {
    fInterpolator->GetKnot(iknot, x, y);

    ofs  << std::fixed << setprecision(5);
    ofs  << "\t<knot>"
         << " <" << xtag << "> " << setfill(' ')
                                 << setw(10) << x << " </" << xtag << ">";
    ofs  << std::scientific << setprecision(10);
    ofs  << " <" << ytag << "> " << setfill(' ')
                                 << setw(10) << y << " </" << ytag << ">"
         << " </knot>" << endl;
  }
  ofs << padding << "</spline>" << endl;
}
//___________________________________________________________________________
void Spline::SaveAsText(string filename, string format) const
{
  ofstream outtxt(filename.c_str());
  if(!outtxt.is_open()) {
    LOG("Spline", pERROR) << "Couldn't create file = " << filename;
    return;
  }
  int nknots = this->NKnots();
  outtxt << nknots << endl;

  double x=0, y=0;
  for(int iknot = 0; iknot < nknots; iknot++) {
    fInterpolator->GetKnot(iknot, x, y);
    char line[1024];
    sprintf(line,format.c_str(),x,y);
    outtxt << line << endl;
  }
  outtxt << endl;
  outtxt.close();
}
//___________________________________________________________________________
void Spline::SaveAsROOT(string filename, string name, bool recreate) const
{
  string spline_name = (name.size()>0 ? name : fName);

  string opt = ( (recreate) ? "RECREATE" : "UPDATE" );

  TFile f(filename.c_str(), opt.c_str());
  if(fInterpolator) fInterpolator->Write(spline_name.c_str());
  f.Close();
}
//___________________________________________________________________________
TGraph * Spline::GetAsTGraph(
          int np, bool scale_with_x, bool in_log, double fx, double fy) const
{
  double xmin = fXMin;
  double xmax = fXMax;

  np = TMath::Max(np,2);

  bool use_log = in_log && (fXMin>0) && (fXMax>0);

  if(use_log) {
    xmin = TMath::Log10(xmin);
    xmax = TMath::Log10(xmax);
  }

  double dx = (xmax - xmin) / (np-1);

  double * x = new double[np];
  double * y = new double[np];

  for(int i=0; i<np; i++) {
      x[i] = ( (use_log) ? TMath::Power(10, xmin+i*dx) : xmin + i*dx );
      y[i] = this->Evaluate( x[i] );

      // scale with x if needed
      if (scale_with_x) y[i] /= x[i];

      // apply x,y scaling
      y[i] *= fy;
      x[i] *= fx;
  }

  TGraph * graph = new TGraph(np, x, y);
  delete[] x;
  delete[] y;
  return graph;
}
//___________________________________________________________________________
void Spline::FindClosestKnot(
              double x, double & xknot, double & yknot, Option_t * opt) const
{
  string option(opt);

  bool pos = (option.find("+") != string::npos);
  bool neg = (option.find("-") != string::npos);

  if(!pos && !neg) return;

  int iknot = fInterpolator->FindX(x);

  double xp=0, yp=0, xn=0, yn=0;
  fInterpolator->GetKnot(iknot,  xn,yn);
  fInterpolator->GetKnot(iknot+1,xp,yp);

  bool p = (TMath::Abs(x-xp) < TMath::Abs(x-xn));

  if(pos&&neg) {
    if(p) { xknot = xp; yknot = yp; }
    else  { xknot = xn; yknot = yn; }
  } else {
    if(pos) { xknot = xp; yknot = yp; }
    if(neg) { xknot = xn; yknot = yn; }
  }
}
//___________________________________________________________________________
bool Spline::ClosestKnotValueIsZero(double x, Option_t * opt) const
{
  double xknot = 0, yknot = 0;
  this->FindClosestKnot(x, xknot, yknot, opt);
  if(utils::math::AreEqual(yknot,0)) return true;
  return false;
}
//___________________________________________________________________________
void Spline::Print(ostream & stream) const
{
  int    nknots = this->NKnots();
  double xmin   = this->XMin();
  double xmax   = this->XMax();     

  stream << endl;
  stream << "** Spline: " << this->Name() << endl;
  stream << "Has " << nknots 
        << " knots in the [" << xmin << ", " << xmax << "] range" << endl;
  double x,y;
  for(int i=0; i<nknots; i++) {  
    this->GetKnot(i,x,y); 
    stream << "-- knot : " << i 
           << " -> (x = " << x << ", y = " << y << ")" << endl; 
  }
}
//___________________________________________________________________________
void Spline::Add(const Spline & spl, double c)
{
  // continue only if the input spline is defined at a wider x-range
  double xmin = this->XMin();
  double xmax = this->XMax();     
  bool ok = spl.IsWithinValidRange(xmin) && spl.IsWithinValidRange(xmax);
  if(!ok) {
    LOG("Spline", pERROR) << "** Can't add splines. X-range mismatch!";
    return;
  }

  int nknots = this->NKnots();
  double * x = new double[nknots];
  double * y = new double[nknots];

  for(int i=0; i<nknots; i++) {  
    this->GetKnot(i,x[i],y[i]); 
    y[i] += (c * spl.Evaluate(x[i]));
  }
  this->ResetSpline();
  this->BuildSpline(nknots,x,y);
  delete [] x;
  delete [] y;
}
//___________________________________________________________________________
void Spline::Multiply(const Spline & spl, double c)
{
  // continue only if the input spline is defined at a wider x-range
  double xmin = this->XMin();
  double xmax = this->XMax();     
  bool ok = spl.IsWithinValidRange(xmin) && spl.IsWithinValidRange(xmax);
  if(!ok) {
    LOG("Spline", pERROR) << "** Can't multiply splines. X-range mismatch!";
    return;
  }

  int nknots = this->NKnots();
  double * x = new double[nknots];
  double * y = new double[nknots];

  for(int i=0; i<nknots; i++) {  
    this->GetKnot(i,x[i],y[i]); 
    y[i] *= (c * spl.Evaluate(x[i]));
  }
  this->ResetSpline();
  this->BuildSpline(nknots,x,y);
  delete [] x;
  delete [] y;
}
//___________________________________________________________________________
void Spline::Divide(const Spline & spl, double c)
{
  // continue only if the input spline is defined at a wider x-range
  double xmin = this->XMin();
  double xmax = this->XMax();     
  bool ok = spl.IsWithinValidRange(xmin) && spl.IsWithinValidRange(xmax);
  if(!ok) {
    LOG("Spline", pERROR) << "** Can't divide splines. X-range mismatch!";
    return;
  }

  int nknots = this->NKnots();
  double * x = new double[nknots];
  double * y = new double[nknots];

  for(int i=0; i<nknots; i++) {  
    this->GetKnot(i,x[i],y[i]); 
    double denom = c * spl.Evaluate(x[i]);
    bool denom_is_zero = TMath::Abs(denom) < DBL_EPSILON;
    if(denom_is_zero) {
        LOG("Spline", pERROR) << "** Refusing to divide spline knot by 0";
        delete [] x;
        delete [] y;
        return;
    }
    y[i] /= denom;
  }
  this->ResetSpline();
  this->BuildSpline(nknots,x,y);
  delete [] x;
  delete [] y;
}
//___________________________________________________________________________
void Spline::Add(double a)
{
  int nknots = this->NKnots();
  double * x = new double[nknots];
  double * y = new double[nknots];

  for(int i=0; i<nknots; i++) {  
    this->GetKnot(i,x[i],y[i]); 
    y[i]+=a;
  }
  this->ResetSpline();
  this->BuildSpline(nknots,x,y);
  delete [] x;
  delete [] y;
}
//___________________________________________________________________________
void Spline::Multiply(double a)
{
  int nknots = this->NKnots();
  double * x = new double[nknots];
  double * y = new double[nknots];

  for(int i=0; i<nknots; i++) {  
    this->GetKnot(i,x[i],y[i]); 
    y[i]*=a;
  }
  this->ResetSpline();
  this->BuildSpline(nknots,x,y);
  delete [] x;
  delete [] y;
}
//___________________________________________________________________________
void Spline::Divide(double a)
{
  bool a_is_zero = TMath::Abs(a) < DBL_EPSILON;
  if(a_is_zero==0) {
    LOG("Spline", pERROR) << "** Refusing to divide spline by 0";
    return;
  }
  int nknots = this->NKnots();
  double * x = new double[nknots];
  double * y = new double[nknots];

  for(int i=0; i<nknots; i++) {  
    this->GetKnot(i,x[i],y[i]); 
    y[i]/=a;
  }
  this->ResetSpline();
  this->BuildSpline(nknots,x,y);
  delete [] x;
  delete [] y;
}
//___________________________________________________________________________
void Spline::InitSpline(void)
{
  LOG("Spline", pDEBUG) << "Initializing spline...";

  fName = "genie-spline";
  fXMin = 0.0;
  fXMax = 0.0;
  fYMax = 0.0;

  fInterpolator = 0;

  fYCanBeNegative = false;

  LOG("Spline", pDEBUG) << "...done initializing spline";
}
//___________________________________________________________________________
void Spline::ResetSpline(void)
{
  if(fInterpolator) delete fInterpolator;
  this->InitSpline();
}
//___________________________________________________________________________
void Spline::BuildSpline(int nentries, double x[], double y[])
{
  LOG("Spline", pDEBUG) << "Building spline...";

  double xmin = x[0];          // minimum x in spline
  double xmax = x[nentries-1]; // maximum x in spline

  fNKnots = nentries;
  fXMin   = xmin;
  fXMax   = xmax;
  fYMax   = y[ TMath::LocMax(nentries, y) ]; // maximum y in spline

  if(fInterpolator) delete fInterpolator;

  fInterpolator = new TSpline3("spl3", x, y, nentries, "0");

  LOG("Spline", pDEBUG) << "...done building spline";
}
//___________________________________________________________________________
