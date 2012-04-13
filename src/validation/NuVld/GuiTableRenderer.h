//_____________________________________________________________________________
/*!

\class    genie::nuvld::GuiTableRenderer

\brief    Responds to GUI events for graphical rendering of DBTable<T>'s

\author   Costas Andreopoulos (Rutherford Lab.)  <costas.andreopoulos \at stfc.ac.uk>

\created  August 26, 2004
*/
//_____________________________________________________________________________

#ifndef _GUI_TABLE_RENDERER_H_
#define _GUI_TABLE_RENDERER_H_

#include <string>
#include <TObject.h>
#include  "ValidationTools/NuVld/DBTable.h"

using std::string;

namespace genie {
namespace nuvld {

const double kRangeInLog = 0.8;

class DBNuXSecTableRow;
class DBElDiffXSecTableRow;
class DBSFTableRow;

class GuiTableRenderer : public TObject {

  friend class NuVldMainFrame;
  
private:
  GuiTableRenderer();
  GuiTableRenderer(TRootEmbeddedCanvas * emb_canvas);
  virtual ~GuiTableRenderer();

  //-- table drawing methods
    
  void DrawXSecTable (genie::nuvld::DBTable<DBNuXSecTableRow> *     table);
  void DrawXSecTable (genie::nuvld::DBTable<DBElDiffXSecTableRow> * table);
  void DrawXSecTable (genie::nuvld::DBTable<DBSFTableRow> *         table);

  //-- table drawing options
  
  void SwitchEmbdeddedCanvas (TRootEmbeddedCanvas * emb_canvas);
  void SetScaleWithEnergy    (bool tf);
  void SetMultigraph         (bool tf);
  void SetCustomXRange       (double xmin, double xmax);
  void ResetCustomXRange     (void);
  void SetPlotVariable       (string plot_var);
  void SetDrawOption         (string draw_option);
  void SetErrorOption        (string err_option);
  void SetExternalLegend     (TLegend * legend);

  void PrintDrawingOptions   (void);
  
  //-- actual drawing methods

  void DrawGraphInCanvas      (genie::nuvld::DBTable<DBNuXSecTableRow> *     table);
  void DrawMultiGraphInCanvas (genie::nuvld::DBTable<DBNuXSecTableRow> *     table);    
  void DrawGraphInCanvas      (genie::nuvld::DBTable<DBElDiffXSecTableRow> * table);
  void DrawMultiGraphInCanvas (genie::nuvld::DBTable<DBElDiffXSecTableRow> * table);
  void DrawGraphInCanvas      (genie::nuvld::DBTable<DBSFTableRow> *         table);
  void DrawMultiGraphInCanvas (genie::nuvld::DBTable<DBSFTableRow> *         table);

  //-- data members

  bool   fScaleWithE;
  bool   fIsMultigraph;
  bool   fHaveCustomXRange;
  double fXmin;        
  double fXmax;
  string fPlotVar;  
  string fDrawOpt;
  string fErrorOpt;
    
  TRootEmbeddedCanvas * fEmbCanvas;
  TLegend *             fLegend;

ClassDef(GuiTableRenderer, 0)
};

} // nuvld namespace
} // genie namespace

#endif // _GUI_TABLE_RENDERER_H_

