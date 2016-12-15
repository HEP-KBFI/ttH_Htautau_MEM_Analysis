# include "ttH_Htautau_MEM_Analysis/MEMAlgo/interface/Processes.h"

# include "ttH_Htautau_MEM_Analysis/MEMAlgo/interface/PyConfig.h"

# include "LHAPDF/LHAPDF.h"

#include "ttH_Htautau_MEM_Analysis/MEMAlgo/interface/CPPProcess_ttH_gg.h"
#include "ttH_Htautau_MEM_Analysis/MEMAlgo/interface/CPPProcess_ttH_qqbar.h"
#include "ttH_Htautau_MEM_Analysis/MEMAlgo/interface/CPPProcess_ttZ_Zonly_gg.h"
#include "ttH_Htautau_MEM_Analysis/MEMAlgo/interface/CPPProcess_ttZ_Zonly_uubar.h"
#include "ttH_Htautau_MEM_Analysis/MEMAlgo/interface/CPPProcess_ttZ_Zonly_Zll_gg.h"
#include "ttH_Htautau_MEM_Analysis/MEMAlgo/interface/CPPProcess_ttZ_Zonly_Zll_uubar.h"
#include "ttH_Htautau_MEM_Analysis/MEMAlgo/interface/CPPProcess_gg_to_ttbar.h"
#include "ttH_Htautau_MEM_Analysis/MEMAlgo/interface/CPPProcess_uubar_to_ttbar.h"

#include "ttH_Htautau_MEM_Analysis/MEMAlgo/interface/Constants.h"



// Main ttH processes
static CPPProcess_ttH_gg process_ttH_gg;
static CPPProcess_ttH_qqbar process_ttH_qqbar;

// Main ttZ processes
static CPPProcess_ttZ_Zonly_gg process_ttZ_Zonly_gg;
static CPPProcess_ttZ_Zonly_uubar process_ttZ_Zonly_uubar;

static CPPProcess_ttZ_Zonly_Zll_gg process_ttZ_Zonly_Zll_gg;
static CPPProcess_ttZ_Zonly_Zll_uubar process_ttZ_Zonly_Zll_uubar;

// Main ttbar processes
static CPPProcess_gg_to_ttbar process_gg_ttbar;
static CPPProcess_uubar_to_ttbar process_qqbar_ttbar;


// Init with python configuration file
void init_ttH_HTauTauMEProcesses( const char *configname ) {

  PyObject *pModule = readPyConfig( configname );
  if ( pModule == 0 ) {
    cerr << "Module " << configname << " not found" << endl;
    cerr << "  or bad python configuration." << endl;
    return;
  }
  
  PyObject *pVariable;
  string paramFileName;
  
  pVariable = PyObject_GetAttrString( pModule, "MGParamCardFile");
  paramFileName = PyString_AsString( pVariable ); 
  
  // ttH
  process_ttH_gg.initProc(paramFileName);
  process_ttH_qqbar.initProc(paramFileName);

  // ttZ
  process_ttZ_Zonly_gg.initProc(paramFileName);
  process_ttZ_Zonly_uubar.initProc(paramFileName);

  process_ttZ_Zonly_Zll_gg.initProc(paramFileName);
  process_ttZ_Zonly_Zll_uubar.initProc(paramFileName);

  // ttbar
  process_gg_ttbar.initProc(paramFileName);
  process_qqbar_ttbar.initProc(paramFileName);

  Py_XDECREF( pVariable );
  closePyConfig( pModule );
}



// ttH
double getME_ttH_gg( vector<double*> &p ) {
  process_ttH_gg.setMomenta( p );
  process_ttH_gg.sigmaKin();
  return
  process_ttH_gg.getMatrixElements()[0];
}

double getME_ttH_qqbar( vector<double*> &p ) {
  process_ttH_qqbar.setMomenta( p );
  process_ttH_qqbar.sigmaKin();
  return
  process_ttH_qqbar.getMatrixElements()[0];
}

// ttZ
double getME_ttZ_Zonly_gg( vector<double*> &p ) {
  process_ttZ_Zonly_gg.setMomenta( p );
  process_ttZ_Zonly_gg.sigmaKin();
  return
  process_ttZ_Zonly_gg.getMatrixElements()[0];
}

double getME_ttZ_Zonly_uubar( vector<double*> &p ) {
  process_ttZ_Zonly_uubar.setMomenta( p );
  process_ttZ_Zonly_uubar.sigmaKin();
  return
  process_ttZ_Zonly_uubar.getMatrixElements()[0];
}

double getME_ttZ_Zonly_Zll_gg( vector<double*> &p ) {
  process_ttZ_Zonly_Zll_gg.setMomenta( p );
  process_ttZ_Zonly_Zll_gg.sigmaKin();
  return
  process_ttZ_Zonly_Zll_gg.getMatrixElements()[0];
}

double getME_ttZ_Zonly_Zll_uubar( vector<double*> &p ) {
  process_ttZ_Zonly_Zll_uubar.setMomenta( p );
  process_ttZ_Zonly_Zll_uubar.sigmaKin();
  return
  process_ttZ_Zonly_Zll_uubar.getMatrixElements()[0];
}


// ttbar
double getME_gg_ttbar( vector<double*> &p ) {
  process_gg_ttbar.setMomenta( p );
  process_gg_ttbar.sigmaKin();
  return
  process_gg_ttbar.getMatrixElements()[0];
}

double getME_qqbar_ttbar( vector<double*> &p ) {
  process_qqbar_ttbar.setMomenta( p );
  process_qqbar_ttbar.sigmaKin();
  return
  process_qqbar_ttbar.getMatrixElements()[0];
}




double get_ttHWeightedME( double x0, double x1, const vector<double*> &p) {
  
  double pdfQ_ = Physics::mtop + 0.5*Physics::mHiggs; 
  double current_ME2, pdf_sum;
  double ME2 = 0;
  
  // The 2 permutations
  vector<double*> p0=p;
  vector<double*> p1=p;  
  p1[0]=p[1];
  p1[1]=p[0];

  
  // Parton codes:
  //    t  b  c  s  u  d  g  d  u  s  c  b  t
  //   -6 -5 -4 -3 -2 -1  0  1  2  3  4  5  6
  double fg0    = LHAPDF::xfx( PDFSet, x0, pdfQ_,  0 );
  double fu0    = LHAPDF::xfx( PDFSet, x0, pdfQ_,  2 );
  double fd0    = LHAPDF::xfx( PDFSet, x0, pdfQ_,  1 );
  double fc0    = LHAPDF::xfx( PDFSet, x0, pdfQ_,  4 );
  double fs0    = LHAPDF::xfx( PDFSet, x0, pdfQ_,  3 );
  double fubar0 = LHAPDF::xfx( PDFSet, x0, pdfQ_, -2 );
  double fdbar0 = LHAPDF::xfx( PDFSet, x0, pdfQ_, -1 );
  double fcbar0 = LHAPDF::xfx( PDFSet, x0, pdfQ_, -4 );
  double fsbar0 = LHAPDF::xfx( PDFSet, x0, pdfQ_, -3 );

  double fg1    = LHAPDF::xfx( PDFSet, x1, pdfQ_,  0 );
  double fu1    = LHAPDF::xfx( PDFSet, x1, pdfQ_,  2 );
  double fd1    = LHAPDF::xfx( PDFSet, x1, pdfQ_,  1 );
  double fc1    = LHAPDF::xfx( PDFSet, x1, pdfQ_,  4 );
  double fs1    = LHAPDF::xfx( PDFSet, x1, pdfQ_,  3 );
  double fubar1 = LHAPDF::xfx( PDFSet, x1, pdfQ_, -2 );
  double fdbar1 = LHAPDF::xfx( PDFSet, x1, pdfQ_, -1 );
  double fcbar1 = LHAPDF::xfx( PDFSet, x1, pdfQ_, -4 );
  double fsbar1 = LHAPDF::xfx( PDFSet, x1, pdfQ_, -3 );
  
  //gg->ttH
  current_ME2 = getME_ttH_gg( p0 );
  pdf_sum     = fg0*fg1;
  ME2        += pdf_sum * current_ME2;

  //qqbar->ttH: 2 perms give same ME
  /*current_ME2 = getME_ttH_qqbar( p0 );
  pdf_sum     = fu0*fubar1+fd0*fdbar1+fc0*fcbar1+fs0*fsbar1;
  ME2        += pdf_sum * current_ME2;
  //
  current_ME2 = getME_ttH_qqbar( p1 );
  pdf_sum     = fu1*fubar0+fd1*fdbar0+fc1*fcbar0+fs1*fsbar0;
  ME2        += pdf_sum * current_ME2;*/

  current_ME2 = getME_ttH_qqbar( p0 );
  pdf_sum     = fu0*fubar1+fd0*fdbar1+fc0*fcbar1+fs0*fsbar1;
  pdf_sum    += fu1*fubar0+fd1*fdbar0+fc1*fcbar0+fs1*fsbar0;
  ME2        += pdf_sum * current_ME2;

  return ME2 / (x0 * x1);

}




double get_ttZZonlyWeightedME( double x0, double x1, const vector<double*> &p) {
  
  double pdfQ_ = Physics::mtop + 0.5*Physics::mZ; 
  double current_ME2, pdf_sum;
  double ME2 = 0;
  
  // The 2 permutations
  vector<double*> p0=p;
  vector<double*> p1=p;  
  p1[0]=p[1];
  p1[1]=p[0];

  
  // Parton codes:
  //    t  b  c  s  u  d  g  d  u  s  c  b  t
  //   -6 -5 -4 -3 -2 -1  0  1  2  3  4  5  6
  double fg0    = LHAPDF::xfx( PDFSet, x0, pdfQ_,  0 );
  double fu0    = LHAPDF::xfx( PDFSet, x0, pdfQ_,  2 );
  double fd0    = LHAPDF::xfx( PDFSet, x0, pdfQ_,  1 );
  double fc0    = LHAPDF::xfx( PDFSet, x0, pdfQ_,  4 );
  double fs0    = LHAPDF::xfx( PDFSet, x0, pdfQ_,  3 );
  double fubar0 = LHAPDF::xfx( PDFSet, x0, pdfQ_, -2 );
  double fdbar0 = LHAPDF::xfx( PDFSet, x0, pdfQ_, -1 );
  double fcbar0 = LHAPDF::xfx( PDFSet, x0, pdfQ_, -4 );
  double fsbar0 = LHAPDF::xfx( PDFSet, x0, pdfQ_, -3 );

  double fg1    = LHAPDF::xfx( PDFSet, x1, pdfQ_,  0 );
  double fu1    = LHAPDF::xfx( PDFSet, x1, pdfQ_,  2 );
  double fd1    = LHAPDF::xfx( PDFSet, x1, pdfQ_,  1 );
  double fc1    = LHAPDF::xfx( PDFSet, x1, pdfQ_,  4 );
  double fs1    = LHAPDF::xfx( PDFSet, x1, pdfQ_,  3 );
  double fubar1 = LHAPDF::xfx( PDFSet, x1, pdfQ_, -2 );
  double fdbar1 = LHAPDF::xfx( PDFSet, x1, pdfQ_, -1 );
  double fcbar1 = LHAPDF::xfx( PDFSet, x1, pdfQ_, -4 );
  double fsbar1 = LHAPDF::xfx( PDFSet, x1, pdfQ_, -3 );

  //gg->ttZ
  current_ME2 = getME_ttZ_Zonly_gg( p0 );
  pdf_sum     = fg0*fg1;
  ME2        += pdf_sum * current_ME2;

  //qqbar->ttZ
  current_ME2 = getME_ttZ_Zonly_uubar( p0 );
  pdf_sum     = fu0*fubar1+fc0*fcbar1+fd0*fdbar1+fs0*fsbar1;
  ME2        += pdf_sum * current_ME2;
  //
  current_ME2 = getME_ttZ_Zonly_uubar( p1 );
  pdf_sum     = fu1*fubar0+fc1*fcbar0+fd1*fdbar0+fs1*fsbar0;
  ME2        += pdf_sum * current_ME2; 

  return ME2 / (x0 * x1);

}






double get_ttZZonlyZllWeightedME( double x0, double x1, const vector<double*> &p) {
  
  double pdfQ_ = Physics::mtop + 0.5*Physics::mZ; 
  double current_ME2, pdf_sum;
  double ME2 = 0;
  
  // The 2 permutations
  vector<double*> p0=p;
  vector<double*> p1=p;  
  p1[0]=p[1];
  p1[1]=p[0];

  
  // Parton codes:
  //    t  b  c  s  u  d  g  d  u  s  c  b  t
  //   -6 -5 -4 -3 -2 -1  0  1  2  3  4  5  6
  double fg0    = LHAPDF::xfx( PDFSet, x0, pdfQ_,  0 );
  double fu0    = LHAPDF::xfx( PDFSet, x0, pdfQ_,  2 );
  double fd0    = LHAPDF::xfx( PDFSet, x0, pdfQ_,  1 );
  double fc0    = LHAPDF::xfx( PDFSet, x0, pdfQ_,  4 );
  double fs0    = LHAPDF::xfx( PDFSet, x0, pdfQ_,  3 );
  double fubar0 = LHAPDF::xfx( PDFSet, x0, pdfQ_, -2 );
  double fdbar0 = LHAPDF::xfx( PDFSet, x0, pdfQ_, -1 );
  double fcbar0 = LHAPDF::xfx( PDFSet, x0, pdfQ_, -4 );
  double fsbar0 = LHAPDF::xfx( PDFSet, x0, pdfQ_, -3 );

  double fg1    = LHAPDF::xfx( PDFSet, x1, pdfQ_,  0 );
  double fu1    = LHAPDF::xfx( PDFSet, x1, pdfQ_,  2 );
  double fd1    = LHAPDF::xfx( PDFSet, x1, pdfQ_,  1 );
  double fc1    = LHAPDF::xfx( PDFSet, x1, pdfQ_,  4 );
  double fs1    = LHAPDF::xfx( PDFSet, x1, pdfQ_,  3 );
  double fubar1 = LHAPDF::xfx( PDFSet, x1, pdfQ_, -2 );
  double fdbar1 = LHAPDF::xfx( PDFSet, x1, pdfQ_, -1 );
  double fcbar1 = LHAPDF::xfx( PDFSet, x1, pdfQ_, -4 );
  double fsbar1 = LHAPDF::xfx( PDFSet, x1, pdfQ_, -3 );

  //gg->ttZ
  current_ME2 = getME_ttZ_Zonly_Zll_gg( p0 );
  pdf_sum     = fg0*fg1;
  ME2        += pdf_sum * current_ME2;

  //qqbar->ttZ
  current_ME2 = getME_ttZ_Zonly_Zll_uubar( p0 );
  pdf_sum     = fu0*fubar1+fc0*fcbar1+fd0*fdbar1+fs0*fsbar1;
  ME2        += pdf_sum * current_ME2;
  //
  current_ME2 = getME_ttZ_Zonly_Zll_uubar( p1 );
  pdf_sum     = fu1*fubar0+fc1*fcbar0+fd1*fdbar0+fs1*fsbar0;
  ME2        += pdf_sum * current_ME2; 

  return ME2 / (x0 * x1);

}





double get_ttbarWeightedME( double x0, double x1, const vector<double*> &p) {
  
  double pdfQ_ = Physics::mtop; 
  double current_ME2, pdf_sum;
  double ME2 = 0;
  
  // The 2 permutations
  vector<double*> p0=p;
  vector<double*> p1=p;  
  p1[0]=p[1];
  p1[1]=p[0];

  
  // Parton codes:
  //    t  b  c  s  u  d  g  d  u  s  c  b  t
  //   -6 -5 -4 -3 -2 -1  0  1  2  3  4  5  6
  double fg0    = LHAPDF::xfx( PDFSet, x0, pdfQ_,  0 );
  double fu0    = LHAPDF::xfx( PDFSet, x0, pdfQ_,  2 );
  double fd0    = LHAPDF::xfx( PDFSet, x0, pdfQ_,  1 );
  double fc0    = LHAPDF::xfx( PDFSet, x0, pdfQ_,  4 );
  double fs0    = LHAPDF::xfx( PDFSet, x0, pdfQ_,  3 );
  double fubar0 = LHAPDF::xfx( PDFSet, x0, pdfQ_, -2 );
  double fdbar0 = LHAPDF::xfx( PDFSet, x0, pdfQ_, -1 );
  double fcbar0 = LHAPDF::xfx( PDFSet, x0, pdfQ_, -4 );
  double fsbar0 = LHAPDF::xfx( PDFSet, x0, pdfQ_, -3 );

  double fg1    = LHAPDF::xfx( PDFSet, x1, pdfQ_,  0 );
  double fu1    = LHAPDF::xfx( PDFSet, x1, pdfQ_,  2 );
  double fd1    = LHAPDF::xfx( PDFSet, x1, pdfQ_,  1 );
  double fc1    = LHAPDF::xfx( PDFSet, x1, pdfQ_,  4 );
  double fs1    = LHAPDF::xfx( PDFSet, x1, pdfQ_,  3 );
  double fubar1 = LHAPDF::xfx( PDFSet, x1, pdfQ_, -2 );
  double fdbar1 = LHAPDF::xfx( PDFSet, x1, pdfQ_, -1 );
  double fcbar1 = LHAPDF::xfx( PDFSet, x1, pdfQ_, -4 );
  double fsbar1 = LHAPDF::xfx( PDFSet, x1, pdfQ_, -3 );

  //gg->ttbar
  current_ME2 = getME_gg_ttbar( p0 );
  pdf_sum     = fg0*fg1;
  ME2        += pdf_sum * current_ME2;


  //qqbar->ttbar: 2 perms give same ME
  /*current_ME2 = getME_qqbar_ttbar( p0 );
  pdf_sum     = fu0*fubar1+fd0*fdbar1+fc0*fcbar1+fs0*fsbar1;
  ME2        += pdf_sum * current_ME2;
  //
  current_ME2 = getME_qqbar_ttbar( p1 );
  pdf_sum     = fu1*fubar0+fd1*fdbar0+fc1*fcbar0+fs1*fsbar0;
  ME2        += pdf_sum * current_ME2;*/

  current_ME2 = getME_qqbar_ttbar( p0 );
  pdf_sum     = fu0*fubar1+fd0*fdbar1+fc0*fcbar1+fs0*fsbar1;
  pdf_sum    += fu1*fubar0+fd1*fdbar0+fc1*fcbar0+fs1*fsbar0;
  ME2        += pdf_sum * current_ME2;
 
  return ME2 / (x0 * x1);

}
