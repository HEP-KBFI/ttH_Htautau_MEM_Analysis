# include "ttH_Htautau_MEM_Analysis/MEMAlgo/interface/RunConfig.h"
# include "ttH_Htautau_MEM_Analysis/MEM/interface/PyConfig.h"

using namespace std;

RunConfig::RunConfig( const char *configname ) {
  maxNbrOfEventsToRead_ = LLONG_MAX;
  eventType_ = BadEventType;
  PyObject *pModule = readPyConfig( configname );
  if ( pModule == 0 ) eventType_ = BadEventType;
  PyObject *pVariable;
  PyObject *pItem;
  int i;
  string str;
  
  configFileName = configname;
  
  // Events
  pVariable = PyObject_GetAttrString( pModule, "maxNbrOfEventsToRead");
  maxNbrOfEventsToRead_ = PyInt_AsLong( pVariable );
  if  (maxNbrOfEventsToRead_ < 0) 
    maxNbrOfEventsToRead_ = LLONG_MAX;
 
  pVariable = PyObject_GetAttrString( pModule, "skipEvents");
  if (pVariable != NULL )
    skipEvents_ = PyInt_AsLong( pVariable );
  else
    skipEvents_ = 0;
  // Input : Integral value file  
  pVariable = PyObject_GetAttrString( pModule, "treeName");
  if (pVariable != NULL )
    treeName_ = PyString_AsString( pVariable );
  pVariable = PyObject_GetAttrString( pModule,"InputFileList");
  if ( pVariable && (PyList_Size(pVariable) >0) ) {
    int nb_items = PyList_Size(pVariable);
    for( i = 0; i < nb_items; ++i) {
      pItem = PyList_GetItem(pVariable, i);
      str = PyString_AsString( pItem );
      if (str.length() > 0)
        inputFileNames_.push_back( str );
    }
  } else 
    cerr << "PyConfig: Bad input filename list" << endl;
  
  //
  // Input type
  pVariable = PyObject_GetAttrString( pModule, "InputType");
  str = PyString_AsString( pVariable );
  //
  if ( str == "PyRun2") {
       eventType_ = PyRun2_EventType;
  }
  else {
    cerr << "PyConfig: Bad input filename type" << endl; 
    eventType_ = BadEventType;
  }
  
  
  
  // MadGraph
  MGParamCardFileName_ = "";
  pVariable = PyObject_GetAttrString( pModule, "MGParamCardFile");
  str = PyString_AsString( pVariable );
  //
  if (str.length() != 0) {
    MGParamCardFileName_ = str;
  }
  
   // LHAPDF
  LHAPDFFileName_ = "";
  pVariable = PyObject_GetAttrString( pModule, "LHAPDFFile");
  str = PyString_AsString( pVariable );
  //
  if (str.length() != 0) {
    LHAPDFFileName_ = str;
  } 

// Run ttZ integration
  runttZ_integration_ = false;
  pVariable = PyObject_GetAttrString( pModule, "runttZ_integration");
  runttZ_integration_ = (pVariable==Py_True);

  // Run ttZ_Zll integration
  runttZ_Zll_integration_ = false;
  pVariable = PyObject_GetAttrString( pModule, "runttZ_Zll_integration");
  runttZ_Zll_integration_ = (pVariable==Py_True);

  // Run ttbar_DL_fakelep integration
  runttbar_DL_fakelep_integration_ = false;
  pVariable = PyObject_GetAttrString( pModule, "runttbar_DL_fakelep_integration");
  runttbar_DL_fakelep_integration_ = (pVariable==Py_True);

  // Run missing jet integration
  run_missing_jet_integration_ = false;
  pVariable = PyObject_GetAttrString( pModule, "run_missing_jet_integration");
  run_missing_jet_integration_ = (pVariable==Py_True);

  // Force missing jet integration
  force_missing_jet_integration_ = false;
  pVariable = PyObject_GetAttrString( pModule, "force_missing_jet_integration");
  force_missing_jet_integration_ = (pVariable==Py_True);

  // Force missing jet integration if no MEM sel perm passing cuts without missing jet
  force_missing_jet_integration_ifnoperm_ = false;
  pVariable = PyObject_GetAttrString( pModule, "force_missing_jet_integration_ifnoperm");
  force_missing_jet_integration_ifnoperm_ = (pVariable==Py_True);

  // Random generator
  flagSameRNG_ = false;
  pVariable = PyObject_GetAttrString( pModule, "flagSameRNG");
  flagSameRNG_ = (pVariable==Py_True); 
  
  // Computation / run
  flagTFLepTau_ = true; flagTFHadTau_ = true;
  flagTFMET_ = true; flagTFJet1_ = true; flagTFJet2_ = true;
  flagTFBJet_leptop_ = true;   flagTFBJet_hadtop_ = true;  flagTFTop_ = true;
  flagTF_fakelep_ = true; flagTF_fakeleptau_ = true;
  flagWME_ = true; flagJac_ = true;
  //
  pVariable = PyObject_GetAttrString( pModule, "flagTFLepTau");
  flagTFLepTau_ = (pVariable==Py_True); 
  pVariable = PyObject_GetAttrString( pModule, "flagTFHadTau");
  flagTFHadTau_ = (pVariable==Py_True);   
  pVariable = PyObject_GetAttrString( pModule, "flagTFMET");
  flagTFMET_    = (pVariable==Py_True);  
  pVariable = PyObject_GetAttrString( pModule, "flagTFJet1");
  flagTFJet1_   = (pVariable==Py_True);
  pVariable = PyObject_GetAttrString( pModule, "flagTFJet2");
  flagTFJet2_   = (pVariable==Py_True); 
  pVariable = PyObject_GetAttrString( pModule, "flagTFBJet_leptop");
  flagTFBJet_leptop_   = (pVariable==Py_True);
  pVariable = PyObject_GetAttrString( pModule, "flagTFBJet_hadtop");
  flagTFBJet_hadtop_   = (pVariable==Py_True);
  pVariable = PyObject_GetAttrString( pModule, "flagTF_fakelep");
  flagTF_fakelep_   = (pVariable==Py_True);
  pVariable = PyObject_GetAttrString( pModule, "flagTF_fakeleptau");
  flagTF_fakeleptau_   = (pVariable==Py_True);
  pVariable = PyObject_GetAttrString( pModule, "flagTFTop");
  flagTFTop_   = (pVariable==Py_True);
  pVariable = PyObject_GetAttrString( pModule, "flagJac");
  flagJac_      = (pVariable==Py_True);  
  pVariable = PyObject_GetAttrString( pModule, "flagWME");
  flagWME_      = (pVariable==Py_True);

  // MET TF
  include_hadrecoil_ = true;
  pVariable = PyObject_GetAttrString( pModule, "include_hadrecoil");
  include_hadrecoil_  = PyInt_AsLong( pVariable ); 

  // Verbose mode  
  pVariable = PyObject_GetAttrString( pModule, "verbose");
  verbose_  = PyInt_AsLong( pVariable ); 
  pVariable = PyObject_GetAttrString( pModule, "logFileName");
  if(pVariable!=0)
    logFileName_ = PyString_AsString( pVariable );

  // MC integration
  nbrOfDim_ttH_ = 0;
  nbrOfDim_ttZ_ = 0;
  nbrOfDim_ttZ_Zll_ = 0;
  nbrOfDim_ttbar_DL_ = 0;

  nbrOfDim_ttH_miss_ = 0;
  nbrOfDim_ttZ_miss_ = 0;
  nbrOfDim_ttZ_Zll_miss_ = 0;

  pVariable = PyObject_GetAttrString( pModule, "nbrOfDimttH");
  nbrOfDim_ttH_ = PyInt_AsLong( pVariable ); 
  pVariable    = PyObject_GetAttrString( pModule, "nbrOfPoints_ttH");
  nbrOfPoints_ttH_ = PyLong_AsLong( pVariable );

  if(runttZ_integration_){
    pVariable = PyObject_GetAttrString( pModule, "nbrOfDimttZ");
    nbrOfDim_ttZ_ = PyInt_AsLong( pVariable );
    pVariable    = PyObject_GetAttrString( pModule, "nbrOfPoints_ttZ");
    nbrOfPoints_ttZ_ = PyLong_AsLong( pVariable );
  }
  if(runttbar_DL_fakelep_integration_){
    pVariable = PyObject_GetAttrString( pModule, "nbrOfDimttbar_DL");
    nbrOfDim_ttbar_DL_ = PyInt_AsLong( pVariable );
    pVariable    = PyObject_GetAttrString( pModule, "nbrOfPoints_ttbar_DL");
    nbrOfPoints_ttbar_DL_ = PyLong_AsLong( pVariable );
  }
  if(runttZ_Zll_integration_){
    pVariable = PyObject_GetAttrString( pModule, "nbrOfDimttZ_Zll");
    nbrOfDim_ttZ_Zll_ = PyInt_AsLong( pVariable );
    pVariable    = PyObject_GetAttrString( pModule, "nbrOfPoints_ttZ_Zll");
    nbrOfPoints_ttZ_Zll_ = PyLong_AsLong( pVariable );
  }


  if(run_missing_jet_integration_){

    pVariable = PyObject_GetAttrString( pModule, "nbrOfDimttH_miss");
    nbrOfDim_ttH_miss_ = PyInt_AsLong( pVariable ); 
    pVariable    = PyObject_GetAttrString( pModule, "nbrOfPoints_ttH_miss");
    nbrOfPoints_ttH_miss_ = PyLong_AsLong( pVariable ); 
 
    if(runttZ_integration_){
      pVariable = PyObject_GetAttrString( pModule, "nbrOfDimttZ_miss");
      nbrOfDim_ttZ_miss_ = PyInt_AsLong( pVariable );
      pVariable    = PyObject_GetAttrString( pModule, "nbrOfPoints_ttZ_miss");
      nbrOfPoints_ttZ_miss_ = PyLong_AsLong( pVariable );
    }
    if(runttZ_Zll_integration_){
      pVariable = PyObject_GetAttrString( pModule, "nbrOfDimttZ_Zll_miss");
      nbrOfDim_ttZ_Zll_miss_ = PyInt_AsLong( pVariable );
      pVariable    = PyObject_GetAttrString( pModule, "nbrOfPoints_ttZ_Zll_miss");
      nbrOfPoints_ttZ_Zll_miss_ = PyLong_AsLong( pVariable );
    }

  }
  

  pVariable    = PyObject_GetAttrString( pModule, "nbrOfPermut_per_jet");
  nbrOfPermut_per_jet_ = PyLong_AsLong( pVariable );


  use_pT_TFJet_ = false;
  pVariable = PyObject_GetAttrString( pModule, "use_pT_TFJet");
  use_pT_TFJet_ = (pVariable==Py_True);

  use_top_compatibility_check_ = false;
  pVariable = PyObject_GetAttrString( pModule, "use_top_compatibility_check");
  use_top_compatibility_check_ = (pVariable==Py_True);
  

  //Integration boundaries
  pVariable = PyObject_GetAttrString( pModule, "CI_TFJet");
  CI_TFJet_ = PyFloat_AsDouble( pVariable );

  pVariable = PyObject_GetAttrString( pModule, "eta_acceptance");
  if (pVariable != NULL )
    eta_acceptance_ = PyFloat_AsDouble( pVariable );
  pVariable = PyObject_GetAttrString( pModule, "jet_radius");
  if (pVariable != NULL )
    jet_radius_ = PyFloat_AsDouble( pVariable );
  pVariable = PyObject_GetAttrString( pModule, "dR_veto_jet_lep");
  if (pVariable != NULL )
    dR_veto_jet_lep_ = PyFloat_AsDouble( pVariable );
  pVariable = PyObject_GetAttrString( pModule, "rel_iso_lep");
  if (pVariable != NULL )
    rel_iso_lep_ = PyFloat_AsDouble( pVariable );
  pVariable = PyObject_GetAttrString( pModule, "pT_cut");
  if (pVariable != NULL )
    pT_cut_ = PyFloat_AsDouble( pVariable );

  if(run_missing_jet_integration_){

    pVariable    = PyObject_GetAttrString( pModule, "phi_missing_jet");
    if (pVariable != NULL ){
      pItem = PyList_GetItem(pVariable, 0);
      phi_missing_jet_Boundaries_[0] =  PyFloat_AsDouble( pItem );
      pItem = PyList_GetItem(pVariable, 1);
      phi_missing_jet_Boundaries_[1] =  PyFloat_AsDouble( pItem );
    }
    
    pVariable    = PyObject_GetAttrString( pModule, "cosTheta_missing_jet");
    if (pVariable != NULL ){
      pItem = PyList_GetItem(pVariable, 0);
      cosTheta_missing_jet_Boundaries_[0] =  PyFloat_AsDouble( pItem );
      pItem = PyList_GetItem(pVariable, 1);
      cosTheta_missing_jet_Boundaries_[1] =  PyFloat_AsDouble( pItem );
    }

  }


  pVariable    = PyObject_GetAttrString( pModule, "phiNu_tlep");
  pItem = PyList_GetItem(pVariable, 0);
  phiNu_tlep_Boundaries_[0] =  PyFloat_AsDouble( pItem );
  pItem = PyList_GetItem(pVariable, 1);
  phiNu_tlep_Boundaries_[1] =  PyFloat_AsDouble( pItem );

  pVariable    = PyObject_GetAttrString( pModule, "cosThetaNu_tlep");
  pItem = PyList_GetItem(pVariable, 0);
  cosThetaNu_tlep_Boundaries_[0] =  PyFloat_AsDouble( pItem );
  pItem = PyList_GetItem(pVariable, 1);
  cosThetaNu_tlep_Boundaries_[1] =  PyFloat_AsDouble( pItem );



  if( runttbar_DL_fakelep_integration_ ){
    
    pVariable    = PyObject_GetAttrString( pModule, "phiNu_ttau");
    pItem = PyList_GetItem(pVariable, 0);
    phiNu_ttau_ttbar_DL_Boundaries_[0] =  PyFloat_AsDouble( pItem );
    pItem = PyList_GetItem(pVariable, 1);
    phiNu_ttau_ttbar_DL_Boundaries_[1] =  PyFloat_AsDouble( pItem );
    
    pVariable    = PyObject_GetAttrString( pModule, "cosThetaNu_ttau");
    pItem = PyList_GetItem(pVariable, 0);
    cosThetaNu_ttau_ttbar_DL_Boundaries_[0] =  PyFloat_AsDouble( pItem );
    pItem = PyList_GetItem(pVariable, 1);
    cosThetaNu_ttau_ttbar_DL_Boundaries_[1] =  PyFloat_AsDouble( pItem );
 
  }

  sqrtS_ = 13000;
  pVariable = PyObject_GetAttrString( pModule, "sqrtS");
  if (pVariable != NULL )
    sqrtS_    = PyFloat_AsDouble( pVariable );   

  
}
