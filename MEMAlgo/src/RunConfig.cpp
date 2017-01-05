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
  

  cout<<endl;
  
  // MadGraph
  MGParamCardFileName_ = "";
  pVariable = PyObject_GetAttrString( pModule, "MGParamCardFile");
  str = PyString_AsString( pVariable );
  //
  if (str.length() != 0) {
    MGParamCardFileName_ = str;
  }
  else{
    cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
    cout<<"!!WARNING: MGParamCardFile missing!!"<<endl;
    cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
  }

   // LHAPDF
  LHAPDFFileName_ = "";
  pVariable = PyObject_GetAttrString( pModule, "LHAPDFFile");
  str = PyString_AsString( pVariable );
  //
  if (str.length() != 0) {
    LHAPDFFileName_ = str;
  } 
  else{
    cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
    cout<<"!!WARNING: LHAPDFFile missing!!!"<<endl;
    cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
  }

  // Run ttZ integration
  runttZ_integration_ = false;
  pVariable = PyObject_GetAttrString( pModule, "runttZ_integration");
  if (pVariable != NULL )
    runttZ_integration_ = (pVariable==Py_True);
  else
    cout<<"WARNING: runttZ_integration missing, ttZ integration skipped by default"<<endl;


  // Run ttZ_Zll integration
  runttZ_Zll_integration_ = false;
  pVariable = PyObject_GetAttrString( pModule, "runttZ_Zll_integration");
  if (pVariable != NULL )
    runttZ_Zll_integration_ = (pVariable==Py_True);
  else
    cout<<"WARNING: runttZ_Zll_integration missing, ttZ Zll integration skipped by default"<<endl;

  // Run ttbar_DL_fakelep integration
  runttbar_DL_fakelep_integration_ = false;
  pVariable = PyObject_GetAttrString( pModule, "runttbar_DL_fakelep_integration");
  if (pVariable != NULL )
    runttbar_DL_fakelep_integration_ = (pVariable==Py_True);
  else
    cout<<"WARNING: runttbar_DL_fakelep_integration missing, ttbar DL fakelep integration skipped by default"<<endl;

  // Run missing jet integration
  run_missing_jet_integration_ = false;
  pVariable = PyObject_GetAttrString( pModule, "run_missing_jet_integration");
  if (pVariable != NULL )
    run_missing_jet_integration_ = (pVariable==Py_True);
  else
    cout<<"WARNING: run_missing_jet_integration missing, integration skipped for events with integration_type=1 by default"<<endl;

  // Force missing jet integration
  force_missing_jet_integration_ = false;
  pVariable = PyObject_GetAttrString( pModule, "force_missing_jet_integration");
  if (pVariable != NULL )
    force_missing_jet_integration_ = (pVariable==Py_True);
  else
    cout<<"WARNING: force_missing_jet_integration missing, integration type defined event by event by default"<<endl;

  // Force missing jet integration if no MEM sel perm passing cuts without missing jet
  force_missing_jet_integration_ifnoperm_ = true;
  pVariable = PyObject_GetAttrString( pModule, "force_missing_jet_integration_ifnoperm");
  if (pVariable != NULL )
    force_missing_jet_integration_ifnoperm_ = (pVariable==Py_True);
  else
    cout<<"WARNING: force_missing_jet_integration_ifnoperm missing: if integration_type=0 and no permutation passes the compatibility checks, integration_type is changed to 1 by default"<<endl;

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
  if (pVariable != NULL )
    flagTFLepTau_ = (pVariable==Py_True); 
  else
    cout<<"WARNING: flagTFLepTau missing, TF for lep tau enabled by default"<<endl;

  pVariable = PyObject_GetAttrString( pModule, "flagTFHadTau");
  if (pVariable != NULL )
    flagTFHadTau_ = (pVariable==Py_True);   
  else
    cout<<"WARNING: flagTFHadTau missing, TF for had tau enabled by default"<<endl;

  pVariable = PyObject_GetAttrString( pModule, "flagTFMET");
  if (pVariable != NULL )  
    flagTFMET_    = (pVariable==Py_True);  
  else
    cout<<"WARNING: flagTFMET missing, TF for MET enabled by default"<<endl;
  

  pVariable = PyObject_GetAttrString( pModule, "flagTFJet1");
  if (pVariable != NULL )  
    flagTFJet1_   = (pVariable==Py_True);
  else
    cout<<"WARNING: flagTFJet1 missing, TF for Jet1 from hadtop enabled by default"<<endl;  

  pVariable = PyObject_GetAttrString( pModule, "flagTFJet2");
  if (pVariable != NULL )  
    flagTFJet2_   = (pVariable==Py_True); 
  else
    cout<<"WARNING: flagTFJet2 missing, TF for Jet2 from hadtop enabled by default"<<endl; 
  
  pVariable = PyObject_GetAttrString( pModule, "flagTFBJet_leptop");
  if (pVariable != NULL )  
    flagTFBJet_leptop_   = (pVariable==Py_True);
  else
    cout<<"WARNING: flagTFBJet_leptop missing, TF for b-jet from leptop enabled by default"<<endl; 

  pVariable = PyObject_GetAttrString( pModule, "flagTFBJet_hadtop");
  if (pVariable != NULL )  
    flagTFBJet_hadtop_   = (pVariable==Py_True);
  else
    cout<<"WARNING: flagTFBJet_hadtop missing, TF for b-jet from hadtop enabled by default"<<endl; 
  
  pVariable = PyObject_GetAttrString( pModule, "flagTF_fakelep");
  if (pVariable != NULL )  
    flagTF_fakelep_   = (pVariable==Py_True);
  else
    cout<<"WARNING: flagTF_fakelep missing, TF for fake lepton from b in ttbar_DL_fakelep weight enabled by default"<<endl; 

  pVariable = PyObject_GetAttrString( pModule, "flagTF_fakeleptau");
  if (pVariable != NULL )  
    flagTF_fakeleptau_   = (pVariable==Py_True);
  else
    cout<<"WARNING: flagTF_fakeleptau missing, TF for fake tau from lepton in ttZ Zll weight enabled by default"<<endl;

  pVariable = PyObject_GetAttrString( pModule, "flagTFTop");
  if (pVariable != NULL )  
    flagTFTop_   = (pVariable==Py_True);
  else
    cout<<"WARNING: flagTFTop missing, TF for top decay enabled by default"<<endl;

  pVariable = PyObject_GetAttrString( pModule, "flagJac");
  if (pVariable != NULL )  
    flagJac_      = (pVariable==Py_True);  
  else
    cout<<"WARNING: flagJac missing, Jacobian set to 1 by default"<<endl;

  pVariable = PyObject_GetAttrString( pModule, "flagWME");
  if (pVariable != NULL )  
    flagWME_      = (pVariable==Py_True);
  else
    cout<<"WARNING: flagWME missing, Matrix Element set to 1 by default"<<endl;

  // MET TF
  include_hadrecoil_ = true;
  pVariable = PyObject_GetAttrString( pModule, "include_hadrecoil");
  if (pVariable != NULL )  
    include_hadrecoil_  = PyInt_AsLong( pVariable ); 
  else
    cout<<"WARNING: include_hadrecoil missing, hadronic recoil included by default"<<endl;

  // Verbose mode  
  pVariable = PyObject_GetAttrString( pModule, "verbose");
  verbose_  = PyInt_AsLong( pVariable ); 
  pVariable = PyObject_GetAttrString( pModule, "logFileName");
  if(pVariable!=0)
    logFileName_ = PyString_AsString( pVariable );

  // MC integration
  nbrOfDim_ttH_ = 5;
  nbrOfDim_ttZ_ = 5;
  nbrOfDim_ttZ_Zll_ = 3;
  nbrOfDim_ttbar_DL_ = 5;

  nbrOfDim_ttH_miss_ = 7;
  nbrOfDim_ttZ_miss_ = 7;
  nbrOfDim_ttZ_Zll_miss_ = 5;

  nbrOfPoints_ttH_ = 2500;
  nbrOfPoints_ttZ_ = 2500;
  nbrOfPoints_ttZ_Zll_ = 500;
  nbrOfPoints_ttbar_DL_ = 10000;

  nbrOfPoints_ttH_miss_ = 20000;
  nbrOfPoints_ttZ_miss_ = 20000;
  nbrOfPoints_ttZ_Zll_miss_ = 5000;

  pVariable = PyObject_GetAttrString( pModule, "nbrOfDimttH");
  if (pVariable != NULL )  
    nbrOfDim_ttH_ = PyInt_AsLong( pVariable ); 
  else
    cout<<"WARNING: nbrOfDimttH missing, set to "<<nbrOfDim_ttH_<<" by default"<<endl;

  pVariable    = PyObject_GetAttrString( pModule, "nbrOfPoints_ttH");
  if (pVariable != NULL )  
    nbrOfPoints_ttH_ = PyLong_AsLong( pVariable );
  else
    cout<<"WARNING: nbrOfPoints_ttH missing, set to "<<nbrOfPoints_ttH_<<" by default"<<endl;

  if(runttZ_integration_){

    pVariable = PyObject_GetAttrString( pModule, "nbrOfDimttZ");
    if (pVariable != NULL )  
      nbrOfDim_ttZ_ = PyInt_AsLong( pVariable );
    else
      cout<<"WARNING: nbrOfDimttZ missing, set to "<<nbrOfDim_ttZ_<<" by default"<<endl;

    pVariable    = PyObject_GetAttrString( pModule, "nbrOfPoints_ttZ");
    if (pVariable != NULL )  
      nbrOfPoints_ttZ_ = PyLong_AsLong( pVariable );
    else
    cout<<"WARNING: nbrOfPoints_ttZ missing, set to "<<nbrOfPoints_ttZ_<<" by default"<<endl;

  }

  if(runttbar_DL_fakelep_integration_){

    pVariable = PyObject_GetAttrString( pModule, "nbrOfDimttbar_DL");
    if (pVariable != NULL )  
      nbrOfDim_ttbar_DL_ = PyInt_AsLong( pVariable );
    else
      cout<<"WARNING: nbrOfDimttbar_DL missing, set to "<<nbrOfDim_ttbar_DL_<<" by default"<<endl;

    pVariable    = PyObject_GetAttrString( pModule, "nbrOfPoints_ttbar_DL");
    if (pVariable != NULL )  
      nbrOfPoints_ttbar_DL_ = PyLong_AsLong( pVariable );
    else
      cout<<"WARNING: nbrOfPoints_ttbar_DL missing, set to "<<nbrOfPoints_ttbar_DL_<<" by default"<<endl;

  }

  if(runttZ_Zll_integration_){
    
    pVariable = PyObject_GetAttrString( pModule, "nbrOfDimttZ_Zll");
    if (pVariable != NULL )  
      nbrOfDim_ttZ_Zll_ = PyInt_AsLong( pVariable );
    else
      cout<<"WARNING: nbrOfDimttZ_Zll missing, set to "<<nbrOfDim_ttZ_Zll_<<" by default"<<endl;
        
    pVariable    = PyObject_GetAttrString( pModule, "nbrOfPoints_ttZ_Zll");
    if (pVariable != NULL )  
      nbrOfPoints_ttZ_Zll_ = PyLong_AsLong( pVariable );
    else
      cout<<"WARNING: nbrOfPoints_ttZ_Zll missing, set to "<<nbrOfPoints_ttZ_Zll_<<" by default"<<endl;

  }


  if(run_missing_jet_integration_){

    pVariable = PyObject_GetAttrString( pModule, "nbrOfDimttH_miss");
    if (pVariable != NULL )  
      nbrOfDim_ttH_miss_ = PyInt_AsLong( pVariable ); 
    else
      cout<<"WARNING: nbrOfDimttH_miss missing, set to "<<nbrOfDim_ttH_miss_<<" by default"<<endl;

    pVariable    = PyObject_GetAttrString( pModule, "nbrOfPoints_ttH_miss");
    if (pVariable != NULL )  
      nbrOfPoints_ttH_miss_ = PyLong_AsLong( pVariable ); 
    else
      cout<<"WARNING: nbrOfPoints_ttH_miss missing, set to "<<nbrOfPoints_ttH_miss_<<" by default"<<endl;
 
    if(runttZ_integration_){
     
      pVariable = PyObject_GetAttrString( pModule, "nbrOfDimttZ_miss");
      if (pVariable != NULL )  
	nbrOfDim_ttZ_miss_ = PyInt_AsLong( pVariable );
      else
	cout<<"WARNING: nbrOfDimttZ_miss missing, set to "<<nbrOfDim_ttZ_miss_<<" by default"<<endl;
      
      pVariable    = PyObject_GetAttrString( pModule, "nbrOfPoints_ttZ_miss");
      if (pVariable != NULL )  
	nbrOfPoints_ttZ_miss_ = PyLong_AsLong( pVariable );
      else
	cout<<"WARNING: nbrOfPoints_ttZ_miss missing, set to "<<nbrOfPoints_ttZ_miss_<<" by default"<<endl;      

    }

    if(runttZ_Zll_integration_){

      pVariable = PyObject_GetAttrString( pModule, "nbrOfDimttZ_Zll_miss");
      if (pVariable != NULL )  
	nbrOfDim_ttZ_Zll_miss_ = PyInt_AsLong( pVariable );
      else
	cout<<"WARNING: nbrOfDimttZ_Zll_miss missing, set to "<<nbrOfDim_ttZ_miss_<<" by default"<<endl;

      pVariable    = PyObject_GetAttrString( pModule, "nbrOfPoints_ttZ_Zll_miss");
      if (pVariable != NULL )  
	nbrOfPoints_ttZ_Zll_miss_ = PyLong_AsLong( pVariable );
      else
	cout<<"WARNING: nbrOfPoints_ttZ_Zll_miss missing, set to "<<nbrOfPoints_ttZ_Zll_miss_<<" by default"<<endl;      


    }

  }
  

  nbrOfPermut_per_jet_ = 4;
  pVariable    = PyObject_GetAttrString( pModule, "nbrOfPermut_per_jet");
  if (pVariable != NULL )  
    nbrOfPermut_per_jet_ = PyLong_AsLong( pVariable );
  else
    cout<<"WARNING: nbrOfPermut_per_jet missing, set to "<<nbrOfPermut_per_jet_<<" by default"<<endl; 


  use_pT_TFJet_ = true;
  pVariable = PyObject_GetAttrString( pModule, "use_pT_TFJet");
  if (pVariable != NULL )  
    use_pT_TFJet_ = (pVariable==Py_True);
  else
    cout<<"WARNING: use_pT_TFJet missing, use pT-parametrized jet TF by default"<<endl; 

  use_top_compatibility_check_ = true;
  pVariable = PyObject_GetAttrString( pModule, "use_top_compatibility_check");
  if (pVariable != NULL )  
    use_top_compatibility_check_ = (pVariable==Py_True);
  else
    cout<<"WARNING: use_top_compatibility_check missing, set to true by default (integration performed even if jets + leptons not compatible with mass constraints)"<<endl; 


  //Integration boundaries
  CI_TFJet_ = 0.95;
  pVariable = PyObject_GetAttrString( pModule, "CI_TFJet");
  if (pVariable != NULL )  
    CI_TFJet_ = PyFloat_AsDouble( pVariable );
  else
    cout<<"WARNING: CI_TFJet missing, 95% CI of jet TF used by default to define integration boundaries"<<endl; 

  eta_acceptance_ = 2.4;
  pVariable = PyObject_GetAttrString( pModule, "eta_acceptance");
  if (pVariable != NULL )
    eta_acceptance_ = PyFloat_AsDouble( pVariable );
  else
    cout<<"WARNING: eta_acceptance missing, eta_max=2.4 used by default for missing jet integration"<<endl; 

  jet_radius_ = 0.4;
  pVariable = PyObject_GetAttrString( pModule, "jet_radius");
  if (pVariable != NULL )
    jet_radius_ = PyFloat_AsDouble( pVariable );
  else
    cout<<"WARNING: jet_radius missing, R=0.4 used by default for missing jet integration"<<endl; 

  dR_veto_jet_lep_ = 0.4;
  pVariable = PyObject_GetAttrString( pModule, "dR_veto_jet_lep");
  if (pVariable != NULL )
    dR_veto_jet_lep_ = PyFloat_AsDouble( pVariable );
  else
    cout<<"WARNING: dR_veto_jet_lep missing, dR=0.4 used by default for missing jet integration"<<endl; 

  rel_iso_lep_ = 0.4;
  pVariable = PyObject_GetAttrString( pModule, "rel_iso_lep");
  if (pVariable != NULL )
    rel_iso_lep_ = PyFloat_AsDouble( pVariable );
  else
    cout<<"WARNING: rel_iso_lep missing, rel_iso=0.4 used by default for missing jet integration"<<endl; 

  pT_cut_ = 25.;
  pVariable = PyObject_GetAttrString( pModule, "pT_cut");
  if (pVariable != NULL )
    pT_cut_ = PyFloat_AsDouble( pVariable );
  else
    cout<<"WARNING: pT_cut missing, pT_min=25 GeV used by default for missing jet integration"<<endl; 


  pVariable = PyObject_GetAttrString( pModule, "PI");
  double pi = PyFloat_AsDouble( pVariable );

  if(run_missing_jet_integration_){

    phi_missing_jet_Boundaries_[0] = 0;
    phi_missing_jet_Boundaries_[1] = 2*pi;

    pVariable    = PyObject_GetAttrString( pModule, "phi_missing_jet");
    if (pVariable != NULL ){
      pItem = PyList_GetItem(pVariable, 0);
      phi_missing_jet_Boundaries_[0] =  PyFloat_AsDouble( pItem );
      pItem = PyList_GetItem(pVariable, 1);
      phi_missing_jet_Boundaries_[1] =  PyFloat_AsDouble( pItem );
    }
    else
      cout<<"WARNING: phi_missing_jet missing, ["<<phi_missing_jet_Boundaries_[0]<<","<<phi_missing_jet_Boundaries_[1]<<"] used by default for missing jet integration"<<endl; 
    
    cosTheta_missing_jet_Boundaries_[0] = -1.;
    cosTheta_missing_jet_Boundaries_[1] = 1.;      
    pVariable    = PyObject_GetAttrString( pModule, "cosTheta_missing_jet");
    if (pVariable != NULL ){
      pItem = PyList_GetItem(pVariable, 0);
      cosTheta_missing_jet_Boundaries_[0] =  PyFloat_AsDouble( pItem );
      pItem = PyList_GetItem(pVariable, 1);
      cosTheta_missing_jet_Boundaries_[1] =  PyFloat_AsDouble( pItem );
    }
   else
     cout<<"WARNING: cosTheta_missing_jet missing, ["<<cosTheta_missing_jet_Boundaries_[0]<<","<<cosTheta_missing_jet_Boundaries_[1]<<"] used by default for missing jet integration"<<endl; 

  }


  phiNu_tlep_Boundaries_[0] = 0;
  phiNu_tlep_Boundaries_[1] = 2*pi;
  pVariable    = PyObject_GetAttrString( pModule, "phiNu_tlep");
  if (pVariable != NULL ){
    pItem = PyList_GetItem(pVariable, 0);
    phiNu_tlep_Boundaries_[0] =  PyFloat_AsDouble( pItem );
    pItem = PyList_GetItem(pVariable, 1);
    phiNu_tlep_Boundaries_[1] =  PyFloat_AsDouble( pItem );
  }
  else
    cout<<"WARNING: phiNu_tlep missing, ["<<phiNu_tlep_Boundaries_[0]<<","<<phiNu_tlep_Boundaries_[1]<<"] used by default"<<endl; 

  cosTheta_missing_jet_Boundaries_[0] = -1.;
  cosTheta_missing_jet_Boundaries_[1] = 1.; 
  pVariable    = PyObject_GetAttrString( pModule, "cosThetaNu_tlep");
  if (pVariable != NULL ){
    pItem = PyList_GetItem(pVariable, 0);
    cosThetaNu_tlep_Boundaries_[0] =  PyFloat_AsDouble( pItem );
    pItem = PyList_GetItem(pVariable, 1);
    cosThetaNu_tlep_Boundaries_[1] =  PyFloat_AsDouble( pItem );
  }
  else
    cout<<"WARNING: cosThetaNu_tlep missing, ["<<cosThetaNu_tlep_Boundaries_[0]<<","<<cosThetaNu_tlep_Boundaries_[1]<<"] used by default"<<endl; 


  if( runttbar_DL_fakelep_integration_ ){
    
    phiNu_ttau_ttbar_DL_Boundaries_[0] = 0;
    phiNu_ttau_ttbar_DL_Boundaries_[1] = 2*pi;
    pVariable    = PyObject_GetAttrString( pModule, "phiNu_ttau");
    if (pVariable != NULL ){
      pItem = PyList_GetItem(pVariable, 0);
      phiNu_ttau_ttbar_DL_Boundaries_[0] =  PyFloat_AsDouble( pItem );
      pItem = PyList_GetItem(pVariable, 1);
      phiNu_ttau_ttbar_DL_Boundaries_[1] =  PyFloat_AsDouble( pItem );
    }
    else
      cout<<"WARNING: phiNu_ttau missing, ["<<phiNu_ttau_ttbar_DL_Boundaries_[0]<<","<<phiNu_ttau_ttbar_DL_Boundaries_[1]<<"] used by default"<<endl; 

    cosThetaNu_ttau_ttbar_DL_Boundaries_[0] = -1.;
    cosThetaNu_ttau_ttbar_DL_Boundaries_[1] = 1.;    
    pVariable    = PyObject_GetAttrString( pModule, "cosThetaNu_ttau");
    if (pVariable != NULL ){
      pItem = PyList_GetItem(pVariable, 0);
      cosThetaNu_ttau_ttbar_DL_Boundaries_[0] =  PyFloat_AsDouble( pItem );
      pItem = PyList_GetItem(pVariable, 1);
      cosThetaNu_ttau_ttbar_DL_Boundaries_[1] =  PyFloat_AsDouble( pItem );
    }
    else
      cout<<"WARNING: cosThetaNu_ttau missing, ["<<cosThetaNu_ttau_ttbar_DL_Boundaries_[0]<<","<<cosThetaNu_ttau_ttbar_DL_Boundaries_[1]<<"] used by default"<<endl; 

  }

  sqrtS_ = 13000;
  pVariable = PyObject_GetAttrString( pModule, "sqrtS");
  if (pVariable != NULL )
    sqrtS_    = PyFloat_AsDouble( pVariable );   
  else
    cout<<"WARNING: sqrtS missing, 13 TeV used by default"<<endl;

  cout<<endl;

  
}
