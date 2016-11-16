# include "ttH_Htautau_MEM_Analysis/MEM/interface/EventReader.h"
# include "ttH_Htautau_MEM_Analysis/MEM/interface/EventReader_impl_PyRun2.h"

template <> 
EventReader<PyRun2EventData_t>::EventReader( const vector<string> &fileNames, long startEvent, long maxNbrOfEventsToRead  ) {
  
  startEvent_ = startEvent;
  maxNbrOfEventsToRead_ = maxNbrOfEventsToRead;
  currentEvent_ = startEvent;

  _integration_type = 0;

  _recolep_sel_px.clear();;
  _recolep_sel_py.clear();
  _recolep_sel_pz.clear();
  _recolep_sel_e.clear();
  _recolep_sel_pdg.clear();

  _recotauh_sel_px.clear();
  _recotauh_sel_py.clear();
  _recotauh_sel_pz.clear();
  _recotauh_sel_e.clear();
  _recotauh_sel_decayMode.clear();
  
  _PFMET = 0;
  _PFMET_phi = 0;
  _PFMETCov00 = 0;
  _PFMETCov01 = 0;
  _PFMETCov10 = 0;
  _PFMETCov11 = 0;

  _recoPFJet_btag_px.clear();
  _recoPFJet_btag_py.clear();
  _recoPFJet_btag_pz.clear();
  _recoPFJet_btag_e.clear();
  
  _n_recoPFJet_untag = 0;
  _recoPFJet_untag_px.clear();
  _recoPFJet_untag_py.clear();
  _recoPFJet_untag_pz.clear();
  _recoPFJet_untag_e.clear();
  

  if ( fileNames.size() != 1 ) {
    if ( fileNames.size() == 0 ) {
      cerr << "Event_impl_PyRun2: error no file" << endl;
    } else {
      cerr << "Event_impl_PyRun2: warning only '" << fileNames[0] 
           << "' is taken into account" << endl;       
    }
    return;
  } else {
    pModule_ = readPyConfig( fileNames[0].c_str() );
    PyObject* pEvents_ = PyObject_GetAttrString( pModule_, "Events");
    // Process array of events
    if ( pEvents_ ) {
      PyObject* pItem;
      unsigned int i;
      nbrEventsInPyFile_ = PyList_Size(pEvents_);
      if ( nbrEventsInPyFile_ == 0 ) { 
        cerr << "Event_impl_PyRun1: No Events in file" << endl;
        return;
      } else {
        evRead_ = new EventRead_t[nbrEventsInPyFile_];
        // Process for each event
        for( i = 0; i < nbrEventsInPyFile_; ++i) {
          pItem = PyList_GetItem(pEvents_, i);
          //
	  evRead_[i]._recolep_sel_px  = getPyDictVectFloat( pItem, "_recolep_sel_px");
	  evRead_[i]._recolep_sel_py  = getPyDictVectFloat( pItem, "_recolep_sel_py");
	  evRead_[i]._recolep_sel_pz  = getPyDictVectFloat( pItem, "_recolep_sel_pz");
	  evRead_[i]._recolep_sel_e   = getPyDictVectFloat( pItem, "_recolep_sel_e");
	  evRead_[i]._recolep_sel_pdg = getPyDictVectInt( pItem,   "_recolep_sel_pdg");
	  //
	  evRead_[i]._recotauh_sel_px        = getPyDictVectFloat( pItem, "_recotauh_sel_px");
	  evRead_[i]._recotauh_sel_py        = getPyDictVectFloat( pItem, "_recotauh_sel_py");
	  evRead_[i]._recotauh_sel_pz        = getPyDictVectFloat( pItem, "_recotauh_sel_pz");
	  evRead_[i]._recotauh_sel_e         = getPyDictVectFloat( pItem, "_recotauh_sel_e");
	  evRead_[i]._recotauh_sel_decayMode = getPyDictVectInt( pItem,   "_recotauh_sel_decayMode");
          //
          evRead_[i]._PFMET      = getPyDictFloat( pItem, "_PFMET");
          evRead_[i]._PFMET_phi  = getPyDictFloat( pItem, "_PFMET_phi");
          evRead_[i]._PFMETCov00 = getPyDictFloat( pItem, "_PFMETCov00");
          evRead_[i]._PFMETCov01 = getPyDictFloat( pItem, "_PFMETCov01");
          evRead_[i]._PFMETCov10 = getPyDictFloat( pItem, "_PFMETCov10");
          evRead_[i]._PFMETCov11 = getPyDictFloat( pItem, "_PFMETCov11");
	  //
	  evRead_[i]._recoPFJet_btag_px  = getPyDictVectFloat( pItem, "_recoPFJet_btag_px");
	  evRead_[i]._recoPFJet_btag_py  = getPyDictVectFloat( pItem, "_recoPFJet_btag_py");
	  evRead_[i]._recoPFJet_btag_pz  = getPyDictVectFloat( pItem, "_recoPFJet_btag_pz");
	  evRead_[i]._recoPFJet_btag_e   = getPyDictVectFloat( pItem, "_recoPFJet_btag_e");
	  //
	  evRead_[i]._n_recoPFJet_untag   = getPyDictInt( pItem, "_n_recoPFJet_untag");
	  evRead_[i]._recoPFJet_untag_px  = getPyDictVectFloat( pItem, "_recoPFJet_untag_px");
	  evRead_[i]._recoPFJet_untag_py  = getPyDictVectFloat( pItem, "_recoPFJet_untag_py");
	  evRead_[i]._recoPFJet_untag_pz  = getPyDictVectFloat( pItem, "_recoPFJet_untag_pz");
	  evRead_[i]._recoPFJet_untag_e   = getPyDictVectFloat( pItem, "_recoPFJet_untag_e");

        }
      }
    }
  }
}




template <> 
void EventReader<PyRun2EventData_t>::copyIn(PyRun2EventData_t &eventData) {
  
  eventData._recolep_sel_px  = _recolep_sel_px;
  eventData._recolep_sel_py  = _recolep_sel_py;
  eventData._recolep_sel_pz  = _recolep_sel_pz;
  eventData._recolep_sel_e   = _recolep_sel_e;
  eventData._recolep_sel_pdg = _recolep_sel_pdg;
  //
  eventData._recotauh_sel_px = _recotauh_sel_px;
  eventData._recotauh_sel_py = _recotauh_sel_py;
  eventData._recotauh_sel_pz = _recotauh_sel_pz;
  eventData._recotauh_sel_e  = _recotauh_sel_e;
  eventData._recotauh_sel_decayMode = _recotauh_sel_decayMode;
  //
  eventData._PFMET = _PFMET;
  eventData._PFMET_phi = _PFMET_phi;
  eventData._PFMETCov00 = _PFMETCov00;
  eventData._PFMETCov01 = _PFMETCov01;
  eventData._PFMETCov10 = _PFMETCov10;
  eventData._PFMETCov11 = _PFMETCov11;
  //
  eventData._recoPFJet_btag_px = _recoPFJet_btag_px;
  eventData._recoPFJet_btag_py = _recoPFJet_btag_py;
  eventData._recoPFJet_btag_pz = _recoPFJet_btag_pz;
  eventData._recoPFJet_btag_e  = _recoPFJet_btag_e;
  //
  eventData._n_recoPFJet_untag  = _n_recoPFJet_untag;
  eventData._recoPFJet_untag_px = _recoPFJet_untag_px;
  eventData._recoPFJet_untag_py = _recoPFJet_untag_py;
  eventData._recoPFJet_untag_pz = _recoPFJet_untag_pz;
  eventData._recoPFJet_untag_e  = _recoPFJet_untag_e;

}


# define PyAssign( field ) (ev->field = pyEv->field)



void getEntry(EventReader<PyRun2EventData_t> *ev, uint64_t currentEvent_) {
  
  int pyEvID = currentEvent_ % ev->nbrEventsInPyFile_;
  EventRead_t *pyEv = &ev->evRead_[ pyEvID ];
  
  PyAssign(_recolep_sel_px); PyAssign(_recolep_sel_py); PyAssign(_recolep_sel_pz); PyAssign(_recolep_sel_e);
  PyAssign(_recolep_sel_pdg);
  //
  PyAssign(_recotauh_sel_px); PyAssign(_recotauh_sel_py); PyAssign(_recotauh_sel_pz); PyAssign(_recotauh_sel_e);
  PyAssign(_recotauh_sel_decayMode);
  //
  PyAssign(_PFMET); PyAssign(_PFMET_phi);
  PyAssign(_PFMETCov00); PyAssign(_PFMETCov01); PyAssign(_PFMETCov10); PyAssign(_PFMETCov11);
  //
  PyAssign(_recoPFJet_btag_px); PyAssign(_recoPFJet_btag_py); PyAssign(_recoPFJet_btag_pz); PyAssign(_recoPFJet_btag_e);
  //
  PyAssign(_n_recoPFJet_untag);
  PyAssign(_recoPFJet_untag_px); PyAssign(_recoPFJet_untag_py); PyAssign(_recoPFJet_untag_pz); PyAssign(_recoPFJet_untag_e);

}


template <>
bool EventReader<PyRun2EventData_t>::readEvent(PyRun2EventData_t &eventData,
					     int &integration_type,
					     TLorentzVector &Lep1_4P,    int &lepton1_Type,
					     TLorentzVector &Lep2_4P,    int &lepton2_Type,    
					     TLorentzVector &HadSys_4P,   int &decayMode, 
					     TLorentzVector &BJet1_4P, TLorentzVector &BJet2_4P,
					     TLorentzVector &Jet1_4P, TLorentzVector &Jet2_4P,
					     vector<TLorentzVector> &Jets_4P,
					     TLorentzVector &recoMET_4P, double recoMETCov[4]) {
  
  getEntry( this, currentEvent_ );
  
  copyIn( eventData );
  
  Lep1_4P.SetPxPyPzE( (_recolep_sel_px[0]), (_recolep_sel_py[0]), (_recolep_sel_pz[0]), (_recolep_sel_e[0]) );
  lepton1_Type = (_recolep_sel_pdg[0]);
  
  Lep2_4P.SetPxPyPzE( (_recolep_sel_px[1]), (_recolep_sel_py[1]), (_recolep_sel_pz[1]), (_recolep_sel_e[1]) );
  lepton2_Type = (_recolep_sel_pdg[1]);
  
  HadSys_4P.SetPxPyPzE( (_recotauh_sel_px)[0], (_recotauh_sel_py)[0], (_recotauh_sel_pz)[0], (_recotauh_sel_e)[0] );
  decayMode = (_recotauh_sel_decayMode)[0];
  
  recoMET_4P.SetPtEtaPhiM( _PFMET, 0,  _PFMET_phi , 0 ); 
  
  // Set MET covariance Matrix component (index order 00, 01, 10, 00 )
  double det = _PFMETCov00 * _PFMETCov11 - _PFMETCov01 * _PFMETCov10;
  //
  recoMETCov[0] =  (_PFMETCov11 / det );
  recoMETCov[1] = -(_PFMETCov01 / det );
  recoMETCov[2] = -(_PFMETCov10 / det );
  recoMETCov[3] =  (_PFMETCov00 / det );
  
  //Two jets with highest CSV score
  BJet1_4P.SetPxPyPzE( (_recoPFJet_btag_px)[0], (_recoPFJet_btag_py)[0], (_recoPFJet_btag_pz)[0], (_recoPFJet_btag_e)[0] );
  BJet2_4P.SetPxPyPzE( (_recoPFJet_btag_px)[1], (_recoPFJet_btag_py)[1], (_recoPFJet_btag_pz)[1], (_recoPFJet_btag_e)[1] );
  
  //Untagged jets
  Jets_4P.clear();
  for(int i=0; i<min(_n_recoPFJet_untag,10); i++){
    TLorentzVector Jet_4P;
    Jet_4P.SetPxPyPzE( (_recoPFJet_untag_px)[i], (_recoPFJet_untag_py)[i], (_recoPFJet_untag_pz)[i], (_recoPFJet_untag_e)[i] );
    Jets_4P.push_back(Jet_4P);
  }


  if( Jets_4P.size()==1 ){
    integration_type = 1;
    Jet1_4P.SetPxPyPzE( 0., 0., 0., 0. );
    Jet2_4P.SetPxPyPzE( 0., 0., 0., 0. );    
  }


  else{
  
    Jet1_4P = Jets_4P[0];
    Jet2_4P = Jets_4P[1];
    double mW = 80.4;
    double delta_Mjj = fabs( (Jet1_4P+Jet2_4P).M() - mW );
    
    for(unsigned int i=0;i<Jets_4P.size()-1;i++){
      for(unsigned int j=i+1;j<Jets_4P.size();j++){
	
	if( fabs( (Jets_4P[i]+Jets_4P[j]).M() - mW ) < delta_Mjj ){
	  Jet1_4P = Jets_4P[i];
	  Jet2_4P = Jets_4P[j];
	  delta_Mjj = fabs( (Jet1_4P+Jet2_4P).M() - mW );
	}
	
      }
    }
    
    if( (Jet1_4P+Jet2_4P).M()>60 && (Jet1_4P+Jet2_4P).M()<100 )
      integration_type = 0;
    
    else{
      integration_type = 1;
      Jet1_4P.SetPxPyPzE( 0., 0., 0., 0. );
      Jet2_4P.SetPxPyPzE( 0., 0., 0., 0. );    
    }  

  }

  cout<<"Event #"<<currentEvent_<<" selected"<<endl;

  return true;


}


template <>
bool EventReader<PyRun2EventData_t>::fillEvent( IntegrationMsg_t *integration, int64_t nbrOfEvents) {

  int integration_type;
  TLorentzVector Lep1_4P, Lep2_4P, HadSys_4P;
  TLorentzVector BJet1_4P, BJet2_4P, Jet1_4P, Jet2_4P;
  vector<TLorentzVector> Jets_4P;
  TLorentzVector RecoMET_4P;
  double RecoMETCov[4];
  int lepton1_Type, lepton2_Type, tauDecayMode;
  int k = 0;  // Number of event stored
  PyRun2EventData_t eventFields;
  

  for( ; (k != nbrOfEvents) && ( ( maxNbrOfEventsToRead_==-1) || (currentEvent_ < (maxNbrOfEventsToRead_ + startEvent_) ) ); currentEvent_++) {
        
    bool goodEvent = readEvent( eventFields,
				integration_type,
                                Lep1_4P, lepton1_Type,
				Lep2_4P, lepton2_Type,
				HadSys_4P, tauDecayMode, 
				BJet1_4P, BJet2_4P,
				Jet1_4P, Jet2_4P,
				Jets_4P,
				RecoMET_4P, RecoMETCov
                              );
    if (goodEvent) {
      
      // Set integration parameters
      // and leave loop

      integration[k].integration_type_ = integration_type;
      //
      integration[k].evLep1_4P_[0] = Lep1_4P.Px(); integration[k].evLep1_4P_[1] = Lep1_4P.Py(); 
      integration[k].evLep1_4P_[2] = Lep1_4P.Pz(); integration[k].evLep1_4P_[3] = Lep1_4P.E();
      integration[k].lepton1_Type_ = lepton1_Type;
      //
      integration[k].evLep2_4P_[0] = Lep2_4P.Px(); integration[k].evLep2_4P_[1] = Lep2_4P.Py(); 
      integration[k].evLep2_4P_[2] = Lep2_4P.Pz(); integration[k].evLep2_4P_[3] = Lep2_4P.E();
      integration[k].lepton2_Type_ = lepton2_Type;
      //
      integration[k].evHadSys_Tau_4P_[0] = HadSys_4P.Px(); integration[k].evHadSys_Tau_4P_[1] = HadSys_4P.Py();
      integration[k].evHadSys_Tau_4P_[2] = HadSys_4P.Pz(); integration[k].evHadSys_Tau_4P_[3] = HadSys_4P.E();
      integration[k].HadtauDecayMode_ =  tauDecayMode;
      //
      integration[k].evBJet1_4P_[0] = BJet1_4P.Px(); integration[k].evBJet1_4P_[1] = BJet1_4P.Py();
      integration[k].evBJet1_4P_[2] = BJet1_4P.Pz(); integration[k].evBJet1_4P_[3] = BJet1_4P.E();
      integration[k].evBJet2_4P_[0] = BJet2_4P.Px(); integration[k].evBJet2_4P_[1] = BJet2_4P.Py();
      integration[k].evBJet2_4P_[2] = BJet2_4P.Pz(); integration[k].evBJet2_4P_[3] = BJet2_4P.E();
      //
      integration[k].n_lightJets_ = min(10,int(Jets_4P.size()));
      for( int i=0; i<integration[k].n_lightJets_; i++){
	integration[k].evJets_4P_[i][0] = Jets_4P[i].Px(); integration[k].evJets_4P_[i][1] = Jets_4P[i].Py();
	integration[k].evJets_4P_[i][2] = Jets_4P[i].Pz(); integration[k].evJets_4P_[i][3] = Jets_4P[i].E();
      }
      if(Jets_4P.size()<10){	
	for(unsigned int i=Jets_4P.size(); i<10; i++){
	  integration[k].evJets_4P_[i][0] = 0; integration[k].evJets_4P_[i][1] = 0;
	  integration[k].evJets_4P_[i][2] = 0; integration[k].evJets_4P_[i][3] = 0;
	}
      }
      
      if(integration_type == 0){
	integration[k].evJet1_4P_[0] = Jet1_4P.Px(); integration[k].evJet1_4P_[1] = Jet1_4P.Py();
	integration[k].evJet1_4P_[2] = Jet1_4P.Pz(); integration[k].evJet1_4P_[3] = Jet1_4P.E();
	integration[k].evJet2_4P_[0] = Jet2_4P.Px(); integration[k].evJet2_4P_[1] = Jet2_4P.Py();
	integration[k].evJet2_4P_[2] = Jet2_4P.Pz(); integration[k].evJet2_4P_[3] = Jet2_4P.E();
      }
      else if(integration_type == 1){
	integration[k].evJet1_4P_[0] = 0; integration[k].evJet1_4P_[1] = 0;
	integration[k].evJet1_4P_[2] = 0; integration[k].evJet1_4P_[3] = 0;
	integration[k].evJet2_4P_[0] = 0; integration[k].evJet2_4P_[1] = 0;
	integration[k].evJet2_4P_[2] = 0; integration[k].evJet2_4P_[3] = 0;	
      }
      //
      integration[k].evRecoMET4P_[0] = RecoMET_4P.Px(); integration[k].evRecoMET4P_[1] = RecoMET_4P.Py();
      integration[k].evRecoMET4P_[2] = RecoMET_4P.Pz(); integration[k].evRecoMET4P_[3] = RecoMET_4P.E();
      //
      integration[k].evV_[0] = RecoMETCov[0]; integration[k].evV_[1] = RecoMETCov[1]; 
      integration[k].evV_[2] = RecoMETCov[2]; integration[k].evV_[3] = RecoMETCov[3];     
    
      k++;
    } 

  } //end loop over the events
 
  return (k != 0);
    
}




template <> 
EventReader<PyRun2EventData_t>::~EventReader( ) {
  delete [] evRead_;
}


