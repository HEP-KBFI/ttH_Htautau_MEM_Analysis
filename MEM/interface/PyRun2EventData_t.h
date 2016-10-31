#ifndef PYRUN2EVENTDATA_T_H
#define	PYRUN2EVENTDATA_T_H

# include <stdint.h>
# include <vector>
# include "PyConfig.h"


using namespace std;


typedef struct {

  // 2 selected lepton
  vector<float> _recolep_sel_px;
  vector<float> _recolep_sel_py;
  vector<float> _recolep_sel_pz;
  vector<float> _recolep_sel_e;
  vector<int> _recolep_sel_pdg;

  // Selected recotauh
  vector<float> _recotauh_sel_px;
  vector<float> _recotauh_sel_py;
  vector<float> _recotauh_sel_pz;
  vector<float> _recotauh_sel_e;
  vector<int> _recotauh_sel_decayMode;
  
  //MET info
  float _PFMET;
  float _PFMET_phi;
  float _PFMETCov00;
  float _PFMETCov01;
  float _PFMETCov10;
  float _PFMETCov11;

  // 2 selected b-jets (highest CSV score)
  vector<float> _recoPFJet_btag_px;
  vector<float> _recoPFJet_btag_py;
  vector<float> _recoPFJet_btag_pz;
  vector<float> _recoPFJet_btag_e;
  
  // Untagged jets: all jets but 2 b-jets, pT-ordered
  int _n_recoPFJet_untag;
  vector<float> _recoPFJet_untag_px;
  vector<float> _recoPFJet_untag_py;
  vector<float> _recoPFJet_untag_pz;
  vector<float> _recoPFJet_untag_e;


} EventRead_t;



class PyRun2EventData_t {

 public :

  int _integration_type;

  vector<float> _recolep_sel_px;
  vector<float> _recolep_sel_py;
  vector<float> _recolep_sel_pz;
  vector<float> _recolep_sel_e;
  vector<int> _recolep_sel_pdg;

  vector<float> _recotauh_sel_px;
  vector<float> _recotauh_sel_py;
  vector<float> _recotauh_sel_pz;
  vector<float> _recotauh_sel_e;
  vector<int> _recotauh_sel_decayMode;
  
  float _PFMET;
  float _PFMET_phi;
  float _PFMETCov00;
  float _PFMETCov01;
  float _PFMETCov10;
  float _PFMETCov11;

  vector<float> _recoPFJet_btag_px;
  vector<float> _recoPFJet_btag_py;
  vector<float> _recoPFJet_btag_pz;
  vector<float> _recoPFJet_btag_e;
  
  int _n_recoPFJet_untag;
  vector<float> _recoPFJet_untag_px;
  vector<float> _recoPFJet_untag_py;
  vector<float> _recoPFJet_untag_pz;
  vector<float> _recoPFJet_untag_e;


  // Python fields
  PyObject *pModule_;
  PyObject *pEvents_;
  uint64_t nbrEventsInPyFile_;
  EventRead_t *evRead_;
  
  PyRun2EventData_t( ) {};
  PyRun2EventData_t( const PyRun2EventData_t *evData); 
};



#endif	/* PYRUN2EVENTDATA_T_H */
