# include "ttH_Htautau_MEM_Analysis/MEM/interface/PyRun2EventData_t.h"

PyRun2EventData_t::PyRun2EventData_t( const PyRun2EventData_t *evData) {

  _recolep_sel_px = evData->_recolep_sel_px;
  _recolep_sel_py = evData->_recolep_sel_py;
  _recolep_sel_pz = evData->_recolep_sel_pz;
  _recolep_sel_e = evData->_recolep_sel_e;
  _recolep_sel_pdg = evData->_recolep_sel_pdg;

  _recotauh_sel_px = evData->_recotauh_sel_px;
  _recotauh_sel_py = evData->_recotauh_sel_py;
  _recotauh_sel_pz = evData->_recotauh_sel_pz;
  _recotauh_sel_e = evData->_recotauh_sel_e;
  _recotauh_sel_decayMode = evData->_recotauh_sel_decayMode;
  
  _PFMET = evData->_PFMET;
  _PFMET_phi =  evData->_PFMET_phi;
  _PFMETCov00 = evData->_PFMETCov00;
  _PFMETCov01 = evData->_PFMETCov01;
  _PFMETCov10 = evData->_PFMETCov10;
  _PFMETCov11 = evData->_PFMETCov11;

  _recoPFJet_btag_px = evData->_recoPFJet_btag_px;
  _recoPFJet_btag_py = evData->_recoPFJet_btag_py;
  _recoPFJet_btag_pz = evData->_recoPFJet_btag_pz;
  _recoPFJet_btag_e = evData->_recoPFJet_btag_e;

  _n_recoPFJet_untag = evData->_n_recoPFJet_untag;
  _recoPFJet_untag_px = evData->_recoPFJet_untag_px;
  _recoPFJet_untag_py = evData->_recoPFJet_untag_py;
  _recoPFJet_untag_pz = evData->_recoPFJet_untag_pz;
  _recoPFJet_untag_e = evData->_recoPFJet_untag_e;
  
  _n_pair_Wtag_recoPFJet_untag = evData->_n_pair_Wtag_recoPFJet_untag;
  _recoPFJet_untag_Wtag_px = evData->_recoPFJet_untag_Wtag_px;
  _recoPFJet_untag_Wtag_py = evData->_recoPFJet_untag_Wtag_py;
  _recoPFJet_untag_Wtag_pz = evData->_recoPFJet_untag_Wtag_pz;
  _recoPFJet_untag_Wtag_e = evData->_recoPFJet_untag_Wtag_e;

}
