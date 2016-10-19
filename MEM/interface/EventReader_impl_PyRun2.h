#ifndef EVENTREADER_IMPL_PYRUN2_H
#define	EVENTREADER_IMPL_PYRUN2_H

#include "ttH_Htautau_MEM_Analysis/MEM/interface/PyRun2EventData_t.h"

template <> 
EventReader<PyRun2EventData_t>::EventReader( const vector<string> &fileNames, long startEvent, long maxNbrOfEventsToRead  );

template <> 
void EventReader<PyRun2EventData_t>::copyIn(PyRun2EventData_t &eventData);

template <>
bool EventReader<PyRun2EventData_t>::readEvent(
	PyRun2EventData_t &eventData,
	int &integration_type,
	TLorentzVector &Lep1_4P,    int &lepton1_Type,
	TLorentzVector &Lep2_4P,    int &lepton2_Type,    
	TLorentzVector &HadSys_4P,   int &decayMode, 
	TLorentzVector &BJet1_4P, TLorentzVector &BJet2_4P,
	TLorentzVector &Jet1_4P, TLorentzVector &Jet2_4P,
	vector<TLorentzVector> &Jets_4P,
	TLorentzVector &recoMET_4P, double recoMETCov[4]); 

template <>
bool EventReader<PyRun2EventData_t>::fillEvent( IntegrationMsg_t *integration, long nbrOfEvents );

# endif
