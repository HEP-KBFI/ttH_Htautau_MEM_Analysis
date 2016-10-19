/* 
 * File:   MCEventReader.h
 * Author: grasseau
 *
 * Created on 16 janvier 2015, 16:11
 */

#ifndef EVENTREADER_H
#define	EVENTREADER_H
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <cstdlib>
#include <TString.h>
#include <algorithm>
#include <map>
#include <sstream>
#include <limits> 
#include "TROOT.h"
#include "TMath.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
# include "ttH_Htautau_MEM_Analysis/MEMAlgo/interface/MGIntegration.h"

using namespace std;

template<typename T>
class EventReader : public T {
 protected:
  uint64_t numberOfEvents_;
  uint64_t currentEvent_;
  uint64_t selectedEvents_;
  uint64_t startEvent_;
  int64_t maxNbrOfEventsToRead_;
  TChain *tchain_;
 public:  
  EventReader() {currentEvent_=0; selectedEvents_=0; maxNbrOfEventsToRead_=0;};
  EventReader( const vector<string>  &fileNames, int64_t startEvent, int64_t maxNbrOfEventsToRead); 
  void copyIn( T &eventData);
 
 bool readEvent( T &eventFields,
		 int& integration_type,
		 TLorentzVector &Lep1_4P,    int &lepton1_Type,
		 TLorentzVector &Lep2_4P,    int &lepton2_Type,    
		 TLorentzVector &HadSys_4P,   int &decayMode, 
		 TLorentzVector &BJet1_4P, TLorentzVector &BJet2_4P,
		 TLorentzVector &Jet1_4P, TLorentzVector &Jet2_4P,
		 vector<TLorentzVector> &Jets_4P,
		 TLorentzVector &recoMET_4P, double recoMETCov[4]) ;
 
  bool fillEvent( IntegrationMsg_t *integration, int64_t nbrOfEvents );
  ~EventReader() {};
};

# include "EventReader_impl_PyRun2.h"
#endif	/* EVENTREADER_H */
