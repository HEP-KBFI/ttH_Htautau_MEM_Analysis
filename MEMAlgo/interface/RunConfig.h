/* 
 * File:   RunConfig.h
 * Author: grasseau
 *
 * Created on 19 mai 2015, 11:46
 */

#ifndef RUNCONFIG_H
#define	RUNCONFIG_H

# include <stdint.h>
# include <string>
# include <vector>
# include <climits>
# include <iostream> 

using namespace std;

typedef enum {
  BadEventType  = -1,
  PyRun2_EventType = 0
} EventType_t;

class RunConfig {
  public :
  EventType_t eventType_;
  int64_t maxNbrOfEventsToRead_;
  int64_t skipEvents_;
  // Inputs
  std::string configFileName;
  string treeName_;
  vector<string> inputFileNames_;
  
  string MGParamCardFileName_;
  string LHAPDFFileName_;  

  int verbose_;
  string logFileName_;

  //Dimensions and number of points
  int    nbrOfDim_ttH_;
  int    nbrOfDim_ttZ_;
  int    nbrOfDim_ttZ_Zll_;
  int    nbrOfDim_ttbar_DL_;

  int    nbrOfDim_ttH_miss_;
  int    nbrOfDim_ttZ_miss_;
  int    nbrOfDim_ttZ_Zll_miss_;

  
  int64_t   nbrOfPoints_ttH_;
  int64_t   nbrOfPoints_ttZ_;
  int64_t   nbrOfPoints_ttZ_Zll_;
  int64_t   nbrOfPoints_ttbar_DL_;

  int64_t   nbrOfPoints_ttH_miss_;
  int64_t   nbrOfPoints_ttZ_miss_;
  int64_t   nbrOfPoints_ttZ_Zll_miss_;

  int    nbrOfPermut_per_jet_;

  // Run bkg integration
  bool runttZ_integration_;
  bool runttZ_Zll_integration_;
  bool runttbar_DL_fakelep_integration_;

  bool run_missing_jet_integration_;
  bool force_missing_jet_integration_;
  bool force_missing_jet_integration_ifnoperm_;

// RNG
  bool flagSameRNG_;
  
  // Computation / run
  bool flagTFLepTau_; 
  bool flagTFHadTau_;
  bool flagTFMET_;   
  bool flagTFJet1_;  
  bool flagTFJet2_;
  bool flagTFBJet_leptop_;  
  bool flagTFBJet_hadtop_;
  bool flagTF_fakelep_;   
  bool flagTF_fakeleptau_; 
  bool flagTFTop_;     
  bool flagJac_; 
  bool flagWME_; 

  // MET TF
  bool include_hadrecoil_;

  // Jets TF
  double CI_TFJet_;
  bool use_pT_TFJet_;
  bool use_top_compatibility_check_;

  //Common boundaries
  double phiNu_tlep_Boundaries_[2];
  double cosThetaNu_tlep_Boundaries_[2];
  double phiNu_ttau_ttbar_DL_Boundaries_[2];
  double cosThetaNu_ttau_ttbar_DL_Boundaries_[2];
  double phi_missing_jet_Boundaries_[2];
  double cosTheta_missing_jet_Boundaries_[2];

  // Missing jet
  double eta_acceptance_;
  double jet_radius_;
  double dR_veto_jet_lep_;
  double rel_iso_lep_;
  double pT_cut_;
  
  double sqrtS_;

public :
  RunConfig(const char *configname);
  
};

#endif	/* RUNCONFIG_H */

