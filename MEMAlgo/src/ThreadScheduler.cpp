# include <fstream>
# include <iostream>
# include <iomanip>
# include <string>
# include <sstream>

# include "LHAPDF/LHAPDF.h"

# include "ttH_Htautau_MEM_Analysis/MEMAlgo/interface/ThreadScheduler.h"
# include "ttH_Htautau_MEM_Analysis/MEMAlgo/interface/MGIntegration.h"
# include "ttH_Htautau_MEM_Analysis/MEMAlgo/interface/Processes.h"



ThreadScheduler::ThreadScheduler():NodeScheduler() {
  
  totalNbrOfQueueAndProcess=0;
  remainingTasks=0;
 
  integration = 0;
}

ThreadScheduler::~ThreadScheduler() {
  delete integration;
}


void ThreadScheduler::initNodeScheduler( const RunConfig *cfg,
                          int mpi_rank
                        ) {
 
  NodeScheduler::initNodeScheduler( cfg, mpi_rank );

  integration = new MGIntegration( cfg );
  
  string logNameFile = cfg->logFileName_;
  logStream_.open( logNameFile );
  integration->setLogStream( logStream_ );  


  //
  // Init Vegas with the greatest number of evaluation points
  //

  // GSL init
  rng = gslRNGInit( );
  saveRNGSeed = gsl_rng_clone( rng );
 
  fctDescr_ttH = { 
    &wrapper_evalttH, (size_t) integration->nbrOfDim_ttH_, integration };
  state_ttH = gsl_monte_vegas_alloc ( integration->nbrOfDim_ttH_ );  

  fctDescr_ttZ = { 
    &wrapper_evalttH, (size_t) integration->nbrOfDim_ttZ_, integration };
  state_ttZ = gsl_monte_vegas_alloc ( integration->nbrOfDim_ttZ_ ); 

  fctDescr_ttZ_Zll = { 
    &wrapper_evalttZ_Zll, (size_t) integration->nbrOfDim_ttZ_Zll_, integration };
  state_ttZ_Zll = gsl_monte_vegas_alloc ( integration->nbrOfDim_ttZ_Zll_ );   

  fctDescr_ttbar_DL_fakelep_tlep = { 
    &wrapper_evalttbar_DL_fakelep_tlep, (size_t) integration->nbrOfDim_ttbar_DL_, integration };
  state_ttbar_DL_fakelep_tlep = gsl_monte_vegas_alloc ( integration->nbrOfDim_ttbar_DL_ );  

  fctDescr_ttbar_DL_fakelep_ttau = { 
    &wrapper_evalttbar_DL_fakelep_ttau, (size_t) integration->nbrOfDim_ttbar_DL_, integration };
  state_ttbar_DL_fakelep_ttau = gsl_monte_vegas_alloc ( integration->nbrOfDim_ttbar_DL_ );  

  fctDescr_ttH_miss = { 
    &wrapper_evalttH, (size_t) integration->nbrOfDim_ttH_miss_, integration };
  state_ttH_miss = gsl_monte_vegas_alloc ( integration->nbrOfDim_ttH_miss_ );  

  fctDescr_ttZ_miss = { 
    &wrapper_evalttH, (size_t) integration->nbrOfDim_ttZ_miss_, integration };
  state_ttZ_miss = gsl_monte_vegas_alloc ( integration->nbrOfDim_ttZ_miss_ ); 

  fctDescr_ttZ_Zll_miss = { 
    &wrapper_evalttZ_Zll, (size_t) integration->nbrOfDim_ttZ_Zll_miss_, integration };
  state_ttZ_Zll_miss = gsl_monte_vegas_alloc ( integration->nbrOfDim_ttZ_Zll_miss_ );  
    
  // LHAPDF
  LHAPDF::initPDFSet( PDFSet, cfg->LHAPDFFileName_.c_str());

  // MadGraph
  init_ttH_HTauTauMEProcesses( cfg->configFileName.c_str() ); 

}



void ThreadScheduler::runNodeScheduler (
                eventList_t evList,
                int nbrOfEvents) {
      
  int k; 
  for ( k=0; k < nbrOfEvents; k++) {

    if(integration->force_missing_jet_integration_)
	    integration->integration_type_ = 1;
    
    integration->setEventParameters( evList[k], integration->force_missing_jet_integration_ );
    integration->copyBoundaries( &evList[k]);

    evList[k].weight_ttH_ = 0; 
    evList[k].weight_ttZ_ = 0; 
    evList[k].weight_ttZ_Zll_ = 0; 
    evList[k].weight_ttbar_DL_fakelep_ = 0; 
    

    //
    // Integrate ttbar_DL_fakelep
    //
    
    if(integration->runttbar_DL_fakelep_integration_){
	    
      //No need to loop on light jets for ttbar_DL_fakelep
      for(int perm=0; perm<integration->nbrOfPermut_per_jet_; perm++){
	
	integration->initVersors(perm);   
	
	if( integration->include_perm_ttbar_DL_fakelep_tlep_[perm] ){
	  
	  integration->lowerValues_[cosThetaNu_tlep_ttbar_DL_id] = integration->cosThetaNu_tlep_Boundaries_[0];
	  integration->upperValues_[cosThetaNu_tlep_ttbar_DL_id] = integration->cosThetaNu_tlep_Boundaries_[1];
	  integration->lowerValues_[phiNu_tlep_ttbar_DL_id] = integration->phiNu_tlep_Boundaries_[0];
	  integration->upperValues_[phiNu_tlep_ttbar_DL_id] = integration->phiNu_tlep_Boundaries_[1];
	  integration->lowerValues_[cosThetaNu_ttau_ttbar_DL_id] = integration->cosThetaNu_ttau_ttbar_DL_Boundaries_[0];
	  integration->upperValues_[cosThetaNu_ttau_ttbar_DL_id] = integration->cosThetaNu_ttau_ttbar_DL_Boundaries_[1];
	  integration->lowerValues_[phiNu_ttau_ttbar_DL_id] = integration->phiNu_ttau_ttbar_DL_Boundaries_[0];
	  integration->upperValues_[phiNu_ttau_ttbar_DL_id] = integration->phiNu_ttau_ttbar_DL_Boundaries_[1];
	  integration->lowerValues_[PTauHad_ttbar_DL_id] = integration->PTauHad_ttbar_DL_Boundaries_[0];
	  integration->upperValues_[PTauHad_ttbar_DL_id] = integration->PTauHad_ttbar_DL_Boundaries_[1];	  
	  
	  integration->tot_Drawsttbar_DL_fakelep_tlep_ = 0;
	  integration->integr_Efficiencyttbar_DL_fakelep_tlep_ = 0;
	  
	  // Compute Integral          
	  if ( integration->flagSameRNG_ )
	    // Copy the initial state to the current RNG State
	    gsl_rng_memcpy ( rng, saveRNGSeed );
	  
	  t = clock();

	  gslIntegrate( &fctDescr_ttbar_DL_fakelep_tlep, integration, integration->nbrOfDim_ttbar_DL_, integration->nbrOfPoints_ttbar_DL_, rng, state_ttbar_DL_fakelep_tlep, true,
			&evList[k].integralttbar_DL_fakelep_tlep_[perm], 
			&evList[k].stderrttbar_DL_fakelep_tlep_[perm],
			&evList[k].chiSquarettbar_DL_fakelep_tlep_[perm]);     
	  
	  t = clock() - t;
	  
	  evList[k].compTimettbar_DL_fakelep_tlep_[perm] = ((float)t)/CLOCKS_PER_SEC;	  
	  evList[k].totalDrawsttbar_DL_fakelep_tlep_[perm] = integration->tot_Drawsttbar_DL_fakelep_tlep_;
	  evList[k].integrationEfficiencyttbar_DL_fakelep_tlep_[perm] = integration->integr_Efficiencyttbar_DL_fakelep_tlep_;
	  evList[k].weight_ttbar_DL_fakelep_ += evList[k].integralttbar_DL_fakelep_tlep_[perm];

	  
	}
	if( integration->include_perm_ttbar_DL_fakelep_ttau_[perm] ){
	  
	  integration->lowerValues_[cosThetaNu_tlep_ttbar_DL_id] = integration->cosThetaNu_tlep_Boundaries_[0];
	  integration->upperValues_[cosThetaNu_tlep_ttbar_DL_id] = integration->cosThetaNu_tlep_Boundaries_[1];
	  integration->lowerValues_[phiNu_tlep_ttbar_DL_id] = integration->phiNu_tlep_Boundaries_[0];
	  integration->upperValues_[phiNu_tlep_ttbar_DL_id] = integration->phiNu_tlep_Boundaries_[1];
	  integration->lowerValues_[cosThetaNu_ttau_ttbar_DL_id] = integration->cosThetaNu_ttau_ttbar_DL_Boundaries_[0];
	  integration->upperValues_[cosThetaNu_ttau_ttbar_DL_id] = integration->cosThetaNu_ttau_ttbar_DL_Boundaries_[1];
	  integration->lowerValues_[phiNu_ttau_ttbar_DL_id] = integration->phiNu_ttau_ttbar_DL_Boundaries_[0];
	  integration->upperValues_[phiNu_ttau_ttbar_DL_id] = integration->phiNu_ttau_ttbar_DL_Boundaries_[1];
	  integration->lowerValues_[PTauHad_ttbar_DL_id] = integration->PTauHad_ttbar_DL_Boundaries_[0];
	  integration->upperValues_[PTauHad_ttbar_DL_id] = integration->PTauHad_ttbar_DL_Boundaries_[1];
	  
	  
	  integration->tot_Drawsttbar_DL_fakelep_ttau_ = 0;
	  integration->integr_Efficiencyttbar_DL_fakelep_ttau_ = 0;
	  
	  // Compute Integral          
	  if ( integration->flagSameRNG_ )
	    // Copy the initial state to the current RNG State
	    gsl_rng_memcpy ( rng, saveRNGSeed );
		
	  t = clock();

	  gslIntegrate( &fctDescr_ttbar_DL_fakelep_ttau, integration, integration->nbrOfDim_ttbar_DL_, integration->nbrOfPoints_ttbar_DL_, rng, state_ttbar_DL_fakelep_ttau, true,
			&evList[k].integralttbar_DL_fakelep_ttau_[perm], 
			&evList[k].stderrttbar_DL_fakelep_ttau_[perm],
			&evList[k].chiSquarettbar_DL_fakelep_ttau_[perm]);     

	  t = clock() - t;
	  
	  evList[k].compTimettbar_DL_fakelep_ttau_[perm] = ((float)t)/CLOCKS_PER_SEC;
	  evList[k].totalDrawsttbar_DL_fakelep_ttau_[perm] = integration->tot_Drawsttbar_DL_fakelep_ttau_;
	  evList[k].integrationEfficiencyttbar_DL_fakelep_ttau_[perm] = integration->integr_Efficiencyttbar_DL_fakelep_ttau_;
	  evList[k].weight_ttbar_DL_fakelep_ += evList[k].integralttbar_DL_fakelep_ttau_[perm];
	  
	}
	
      }
      
    }    
    


    //No missing jet
    if(integration->integration_type_ == 0 && !integration->force_missing_jet_integration_){

      for(int perm=0; perm<integration->nbrOfPermut_; perm++){
	
	integration->initVersors(perm);          
	
	//
	// Integrate ttH
	//	

	if( integration->include_perm_ttH_[perm] ){
	  
	  integration->signalME_ = true;
	  integration->m_TauTau_2_ = pow(integration->mTauTau_ttH_[perm],2);
	  integration->lowerValues_[PTauLep_id] = integration->PTauLep_ttH_Lower_[perm];
	  integration->upperValues_[PTauLep_id] = integration->PTauLep_ttH_Upper_[perm];
	  integration->lowerValues_[cosThetaTauLepTauHad_id] = integration->cosTheta_diTau_ttH_Lower_[perm];
	  integration->upperValues_[cosThetaTauLepTauHad_id] = integration->cosTheta_diTau_ttH_Upper_[perm];
	  integration->lowerValues_[EQuark1_id] = integration->EQuark1_Lower_[perm];
	  integration->upperValues_[EQuark1_id] = integration->EQuark1_Upper_[perm];
	  integration->lowerValues_[cosThetaNu_tlep_id] = integration->cosThetaNu_tlep_Boundaries_[0];
	  integration->upperValues_[cosThetaNu_tlep_id] = integration->cosThetaNu_tlep_Boundaries_[1];
	  integration->lowerValues_[phiNu_tlep_id] = integration->phiNu_tlep_Boundaries_[0];
	  integration->upperValues_[phiNu_tlep_id] = integration->phiNu_tlep_Boundaries_[1];
	  

	  integration->integr_EfficiencyttH_ = 0;
	  integration->tot_DrawsttH_ = 0;
	  
	  // Compute Integral          
	  // Copy the initial state to the current RNG State
	  if ( integration->flagSameRNG_ )
	    gsl_rng_memcpy ( rng, saveRNGSeed);
	  
	  t = clock();
	  
	  gslIntegrate( &fctDescr_ttH, integration,  integration->nbrOfDim_ttH_, integration->nbrOfPoints_ttH_, rng, state_ttH, true,
			&evList[k].integralttH_[perm], 
			&evList[k].stderrttH_[perm],
			&evList[k].chiSquarettH_[perm]);
	  
	  
	  t = clock() - t;
	  
	  evList[k].compTimettH_[perm] = ((float)t)/CLOCKS_PER_SEC;
	  evList[k].totalDrawsttH_[perm] = integration->tot_DrawsttH_;
	  evList[k].integrationEfficiencyttH_[perm] = integration->integr_EfficiencyttH_;
	  evList[k].weight_ttH_ += evList[k].integralttH_[perm];

	}
	
	
	//
	// Integrate ttZ
	//
	
	if(integration->runttZ_integration_){
	  
	  if( integration->include_perm_ttZ_[perm] ){
	    integration->signalME_ = false;
	    integration->m_TauTau_2_ = pow(integration->mTauTau_ttZ_[perm],2);
	    integration->lowerValues_[PTauLep_id] = integration->PTauLep_ttZ_Lower_[perm];
	    integration->upperValues_[PTauLep_id] = integration->PTauLep_ttZ_Upper_[perm];
	    integration->lowerValues_[cosThetaTauLepTauHad_id] = integration->cosTheta_diTau_ttZ_Lower_[perm];
	    integration->upperValues_[cosThetaTauLepTauHad_id] = integration->cosTheta_diTau_ttZ_Upper_[perm];
	    integration->lowerValues_[EQuark1_id] = integration->EQuark1_Lower_[perm];
	    integration->upperValues_[EQuark1_id] = integration->EQuark1_Upper_[perm];
	    integration->lowerValues_[cosThetaNu_tlep_id] = integration->cosThetaNu_tlep_Boundaries_[0];
	    integration->upperValues_[cosThetaNu_tlep_id] = integration->cosThetaNu_tlep_Boundaries_[1];
	    integration->lowerValues_[phiNu_tlep_id] = integration->phiNu_tlep_Boundaries_[0];
	    integration->upperValues_[phiNu_tlep_id] = integration->phiNu_tlep_Boundaries_[1];
	    
	    integration->integr_EfficiencyttZ_ = 0;
	    integration->tot_DrawsttZ_ = 0;
	    
	    // Compute Integral          
	    if ( integration->flagSameRNG_ )
	      // Copy the initial state to the current RNG State
	      gsl_rng_memcpy ( rng, saveRNGSeed );
	    
	    t = clock();
	    
	    gslIntegrate( &fctDescr_ttZ, integration, integration->nbrOfDim_ttZ_, integration->nbrOfPoints_ttZ_, rng, state_ttZ, true,
			  &evList[k].integralttZ_[perm], 
			  &evList[k].stderrttZ_[perm],
			  &evList[k].chiSquarettZ_[perm]);     
	    
	    t = clock() - t;

	    evList[k].compTimettZ_[perm] =  ((float)t)/CLOCKS_PER_SEC;
	    evList[k].integrationEfficiencyttZ_[perm] = integration->integr_EfficiencyttZ_;
	    evList[k].totalDrawsttZ_[perm] = integration->tot_DrawsttZ_;
	    evList[k].weight_ttZ_ += evList[k].integralttZ_[perm];

	  }
	
	}


	//
	// Integrate ttZ Zll
	//
	
	if(integration->runttZ_Zll_integration_){
	  
	  if( integration->include_perm_ttZ_Zll_[perm] ){
	    
	    integration->lowerValues_[EQuark1_ttZ_Zll_id] = integration->EQuark1_Lower_[perm];
	    integration->upperValues_[EQuark1_ttZ_Zll_id] = integration->EQuark1_Upper_[perm];
	    integration->lowerValues_[cosThetaNu_tlep_ttZ_Zll_id] = integration->cosThetaNu_tlep_Boundaries_[0];
	    integration->upperValues_[cosThetaNu_tlep_ttZ_Zll_id] = integration->cosThetaNu_tlep_Boundaries_[1];
	    integration->lowerValues_[phiNu_tlep_ttZ_Zll_id] = integration->phiNu_tlep_Boundaries_[0];
	    integration->upperValues_[phiNu_tlep_ttZ_Zll_id] = integration->phiNu_tlep_Boundaries_[1];	
	    
	    integration->tot_DrawsttZ_Zll_ = 0;
	    integration->integr_EfficiencyttZ_Zll_ = 0;
	    
	    // Compute Integral          
	    if ( integration->flagSameRNG_ )
	      // Copy the initial state to the current RNG State
	      gsl_rng_memcpy ( rng, saveRNGSeed );
	    
	    t = clock();
	    
	    gslIntegrate( &fctDescr_ttZ_Zll, integration, integration->nbrOfDim_ttZ_Zll_, integration->nbrOfPoints_ttZ_Zll_, rng, state_ttZ_Zll, true,
			  &evList[k].integralttZ_Zll_[perm], 
			  &evList[k].stderrttZ_Zll_[perm],
			  &evList[k].chiSquarettZ_Zll_[perm]);     
	    
	    t = clock()-t;

	    evList[k].compTimettZ_Zll_[perm] =  ((float)t)/CLOCKS_PER_SEC;
	    evList[k].totalDrawsttZ_Zll_[perm] = integration->tot_DrawsttZ_Zll_;
	    evList[k].integrationEfficiencyttZ_Zll_[perm] = integration->integr_EfficiencyttZ_Zll_;
	    evList[k].weight_ttZ_Zll_ += evList[k].integralttZ_Zll_[perm];
	    
	  }
	  
	}	
	          
      } // end perm loop 
    
    }





    //Missing jet
    else if((integration->integration_type_ == 1 && integration->run_missing_jet_integration_) || (integration->integration_type_ == 0 && integration->force_missing_jet_integration_)){
      
      //Loop over permutations
      for(int perm=0; perm<integration->nbrOfPermut_; perm++){
	
	integration->initVersors_miss(perm);          
	
	//
	// Integrate ttH
	//

	if( integration->include_perm_ttH_[perm] ){
	  
	  integration->signalME_ = true;
	  integration->m_TauTau_2_ = pow(integration->mTauTau_ttH_[perm],2);
	  integration->lowerValues_[PTauLep_id] = integration->PTauLep_ttH_Lower_[perm];
	  integration->upperValues_[PTauLep_id] = integration->PTauLep_ttH_Upper_[perm];
	  integration->lowerValues_[cosThetaTauLepTauHad_id] = integration->cosTheta_diTau_ttH_Lower_[perm];
	  integration->upperValues_[cosThetaTauLepTauHad_id] = integration->cosTheta_diTau_ttH_Upper_[perm];
	  integration->lowerValues_[EQuark1_id] = integration->EQuark1_Lower_[perm];
	  integration->upperValues_[EQuark1_id] = integration->EQuark1_Upper_[perm];
	  integration->lowerValues_[cosThetaNu_tlep_id] = integration->cosThetaNu_tlep_Boundaries_[0];
	  integration->upperValues_[cosThetaNu_tlep_id] = integration->cosThetaNu_tlep_Boundaries_[1];
	  integration->lowerValues_[phiNu_tlep_id] = integration->phiNu_tlep_Boundaries_[0];
	  integration->upperValues_[phiNu_tlep_id] = integration->phiNu_tlep_Boundaries_[1];
	  integration->lowerValues_[cosTheta_missing_jet_id] = integration->cosTheta_missing_jet_Boundaries_[0];
	  integration->upperValues_[cosTheta_missing_jet_id] = integration->cosTheta_missing_jet_Boundaries_[1];
	  integration->lowerValues_[phi_missing_jet_id] = integration->phi_missing_jet_Boundaries_[0];
	  integration->upperValues_[phi_missing_jet_id] = integration->phi_missing_jet_Boundaries_[1];	 	 
	  
	  integration->tot_DrawsttH_ = 0;
	  integration->integr_EfficiencyttH_ = 0;
	  
	  // Compute Integral          
	  if ( integration->flagSameRNG_ )
	  // Copy the initial state to the current RNG State
	    gsl_rng_memcpy ( rng, saveRNGSeed );
	  
	  t = clock();
	  
	  gslIntegrate( &fctDescr_ttH_miss, integration, integration->nbrOfDim_ttH_miss_, integration->nbrOfPoints_ttH_miss_, rng, state_ttH_miss, true,
			&evList[k].integralttH_[perm], 
			&evList[k].stderrttH_[perm],
			&evList[k].chiSquarettH_[perm]);     
	  
	  t = clock() - t;	 

	  evList[k].compTimettH_[perm] =  ((float)t)/CLOCKS_PER_SEC;
	  evList[k].totalDrawsttH_[perm] = integration->tot_DrawsttH_;
	  evList[k].integrationEfficiencyttH_[perm] = integration->integr_EfficiencyttH_;
	  evList[k].weight_ttH_ += evList[k].integralttH_[perm]; 
	  
	}


	
	//
	// Integrate ttZ
	//
	
	if(integration->runttZ_integration_){
	  
	  if( integration->include_perm_ttZ_[perm] ){
	    integration->signalME_ = false;
	    integration->m_TauTau_2_ = pow(integration->mTauTau_ttZ_[perm],2);
	    integration->lowerValues_[PTauLep_id] = integration->PTauLep_ttZ_Lower_[perm];
	    integration->upperValues_[PTauLep_id] = integration->PTauLep_ttZ_Upper_[perm];
	    integration->lowerValues_[cosThetaTauLepTauHad_id] = integration->cosTheta_diTau_ttZ_Lower_[perm];
	    integration->upperValues_[cosThetaTauLepTauHad_id] = integration->cosTheta_diTau_ttZ_Upper_[perm];
	    integration->lowerValues_[EQuark1_id] = integration->EQuark1_Lower_[perm];
	    integration->upperValues_[EQuark1_id] = integration->EQuark1_Upper_[perm];
	    integration->lowerValues_[cosThetaNu_tlep_id] = integration->cosThetaNu_tlep_Boundaries_[0];
	    integration->upperValues_[cosThetaNu_tlep_id] = integration->cosThetaNu_tlep_Boundaries_[1];
	    integration->lowerValues_[phiNu_tlep_id] = integration->phiNu_tlep_Boundaries_[0];
	    integration->upperValues_[phiNu_tlep_id] = integration->phiNu_tlep_Boundaries_[1];
	    integration->lowerValues_[cosTheta_missing_jet_id] = integration->cosTheta_missing_jet_Boundaries_[0];
	    integration->upperValues_[cosTheta_missing_jet_id] = integration->cosTheta_missing_jet_Boundaries_[1];
	    integration->lowerValues_[phi_missing_jet_id] = integration->phi_missing_jet_Boundaries_[0];
	    integration->upperValues_[phi_missing_jet_id] = integration->phi_missing_jet_Boundaries_[1];
	    
	    integration->tot_DrawsttZ_ = 0;
	    integration->integr_EfficiencyttZ_ = 0;
	    
	    // Compute Integral          
	    if ( integration->flagSameRNG_ )
	      // Copy the initial state to the current RNG State
	      gsl_rng_memcpy ( rng, saveRNGSeed );
	    
	    t = clock();	    
	  
	    gslIntegrate( &fctDescr_ttZ_miss, integration, integration->nbrOfDim_ttZ_miss_, integration->nbrOfPoints_ttZ_miss_, rng, state_ttZ_miss, true,
			  &evList[k].integralttZ_[perm], 
			  &evList[k].stderrttZ_[perm],
			  &evList[k].chiSquarettZ_[perm]);     
	    
	    t = clock() - t;	    

	    evList[k].compTimettZ_[perm] =  ((float)t)/CLOCKS_PER_SEC;
	    evList[k].totalDrawsttZ_[perm] = integration->tot_DrawsttZ_;
	    evList[k].integrationEfficiencyttZ_[perm] = integration->integr_EfficiencyttZ_;
	    evList[k].weight_ttZ_ += evList[k].integralttZ_[perm];
	    
	  }
	  
	}
	

	

	//
	// Integrate ttZ Zll
	//
	
	if(integration->runttZ_Zll_integration_){
	  
	  if( integration->include_perm_ttZ_Zll_[perm] ){
	    
	    integration->lowerValues_[EQuark1_ttZ_Zll_id] = integration->EQuark1_Lower_[perm];
	    integration->upperValues_[EQuark1_ttZ_Zll_id] = integration->EQuark1_Upper_[perm];
	    integration->lowerValues_[cosThetaNu_tlep_ttZ_Zll_id] = integration->cosThetaNu_tlep_Boundaries_[0];
	    integration->upperValues_[cosThetaNu_tlep_ttZ_Zll_id] = integration->cosThetaNu_tlep_Boundaries_[1];
	    integration->lowerValues_[phiNu_tlep_ttZ_Zll_id] = integration->phiNu_tlep_Boundaries_[0];
	    integration->upperValues_[phiNu_tlep_ttZ_Zll_id] = integration->phiNu_tlep_Boundaries_[1];
	    integration->lowerValues_[cosTheta_missing_jet_ttZ_Zll_id] = integration->cosTheta_missing_jet_Boundaries_[0];
	    integration->upperValues_[cosTheta_missing_jet_ttZ_Zll_id] = integration->cosTheta_missing_jet_Boundaries_[1];
	    integration->lowerValues_[phi_missing_jet_ttZ_Zll_id] = integration->phi_missing_jet_Boundaries_[0];
	    integration->upperValues_[phi_missing_jet_ttZ_Zll_id] = integration->phi_missing_jet_Boundaries_[1];
	    
	    integration->tot_DrawsttZ_Zll_ = 0;
	    integration->integr_EfficiencyttZ_Zll_ = 0;
	    
	    // Compute Integral          
	    if ( integration->flagSameRNG_ )
	      // Copy the initial state to the current RNG State
	      gsl_rng_memcpy ( rng, saveRNGSeed );
	  
	    t = clock();	    
	    
	    gslIntegrate( &fctDescr_ttZ_Zll_miss, integration, integration->nbrOfDim_ttZ_Zll_miss_, integration->nbrOfPoints_ttZ_Zll_miss_, rng, state_ttZ_Zll_miss, true,
			  &evList[k].integralttZ_Zll_[perm], 
			  &evList[k].stderrttZ_Zll_[perm],
			  &evList[k].chiSquarettZ_Zll_[perm]);     

	    t = clock() - t;	  

	    evList[k].compTimettZ_Zll_[perm] =  ((float)t)/CLOCKS_PER_SEC;
	    evList[k].totalDrawsttZ_Zll_[perm] = integration->tot_DrawsttZ_Zll_;
	    evList[k].integrationEfficiencyttZ_Zll_[perm] = integration->integr_EfficiencyttZ_Zll_;
	    evList[k].weight_ttZ_Zll_ += evList[k].integralttZ_Zll_[perm];
	    
	  }
	  
	}
	
      } //End permutation  
      
    
    }


  }

  logStream_.close();

}


void ThreadScheduler::completeNodeScheduler( eventList_t evList ) {

}


void ThreadScheduler::terminateNodeScheduler( ) {
  
}
