/* 
 * File:   deviceScheduler.h
 * Author: grasseau
 *
 * Created on 27 août 2014, 23:25
 */

#ifndef THREADSCHEDULER_H
#define	THREADSCHEDULER_H

#include <time.h>

# include "ttH_Htautau_MEM_Analysis/MEMAlgo/interface/NodeScheduler.h"
# include "ttH_Htautau_MEM_Analysis/MEMAlgo/interface/MGIntegration.h"

# if (USE_GSL_LIB == 1)
# include "gsl/gsl_math.h"
# include "gsl/gsl_monte_vegas.h"
# include "gsl/gsl_pow_int.h"
# endif

class ThreadScheduler: public NodeScheduler {

 // GG XXX
 public:
 
  clock_t t;

  int totalNbrOfQueueAndProcess;

  long int remainingTasks;
  
  MGIntegration *integration;
  
  // GSL
  gsl_monte_function fctDescr_ttH;
  gsl_monte_vegas_state *state_ttH;

  gsl_monte_function fctDescr_ttZ;
  gsl_monte_vegas_state *state_ttZ;
  gsl_monte_function fctDescr_ttZ_Zll;
  gsl_monte_vegas_state *state_ttZ_Zll;

  gsl_monte_function fctDescr_ttbar_DL_fakelep_tlep;
  gsl_monte_vegas_state *state_ttbar_DL_fakelep_tlep;
  gsl_monte_function fctDescr_ttbar_DL_fakelep_ttau;
  gsl_monte_vegas_state *state_ttbar_DL_fakelep_ttau;


  gsl_rng *rng;
  gsl_rng *saveRNGSeed;
  
  ThreadScheduler();
  ~ThreadScheduler();
 
  virtual void initNodeScheduler( const RunConfig *config, int mpi_rank );
  
  void runNodeScheduler ( eventList_t evMsgList, int nbrOfEvents);

  void completeNodeScheduler( eventList_t evList );

  void terminateNodeScheduler( ); 


};

#ifdef	__cplusplus
extern "C" {
#endif
  

#ifdef	__cplusplus
}
#endif

#endif	/* THREADSCHEDULER_H */

