/* 
 * File:   deviceScheduler.h
 * Author: grasseau
 *
 * Created on 27 août 2014, 23:25
 */

#ifndef NODESCHEDULER_H
#define	NODESCHEDULER_H

# include <fstream>
# include <stdexcept>

# include "ttH_Htautau_MEM_Analysis/MEMAlgo/interface/config.h"
# include "ttH_Htautau_MEM_Analysis/MEMAlgo/interface/RunConfig.h"
# include "ttH_Htautau_MEM_Analysis/MEMAlgo/interface/MGIntegration.h"


const int nbrOfEventsPerBlock = NbrEventPerBlock;
typedef IntegrationMsg_t eventList_t[nbrOfEventsPerBlock];

// Queue Status
typedef struct {
  
  const int ComputationCompleted = 1;
  const int ComputationStarted   = 0;
  const int ComputationNotActive = -1;

  int     platform;
  int     device; 
  // Index of the processed events (startEventIndex, nbrOfEvents) in 
  // the event list vegasState
  int startEventIndex;
  int nbrOfEvents;
  int     status;
} QueueStatus_t;

class NodeScheduler {


 // GG XXX
 public:
  typedef enum {Vegas, EvalOnGrid} ComputationType_t;
  struct timespec delay;
  
  const RunConfig *config;

  int nbrOfQueues;
  // 
  // Devices Status
  //
  QueueStatus_t *queueStatus;
  
  // For Grid evaluation
  unsigned long int nbrPointsPerAxe;
  ComputationType_t computationType;
  
  // One log file per process
  ofstream  logStream_; 
  
  NodeScheduler() { nbrPointsPerAxe = 2; computationType = Vegas; 
                    delay.tv_sec = 0L;  //.tv_sec
                    delay.tv_sec = 50L; // tv_nsec 
                  };
                  
  virtual ~NodeScheduler() {};

  virtual void initNodeScheduler( const RunConfig *config, int mpi_rank );

  virtual void runNodeScheduler ( eventList_t evMsgList, int nbrOfEvents)= 0;

  virtual void completeNodeScheduler( eventList_t evList ) = 0;

  virtual void terminateNodeScheduler( ) = 0;


  // Integration selection
  // EvalOnGrid
  ComputationType_t getComputationType() { return computationType;}
  void setComputationType(ComputationType_t ct) { computationType=ct;}

};

#ifdef	__cplusplus
extern "C" {
#endif

#ifdef	__cplusplus
}
#endif

#endif	/* NODESCHEDULER_H */

