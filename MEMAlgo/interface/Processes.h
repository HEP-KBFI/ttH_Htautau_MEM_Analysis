/* 
 * File:   Processes.h
 * Author: grasseau
 *
 * Created on 1 décembre 2014, 19:15
 */

#ifndef PROCESSES_H
#define	PROCESSES_H


# include <vector>
using namespace std;

//static const int PDFSet = 1;
static const int PDFSet = 0;

double getME_ttH_gg( vector<double*> &p );
double getME_ttH_qqbar( vector<double*> &p );

double getME_ttZ_Zonly_gg( vector<double*> &p );
double getME_ttZ_Zonly_uubar( vector<double*> &p );

double getME_ttZ_Zonly_Zll_gg( vector<double*> &p );
double getME_ttZ_Zonly_Zll_uubar( vector<double*> &p );

double getME_gg_ttbar( vector<double*> &p );
double getME_qqbar_ttbar( vector<double*> &p );

// ME
void init_ttH_HTauTauMEProcesses( const char *configname );

double get_ttHWeightedME( double x0, double x1, const vector<double*> &p );
double get_ttZZonlyWeightedME( double x0, double x1, const vector<double*> &p );
double get_ttZZonlyZllWeightedME( double x0, double x1, const vector<double*> &p );
double get_ttbarWeightedME( double x0, double x1, const vector<double*> &p );

#endif	/* PROCESSES_H */

