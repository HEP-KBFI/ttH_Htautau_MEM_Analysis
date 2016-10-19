//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.2.2, 2014-11-06
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#ifndef MG5_Sigma_sm_ttH_uux_ttxtamtap_H
#define MG5_Sigma_sm_ttH_uux_ttxtamtap_H

#include <complex> 
#include <vector> 

#include "ttH_Htautau_MEM_Analysis/MEMAlgo/interface/Parameters_sm.h"

using namespace std; 

//==========================================================================
// A class for calculating the matrix elements for
// Process: u u~ > t t~ h WEIGHTED=4
// *   Decay: h > ta- ta+ WEIGHTED=2
// Process: c c~ > t t~ h WEIGHTED=4
// *   Decay: h > ta- ta+ WEIGHTED=2
// Process: d d~ > t t~ h WEIGHTED=4
// *   Decay: h > ta- ta+ WEIGHTED=2
// Process: s s~ > t t~ h WEIGHTED=4
// *   Decay: h > ta- ta+ WEIGHTED=2
//--------------------------------------------------------------------------

class CPPProcess_ttH_qqbar
{
  public:

    // Constructor.
    CPPProcess_ttH_qqbar() {}

    // Initialize process.
    virtual void initProc(string param_card_name); 

    // Calculate flavour-independent parts of cross section.
    virtual void sigmaKin(); 

    // Evaluate sigmaHat(sHat).
    virtual double sigmaHat(); 

    // Info on the subprocess.
    virtual string name() const {return "u u~ > t t~ ta- ta+ (sm)";}

    virtual int code() const {return 0;}

    const vector<double> & getMasses() const {return mME;}

    // Get and set momenta for matrix element evaluation
    vector < double * > getMomenta(){return p;}
    void setMomenta(vector < double * > & momenta){p = momenta;}
    void setInitial(int inid1, int inid2){id1 = inid1; id2 = inid2;}

    // Get matrix element vector
    const double * getMatrixElements() const {return matrix_element;}

    // Constants for array limits
    static const int ninitial = 2; 
    static const int nexternal = 6; 
    static const int nprocesses = 2; 

  private:

    // Private functions to calculate the matrix element for all subprocesses
    // Calculate wavefunctions
    void calculate_wavefunctions(const int perm[], const int hel[]); 
    static const int nwavefuncs = 10; 
    std::complex<double> w[nwavefuncs][18]; 
    static const int namplitudes = 2; 
    std::complex<double> amp[namplitudes]; 
    double matrix_uux_ttxh_h_tamtap(); 

    // Store the matrix element value from sigmaKin
    double matrix_element[nprocesses]; 

    // Color flows, used when selecting color
    double * jamp2[nprocesses]; 

    // Pointer to the model parameters
    Parameters_sm * pars; 

    // vector with external particle masses
    vector<double> mME; 

    // vector with momenta (to be changed each event)
    vector < double * > p; 
    // Initial particle ids
    int id1, id2; 

}; 


#endif  // MG5_Sigma_sm_uux_ttxtamtap_H
