/* 
 * File:   Constants.h
 * Author: grasseau
 *
 * Created on 6 mars 2015, 19:47
 */

#ifndef CONSTANTS_H
#define	CONSTANTS_H

namespace EventDescr {
  const int electronType = 0;
  const int muonType = 1;
}; 

namespace Physics {
  const double mtop   = 173.21;
  const double mHiggs = 125.09;
  const double mTau   = 1.776;
  const double mW     = 80.385;
  const double mZ     = 91.187;
  const double mb     = 4.7;
  const double mTau2 = mTau*mTau;
  const double mTauTauH2 = mHiggs*mHiggs;
  const double mTauTauZ2 = mZ*mZ;
  const double mElectron = 0.000511; // in Gev
  const double mMuon     = 0.105;    // in Gev
}
#endif	/* CONSTANTS_H */

