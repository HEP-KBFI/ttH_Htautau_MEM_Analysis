/* 
 * File:   MGIntegration-def.h
 * Author: grasseau
 *
 * Created on 4 juin 2015, 17:17
 */

#ifndef MGINTEGRATION_DEF_H
#define	MGINTEGRATION_DEF_H

// Argument index of the VEGAS integration 
# define mTauTau2_id       0 //squared invariant mass of the di-tau system
# define PTauLep_id        1 //module of the momentum of the tau
# define PTauHad_id        2 //module of the momentum of the anti-tau
# define cosThetaTauLepTauHad_id 2 //cos of the TauLep/TauHad angle
# define POutQuark1_id     3 //module of the momentum of the final state Quark1
# define POutQuark2_id     4 //module of the momentum of the final state Quark2
# define cosThetaTauLep_id 5 //cosine between the lepton and the associated tau
# define phiTauLep_id      6 //phi angle between the lepton and the associated tau

// Verbose values
# define ResultLevel      1  // Display integration results
# define IntegrationLevel 2  // Display integration parameters
# define IntegrandLevel   3  // Display Integrand information

// Dimension max of the integration
# define DimensionMax 7

#endif	/* MGINTEGRATION_DEF_H */

