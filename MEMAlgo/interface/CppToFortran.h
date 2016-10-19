/* 
 * File:   LHAPDF-interface.h
 * Author: grasseau
 *
 * Created on 28 mai 2015, 13:11
 */

#ifndef MEM_MEMAlgo_LHAPDF_INTERFACE_H
#define	MEM_MEMAlgo_LHAPDF_INTERFACE_H

# include "lhapdf.h"

#ifdef	__cplusplus
extern "C" {
#endif

void pdf_initpdfsetbyname_(char *name, int *len);
void numberpdf_(int *N);
void initpdf_(int *i);
double alphaspdf_(double *QMZ);
void evolvepdf_(double *x, double *Q, double *f);
void evolvepdfm_(int *iset, double *x, double *Q, double *f);
void getdesc_();
void getdescm_( int *nset);
void pdf_initpdfsetbyname_(char *name, int *len);


#ifdef	__cplusplus
}
#endif

#endif	/* LHAPDF_INTERFACE_H */

