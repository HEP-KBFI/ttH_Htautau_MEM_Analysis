/* 
 * File:   LorentzVector.h
 * Author: grasseau
 *
 * Created on 3 juin 2015, 18:13
 */

/*
 * ROOT Lorentz vector arithmetics for OCL kernels
 */

#ifndef LORENTZVECTOR_H
#define	LORENTZVECTOR_H

# include "ttH_Htautau_MEM_Analysis/MEMAlgo/interface/config.h"
# include "ttH_Htautau_MEM_Analysis/MEMAlgo/interface/Vector3.h"

# if (KERNEL_PRINTF != 0 )
# define LorentzVector_print( title, lv ) ( \
  printf("%s: %e %e %e %e\n", title, lv.P.x, lv.P.y, lv.P.z, lv.E ) )
# else
# define LorentzVector_print( title, lv ) {double ___a = lv.E;};
# endif

typedef struct {
  Vector3_t P;
  double E;
} LorentzVector_t;

Static inline LorentzVector_t
LorentzVector_constr_global( __global LorentzVector_t *v) {
  LorentzVector_t lv;  
  lv.P.x = v->P.x;  lv.P.y = v->P.y;  lv.P.z = v->P.z; lv.E = v->E;
  return lv;
}

Static inline LorentzVector_t
LorentzVector_constr( Vector3_t *v, double E ) {
  LorentzVector_t lv;  
  lv.P.x = v->x;  lv.P.y = v->y;  lv.P.z = v->z; lv.E = E;
  return lv;
}

Static inline LorentzVector_t
LorentzVector_constr_double( double Px, double Py, double Pz, double E ) {
  LorentzVector_t lv;
  lv.P.x = Px;  lv.P.y = Py;  lv.P.z = Pz; lv.E = E;
  return lv;
}

Static inline void
LorentzVector_set( LorentzVector_t *p, 
        double x, double y, double z, double E) {
  p->P.x = x;  p->P.y = y;  p->P.z = z; p->E = E;
}  

Static inline void
LorentzVector_set_DoubleAddr( LorentzVector_t *p, 
        const double *q) {
  p->P.x = q[0];  p->P.y = q[1];  p->P.z = q[2]; p->E = q[3];
} 

Static inline void
LorentzVector_set_LorentzVector( LorentzVector_t *p, Const LorentzVector_t *src ) {
  Vector3_set_Vector3( &p->P, &src->P );
  p->E = src->E;
}

Static inline void
LorentzVector_copyTo_Vector3( Const LorentzVector_t *p, Vector3_t *dst ) {
  Vector3_copyTo_Vector3( &p->P, dst);
}

Static inline Vector3_t
LorentzVector_getVector3( Const LorentzVector_t *p ) {
  return p->P;
}

Static inline double
LorentzVector_getE( Const LorentzVector_t *p ) {
  return p->E;
}

Static inline double
LorentzVector_getP( Const LorentzVector_t *p ) {
  return sqrt( Vector3_mag2( &p->P ) );
}

Static inline double
LorentzVector_phi( Const LorentzVector_t *lv ) {
  //return the  azimuth angle. returns phi from -pi to pi
  return Vector3_phi( &lv->P );
}

Static inline double
LorentzVector_eta( Const LorentzVector_t *p ) {
  return Vector3_eta( &p->P );
}

Static inline void
LorentzVector_setPtEtaPhiE( LorentzVector_t *lv, double pt, double eta, double phi, double e) {
   pt = fabs( pt );
   lv->P.x = pt*cos(phi);
   lv->P.y = pt*sin(phi);
   lv->P.z = pt*sinh(eta);
   lv->E = e;
}

Static inline LorentzVector_t
LorentzVector_sum_4( LorentzVector_t *p1, LorentzVector_t *p2, 
        LorentzVector_t *p3, LorentzVector_t *p4 ) {

  LorentzVector_t sum;
  sum.P.x = p1->P.x + p2->P.x + p3->P.x + p4->P.x;
  sum.P.y = p1->P.y + p2->P.y + p3->P.y + p4->P.y;
  sum.P.z = p1->P.z + p2->P.z + p3->P.z + p4->P.z;
  sum.E = p1->E + p2->E + p3->E + p4->E;
  return sum;  
}

Static inline LorentzVector_t
LorentzVector_sum_5( LorentzVector_t *p1, LorentzVector_t *p2, 
        LorentzVector_t *p3, LorentzVector_t *p4, LorentzVector_t *p5) {

  LorentzVector_t sum;
  sum.P.x = p1->P.x + p2->P.x + p3->P.x + p4->P.x + p5->P.x;
  sum.P.y = p1->P.y + p2->P.y + p3->P.y + p4->P.y + p5->P.y;
  sum.P.z = p1->P.z + p2->P.z + p3->P.z + p4->P.z + p5->P.z;
  sum.E   = p1->E + p2->E + p3->E + p4->E + p5->E;
  return sum;  
}

Static inline double
LorentzVector_Pt( Const LorentzVector_t *p) {
   return Vector3_perp( &p->P );
}

Static inline Vector3_t
LorentzVector_getBoostVector( Const LorentzVector_t *p) {
  Vector3_t v;
  double E = p->E;
  v.x = p->P.x / E ;  v.y = p->P.y / E ;   v.z = p->P.z / E ;
  return v;
}

Static inline void 
LorentzVector_boost(LorentzVector_t *lv, Vector3_t *v) {
   //Boost this Lorentz vector
   double b2 = v->x*v->x + v->y*v->y + v->z*v->z;
   double gamma = 1.0 / sqrt(1.0 - b2);
   double bp = v->x*lv->P.x + v->y*lv->P.y + v->z*lv->P.z;
   double gamma2 = b2 > 0 ? (gamma - 1.0)/b2 : 0.0;
   /*
   SetX(X() + gamma2*bp*bx + gamma*bx*T());
   SetY(Y() + gamma2*bp*by + gamma*by*T());
   SetZ(Z() + gamma2*bp*bz + gamma*bz*T());
   */
   double c = gamma2*bp + gamma*lv->E;
   lv->P.x = lv->P.x + c * v->x;
   lv->P.y = lv->P.y + c * v->y;
   lv->P.z = lv->P.z + c * v->z;
   
   //SetT(gamma*(T() + bp));
   lv->E = gamma * (lv->E + bp);
}

#endif	/* LORENTZVECTOR_H */
