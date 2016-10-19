/* 
 * File:   Vector3.h
 * Author: grasseau (grasseau@llr.in2p3.fr)
 *
 * Created on 3 juin 2015, 18:41
 */

/*
 * ROOT vector arithmetics for OCL kernels
 */

#ifndef VECTOR3_H
#define	VECTOR3_H

# include "ttH_Htautau_MEM_Analysis/MEMAlgo/interface/config.h"

typedef struct {
  double x;
  double y;
  double z;
} Vector3_t;

Static inline Vector3_t
Vector3_constr( double x, double y, double z ) {
  Vector3_t v;
  v.x = x; v.y = y; v.z = z;
  return v;
}

Static inline void
Vector3_set_Vector3( Vector3_t *p, Const Vector3_t *src ) {
  p->x = src->x;  p->y = src->y;  p->z = src->z;
}

Static inline void
Vector3_copyTo_Vector3( Const Vector3_t *p, Vector3_t *dst ) {
  dst->x = p->x; dst->y = p->y; dst->z = p->z;
}

Static inline double
Vector3_mag2( Const Vector3_t *p ) {
  double  mod2 = p->x*p->x
              + p->y*p->y
              + p->z*p->z;
 
  return mod2;
}

Static inline double
Vector3_mag( Const Vector3_t *p ) {
  return sqrt( Vector3_mag2( p ) );
}

Static inline Vector3_t
Vector3_unit( Const Vector3_t *p ) {
  double  c = Vector3_mag2( p );
  c = (c != 0) ? 1.0 / sqrt( c) : 1.0; 
  Vector3_t unit = { p->x*c, p->y*c, p->z*c};
  return unit;
}

Static inline Vector3_t
Vector3_mult_double( Const Vector3_t *p, double c ) {
  Vector3_t vect = { p->x*c, p->y*c, p->z*c};
  return vect;
}

Static inline double 
Vector3_perp2(Const Vector3_t *p) {
  return p->x*p->x + p->y*p->y;
}

Static inline double 
Vector3_perp(Const Vector3_t *p) {
  return sqrt( Vector3_perp2( p ) );
}

Static inline double 
Vector3_dot(Const Vector3_t *p, Const Vector3_t *q) {
  return 
    p->x*q->x + p->y*q->y + p->z*q->z;
}

Static inline Vector3_t
Vector3_cross(Const Vector3_t *p, Const Vector3_t *q) {
  Vector3_t u = Vector3_constr( p->y*q->z - q->y*p->z, 
                 p->z*q->x - q->z*p->x, 
                 p->x*q->y - q->x*p->y
               );

  return u;
}

Static inline double
Vector3_angle( Const Vector3_t *p, Const Vector3_t *q) {
   double ptot2 = Vector3_mag2( p ) * Vector3_mag2( q );
   if(ptot2 <= 0) {
      return 0.0;
   } else {
      double arg = Vector3_dot(p, q) / sqrt(ptot2);
      if(arg >  1.0) arg =  1.0;
      if(arg < -1.0) arg = -1.0;
      return acos(arg);
   }
}

Static inline Vector3_t
Vector3_sum_3Vector3( 
    double a, Const Vector3_t *va, 
    double b, Const Vector3_t *vb, 
    double c, Const Vector3_t *vc ) { 
  Vector3_t s = {
    a*va->x + b*vb->x + c*vc->x,
    a*va->y + b*vb->y + c*vc->y,
    a*va->z + b*vb->z + c*vc->z
  };
  return s;
}

Static inline double
Vector3_cosTheta( Const Vector3_t *v ) {
   double ptot = Vector3_mag( v );
   return ptot == 0.0 ? 1.0 : v->z / ptot;
}

Static inline double
Vector3_pseudoRapidity( Const Vector3_t *v) {
   //Double_t m = Mag();
   //return 0.5*log( (m+fZ)/(m-fZ) );
   // guard against Pt=0
   double cosTheta = Vector3_cosTheta( v );
   if (cosTheta*cosTheta < 1) 
     return -0.5* log( (1.0-cosTheta)/(1.0+cosTheta) );
// GG XXX  Warning("PseudoRapidity","transvers momentum = 0! return +/- 10e10");
   if (v->z > 0) return  10e10;
   else          return -10e10;
}

Static inline double
Vector3_phi( Const Vector3_t *v ) {
   //return the  azimuth angle. returns phi from -pi to pi
   return (v->x == 0.0 && v->y == 0.0) ? 0.0 : atan2(v->y,v->x);
}

Static inline double
Vector3_eta( Const Vector3_t *v) {
  return Vector3_pseudoRapidity(v);
}

#endif	/* VECTOR3_H */

