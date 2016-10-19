/* 
 * File:   lhapdf.h
 * Author: grasseau
 *
 * Created on 26 mars 2014, 15:23
 */

# include "ttH_Htautau_MEM_Analysis/MEMAlgo/interface/config.h"

#ifndef MEM_MEMAlgo_LHAPDF_C_H
#define	MEM_MEMAlgo_LHAPDF_C_H
/* parmsetup.inc
! nopmax is the maximum number of PDF parameters one can use.
! noemax is the maximum number of PDF's in a list. 
! npfmax is the maximum number of functional parameters
! nofmax is the maximum number of functionals
! linemax is the maximum number of lines in pdf description text
! nmxset is the max number of sets that can be initialised at one time
*/
// Index + 1 in Fortran
# define NoeMAX  1001
# define NopMAX  41
# define NpfMAX  11
# define NofMAX  11
# define LineMAX 21
# define NmxSET  4 

# define MXX   205 // Fortran 204
# define MXQ    26 // Fortran  25
# define MXF     6 // Fortran   6
# define MaxVal  3 // Fortran   3
# define NHESS  53 // Fortran  52
// parameter (MXPQX=(MXF+1+MaxVal) * MXQ * MXX)
# define MXPQX  ((MXF+1+MaxVal) * (MXQ-1) * (MXX-1) + 1)
# define NqVEC  4 // Fortran 4

// GG XXX fuse
# define Error( status, a, b) (status)->Error = (a); (status)->errValue = (b) 
# define Warning( status, a, b) (status)->Warning = (a); (status)->warnValue = (b) 
# define Msg(a,b) (a)

# define Square(a) ((a)*(a))

# if (KERNEL_PRINTF != 0 )

# define PRINT2D( str, array, N, iset) \
 printf(" %s [%d, 0..%d] %f %f % f ... %f %f\n",  str, iset, N-1, \
   array[iset][0], array[iset][1], array[iset][2], array[iset][N-2], array[iset][N-1] )

/*
// # define PRINT2D( str, array, N, iset) \
//  printf(" %s [%d, 0..%d] %f %f %f\n",  str, iset, N-1, \
//    array[iset][0], array[iset][1], array[iset][N-1] )
*/

# define PRINT3D( str, array, N, k, iset) \
 printf(" %s [%d, %d, 0..%d] %f %f %f\n",  str, iset, k, N-1, \
   array[iset][k][0], array[iset][k][1], array[iset][k][N-1] )
//   array[0][j][k], array[1][j][k], array[N-1][j][k] )

# define PRINT1D( str, array, N) \
 printf(" %s [%d] %f %f %f\n",  str, 0, \
   array[0], array[1], array[N-1] )

# define PRINT1DInt( str, array, N) \
 printf(" %s [0..%d] %d %d %d %d\n",  str, N-1, \
   array[0], array[1], array[2], array[N-1] )

# define PRINT1DStr( str, array, N) \
 printf(" %s [%d] : \n      %s\n      %s\n      %s\n",  str, 0, \
   array[0], array[1], array[N-1] )

# endif // # if defined(__OPENCL_VERSION__)

#ifdef	__cplusplus
extern "C" {
#endif
  
# define Bad_PDF_set  1
# define X_out_of_range_in_CtLhCtq65Pdf 2
# define Q_out_of_range_in_CtLhCtq65Pdf 3
# define Warning_Iparton_out_of_range_in_CtLhCtq65Pdf 4
# define Severe_error_x_negatif_in_CtLhPartonX65  10
# define Severe_error_x_gt_1_in_CtLhPartonX65 11
# define Warning_Out_of_range_X  12
# define Warning_Out_of_range_Q  13
# define Fatal_CtLhPolint4_call  20
  
typedef struct{
  // COMMON
  struct { 
    //character*16 name(nmxset)
    char name[NmxSET][16];
    int nmem[NmxSET];
    int ndef[NmxSET];
    int mem;
  } NAME;
  // common/SET/iset,iimem
  struct {
    int iset;
    int iimem;
  } SET;
  struct {
    char   lhaparm[21][21];
    double lhavalue[21];
    char   commoninitflag[21];
  } lhacontrol;
   struct {
     int Nx[NmxSET], Nt[NmxSET], NfMx[NmxSET];   
   } CtqPar2;
   struct {
     int Nfl[NmxSET], Iorder[NmxSET];
     double Alambda[NmxSET];
   } QCDtable; 
   struct {
     // GG XV(0:MXX,nmxset), TV(0:MXQ,nmxset), UPD(0:nhess,MXPQX,nmxset)                     
     double Al[NmxSET];
     double XV[NmxSET][MXX];
     double TV[NmxSET][MXQ];
     double UPD[NmxSET][MXPQX][NHESS];
   } CtqPar1nhess65;
   struct {
     double Qini[NmxSET];
     double Qmax[NmxSET];
     double Xmin[NmxSET];
   } XQrange;
   struct {
     double xlast; 
     double qlast;
     int nxsave;
    // Fortran xvpow(0:mxx) 
     double XVpow[MXX+1];
     int Jq;
     int JLq;
     double tdet;
     double tmp1, tmp2;
     double tt;
   } ctq65co;
         
  // SAVE
  struct { 
   int member[NmxSET];
   // save alpha.f
   int iset; // GG XXX conflit avec SET.iset ?
  } ALPHA;
  
  // Save QCDparams
  struct {
  // GG parmXmin(nmxset,0:noemax)
  double parmXmin [NmxSET][NoeMAX];
  double parmXmax [NmxSET][NoeMAX];
  double parmQ2min[NmxSET][NoeMAX];
  double parmQ2max[NmxSET][NoeMAX];
  } QCDparam;
  
  // DATA
  struct {
    double OneP; // GG XXX data =1.00001; 
    //!**** choice of interpolation variable 
    double xpow; // GG XXX data = 0.30; 
    int ixprint,iqprint; // Data = 0, 0
  } Data;
  // Inv int enableWarnMsg;
  int Warning;
  double warnValue;
  int Error;
  double errValue;
} pdf_t;

 // alpha.f
inline int PDF_getnset( __global Const pdf_t *_this) {   
  return _this->ALPHA.iset;
}
        
inline void PDF_setnset( __global pdf_t *_this, int nset) {
  _this->ALPHA.iset = nset;
 }
 
inline int PDF_getnmem( __global Const pdf_t* _this, int nset) {
  return _this->ALPHA.member[nset];
}

inline void PDF_setnmem( __global pdf_t* _this, int nset, int nmem) {
  _this->ALPHA.member[nset] = nmem;
}
// LHpdflib.F
inline __global Const char* PDF_getlhaparm( __global Const pdf_t* _this, int nparm) {
  return &_this->lhacontrol.lhaparm[nparm][0];
}

// QCDparams.f
inline void PDF_getMinMaxM( __global Const pdf_t* _this, int nset,
          double* xmin, double* xmax, double* q2min, double* q2max) {
  int mem = _this->SET.iimem;
  *xmin  = _this->QCDparam.parmXmin [nset][mem];
  *xmax  = _this->QCDparam.parmXmax [nset][mem];
  *q2min = _this->QCDparam.parmQ2min[nset][mem];
  *q2max = _this->QCDparam.parmQ2max[nset][mem];
}

# if ! defined(__OPENCL_VERSION__)
void PDF_fillNAME_( char *name[], int* nmem, int* ndef, int *mem);
pdf_t* PDF_new( );
# endif

# if (KERNEL_PRINTF != 0 )
// # if ! defined(__OPENCL_VERSION__)
void PDF_dump(__global pdf_t *pdf);
# endif

# if ! defined(__OPENCL_VERSION__)
void PDF_copy( Const pdf_t *src, pdf_t *dest );
# endif   // # if defined(__OPENCL_VERSION__)

void PDF_dump( __global pdf_t *pdf);
  
#ifdef	__cplusplus
}
#endif

#endif	/* LHAPDF_C_H */

