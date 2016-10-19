/* 
 * File:   config.h
 * Author: grasseau
 *
 * Created on 4 mai 2015, 15:12
 */

#ifndef MEM_MEMAlgo_CONFIG_H
#define	MEM_MEMAlgo_CONFIG_H

#ifdef	__cplusplus
extern "C" {
#endif

// Use ROOT environment
# define USE_ROOT 1

// Use GSL library
# define USE_GSL_LIB 1
//
//  Kernel Debug / Tracer
//
// KERNEL_PRINTF 
// + 0  No print
// + 2  Inside box loop
// + 1  Initialize/Finalize Chi-Squared
// + 10 All
# define LEVELNoPrint       0
# define LEVELChiSqLoop     1
# define LEVELBoxLoop       2
# define LEVELAll          10
// Default
// KERNEL_PRINTF flag
# define USE_KERNEL_PRINTF 0 // 0/1
# if defined(cl_amd_printf)
#   pragma OPENCL EXTENSION cl_amd_printf :enable
#   define KERNEL_PRINTF 0 
# elif defined(cl_intel_printf)
#   pragma OPENCL EXTENSION cl_intel_printf : enable
#   define KERNEL_PRINTF USE_KERNEL_PRINTF
# elif ! defined(__OPENCL_VERSION__)
//#   define KERNEL_PRINTF USE_KERNEL_PRINTF
#   define KERNEL_PRINTF 1
# else 
#   define KERNEL_PRINTF 0
# endif


// Used for check sums: sum_fval, sum_nbrCalls
# define TRACER 0   // 0/1

//
// Kernel definitions
//
# if defined(__OPENCL_VERSION__)
//
// For OCL compiler
//
#  define cl_char   char
#  define cl_short  short
#  define cl_int    int
#  define cl_uint   unsigned int
#  define cl_long   long
#  define cl_double double
#  define cl_ulong  unsigned long
#  define Static
// Pb with const <-> __constant
#  define Const
#  pragma OPENCL EXTENSION cl_khr_fp64: enable
#  define PowInt pown
# else
//
// For classical compiler
//  
// In the OCL API file
#  include <string.h>
#  include <stdio.h>
#  include <math.h>
#  include <stdlib.h>
#  include <stdbool.h>
// #  include <CL/cl.h>
//#  include "Analysis/MEMAlgo/interface/KernelEmu.h"
#  define __kernel 
#  define __global 
#  define __local
#  define __private
#  define __constant 
#  define Static static
#  define Const  const
#  define PowInt pow
# endif

// Computig modes
// GG XXX Future uses
# define VegasOneEventMultiDevices    1
# define VegasMultiEventsMultiDevices 2

//
//   VEGAS Configuration
//
# define MaxNumberOfDIMs 8
# define MaxNumberOfBins 50          // Must be even 
# define NbrMaxVegasChi2Iterations 5
// GG XXX Not used
# define VEGASWarmUpNbrOfPoints 10000
// GG XXX should be maximum value
// Inv ??? #  define VEGASNbrOfPoints       64
# define MAXNbrOfVegasBoxes 10000
//# define VEGASNbrOfPoints     20002
//# define VEGASNbrOfPoints       65536 
//# define VEGASNbrOfPoints       586
//# define VEGASNbrOfPoints       2048   
// Used for reduction processes:
//   NbreGroups < DEVICEMaxNbrOfWorkGroups
# define DEVICEMaxNbrOfWorkGroups                 2048
// Max. number of variables to reduce
# define DEVICEMaxNbrVariablesForWorkGroupStorage 6
# define DEVICEMaxDebugBuffer 1024
//
// Chi-squared computation
//
# define USE_ORIGINAL_CHISQ_FORMULA 0
//
// RNG
//
# define RAND48
//# define RANDCYCLE
//
// GG XXX Only valid in Kernels
# define RESTART_SEED_AT_IITERATION 0

//
// Function to integrate
//
# define VBFFCT
//# define SinUOverUFCT
//# define PyramidFCT
# if defined(TestCubeFCT)
#   define FCTDims 3
#   define XLowerDomain 0
#   define XUpperDomain M_PI
#   define OCL_EvaluateFct( x, dims, parameters) testCube( (x), (dims), (parameters) )
#   define VEGASFCTName testCube
# elif defined(SinUOverUFCT)
#   define XLowerDomain 0
#   define XUpperDomain (2*M_PI)
#   define FCTDims 5
#   define OCL_EvaluateFct( x, dims, parameters) sinu_over_u( (x), (dims), (parameters) )
#   define VEGASFCTName sinu_over_u
# elif defined(PyramidFCT)
#   define XLowerDomain 0
#   define XUpperDomain (1.0)
#   define FCTDims 5
#   define OCL_EvaluateFct( x, dims, parameters) pyramid( (x), (dims), (parameters) )
#   define VEGASFCTName pyramid
# elif defined(VBFFCT)
#   define XLowerDomain 0.9
#   define XUpperDomain 1.1
#   define FCTDims 5
#   define OCL_EvaluateFct( x, dims, parameters) vbf( (x), (dims), (parameters) )
#   define VEGASFCTName vbf
//  For finding kernel sources
// ??? #   define USE_VBFKernels
// # else
// # error   
# endif

// MPI  
# define MPISchedulerEnabled 1  // 0/1 Disable/Enable
# define NbrEventPerBlock  16


// GG XXX obsolete options
# if 0
#  define PlatformToForce 0
#  define ForceToDevice   1
# endif
//
# define ForceToDeviceONE 0
# define MPIForceToDevice  0 // 0/1
# if 0
# define SplitCPUs 16
# endif      

//
// Debugging
//
// Setting BlockingOrders to 1, 
//  all command/orders to the command queue are executed 
//  in blocking order mode(oclWaitForEvents()) 
// Setting BlockingOrders to 0, non-blocking mode  
# define BlockingOrders 0

// GG XXX Not used
# define InlineExtern
# define InlineStatic extern

#ifdef	__cplusplus
}
#endif
          
#endif	/* CONFIG_H */

