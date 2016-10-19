#ifndef MEM_MEMAlgo_Utility_H
#define MEM_MEMAlgo_Utility_H

# if defined(__OPENCL_VERSION__)
#   pragma OPENCL EXTENSION cl_khr_fp64 : enable
# else
#   include <math.h>
# endif

# define False 0
# define True 1

# ifndef M_PI
# define M_PI 3.14159265358979323846
# endif

typedef char MG_bool_t;


# define Min(a, b) ((a) < (b) ? (a) : (b))
# define Max(a, b) ((a) > (b) ? (a) : (b))

# define Sgn( a, b) ( ((b)<0) ? -fabs(a) : fabs(a) )

# define Modulus( x, y, z) ( sqrt( (x)*(x) + (y)*(y) + (z)*(z) ) )

# define SetArray( type, value, start_dst, end_dst )  for( type *ptr_dst = (start_dst), *ptr_end = (end_dst); ptr_dst != ptr_end; *ptr_dst++ = value)
# define CopyArray( type, start_src, start_dst, end_dst )  for( type *ptr_dst = (start_dst), *ptr_end = (end_dst), *ptr_src=(start_src); ptr_dst != ptr_end; *ptr_dst++ = *ptr_src++)
# define CopyArrayT( type_src, type_dest, start_src, start_dst, end_dst ) { type_src *ptr_src = (start_src); \
   for( type_dest *ptr_dst = (start_dst), *ptr_end = (end_dst); ptr_dst != ptr_end; *ptr_dst++ = *ptr_src++);}

# endif // Utility_H
