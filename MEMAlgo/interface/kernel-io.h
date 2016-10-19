/* 
 * File:   io.h
 * Author: grasseau
 *
 * Created on 4 juin 2014, 12:32
 */

#ifndef MEM_MEMAlgo_IO_H
#define	MEM_MEMAlgo_IO_H

# include "ttH_Htautau_MEM_Analysis/MEMAlgo/interface/config.h"

// # define OCL_IO_BUFFER_SIZE 1048576
# define OCL_IO_BUFFER_SIZE 524288
# define OCL_IO_WARNING_SIZE (OCL_IO_BUFFER_SIZE - 256 )

# ifdef __cplusplus
  extern "C" {
# endif

typedef struct __attribute__((packed)) OCL_FILE_s {
  cl_char buffer[OCL_IO_BUFFER_SIZE];
  cl_long pos;    
  // GG XXX __attribute__((aligned(8)));
} OCL_FILE;

# ifdef __cplusplus
  }
# endif

#endif	/* IO_H */

