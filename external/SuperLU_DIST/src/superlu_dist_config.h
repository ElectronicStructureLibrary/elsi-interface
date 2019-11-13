/* superlu_dist_config.h.in */

#define XSDK_INDEX_SIZE 32
#define HAVE_PARMETIS TRUE
#define SLU_HAVE_LAPACK TRUE

#if (XSDK_INDEX_SIZE == 64)
#define _LONGINT 1

#endif
