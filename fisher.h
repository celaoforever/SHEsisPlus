

#ifndef __FISHER_H__
#define __FISHER_H__

#ifndef R_EXT_MEMORY_H_
#define R_EXT_MEMORY_H_
#ifdef __cplusplus
extern "C" {
#endif

char* vmaxget(void);
void vmaxset(char*);

void R_gc(void);

char* R_alloc(long, int);
char* S_alloc(long, int);
char* S_realloc(char*, long, long, int);

#ifdef __cplusplus
}
#endif

#endif /* R_EXT_MEMORY_H_ */

#ifndef R_EXT_BOOLEAN_H_
#define R_EXT_BOOLEAN_H_

#undef FALSE
#undef TRUE

#ifdef __cplusplus
extern "C" {
#endif
typedef enum { FALSE = 0, TRUE /*, MAYBE */ } Rboolean;

#ifdef __cplusplus
}
#endif

#endif /* R_EXT_BOOLEAN_H_ */

#ifndef R_EXT_CONSTANTS_H_
#define R_EXT_CONSTANTS_H_

#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884197169399375
#endif

#define PI M_PI
#define SINGLE_EPS FLT_EPSILON
#define SINGLE_BASE FLT_RADIX
#define SINGLE_XMIN FLT_MIN
#define SINGLE_XMAX FLT_MAX
#define DOUBLE_DIGITS DBL_MANT_DIG
#define DOUBLE_EPS DBL_EPSILON
#define DOUBLE_XMAX DBL_MAX
#define DOUBLE_XMIN DBL_MIN

#endif

namespace SHEsis {
// Fisher's exact test

void fexact(int* nrow, int* ncol, double* table, int* ldtabl, double* expect,
            double* percnt, double* emin, double* prt, double* pre,
            int* workspace);
}

#endif
