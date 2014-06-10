#ifndef _BRIDGETEST_H_
#define _BRIDGETEST_H_


#include "SampleRFZ_Mangling.h"

struct thecstruct {
int m, n;
float r;
double d;
};



#ifdef __cplusplus
extern "C" {
#endif

	typedef void (*CallbackType)(int, double*&);

	void SAMPLERFZ_GLOBAL(sampler, SAMPLER)(int* pgnum, int* nsteps, thecstruct* mycstruct, CallbackType callback);

#ifdef __cplusplus
}
#endif



#endif /*_BRIDGETEST_H_*/
