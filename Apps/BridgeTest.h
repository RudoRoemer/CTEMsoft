#ifndef _BRIDGETEST_H_
#define _BRIDGETEST_H_


#include "SampleRFZ_Mangling.h"

#ifdef __cplusplus
extern "C" {
#endif

	typedef void (*CallbackType)(int, double*&);

	void SAMPLERFZ_GLOBAL(sampler, SAMPLER)(int* pgnum, int* nsteps, CallbackType callback);

#ifdef __cplusplus
}
#endif



#endif /*_BRIDGETEST_H_*/
