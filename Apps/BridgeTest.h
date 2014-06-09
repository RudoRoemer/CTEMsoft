#ifndef _BRIDGETEST_H_
#define _BRIDGETEST_H_

// 	struct nmltypefile {
//  		double sig;
//  	};

#ifdef __cplusplus
extern "C" {
#endif

	typedef void (*CallbackType)(int, double*&);
	void sampler_(int* pgnum, int* nsteps, CallbackType callback);

    //void domcsimulation_(nmltypefile* nmlfile);



#ifdef __cplusplus
}
#endif



#endif /*_BRIDGETEST_H_*/
