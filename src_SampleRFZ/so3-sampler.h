//
//  so3-sampler.h
//  
//
//  Created by Saransh Singh on 6/17/14.
//
//

#ifndef _so3_sampler_h
#define _so3_sampler_h

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>


// these are a bunch of constants used for the mapping; they are all in double precision

struct LambertParameters{
    double Pi;          // pi
    double iPi;         // 1/pi
    double sPi;         // sqrt(pi)
    double srt;         // sqrt(3)/2
    double pref;        // sqrt(6/pi)
// the following constants are used for the cube to quaternion hemisphere mapping
    double a;           // pi^(5/6)/6^(1/6)
    double ap;          // pi^(2/3)
    double sc;          // a/ap
    double beta;        // pi^(5/6)/6^(1/6)/2
    double R1;          // (3pi/4)^(1/3)
    double r2;          // sqrt(2)
    double r22;         // 1/sqrt(2)
    double pi12;        // pi/12
    double pi8;         // pi/8
    double prek;        // R1 2^(1/4)/beta
    double r24;         // sqrt(24)
    double rvmax2;      // square of max rodrigues vector length
    double BP[6];       // used for fundamental zone determination
};

// type definition for linked list of Rodrigues points

struct FZpointd{
    double rod[3];          // Rodrigues point
    struct FZpointd *next;  // link to next point
};

// Return type for SampleRFZ function

struct SampleRFZ_return{
    struct FZpointd *FZpointer;
    int FZcnt;
};

// Return type for LambertCubeToBall function

struct LambertCubeToBall_return{
    double rod[3];  // rodrigues vector
    int ierr;       // error flags
};

double max2(double a, double b);

double maxval(double a[3]);

double *fabsarr(double a[3]);

struct LambertParameters InitLambertParameters(void);

int GetPyramid(double xyz[3]);

struct LambertCubeToBall_return LambertCubeToBall(double lxyz[3]);

double *ho2ro(double h[3]);

double *cu2ro(double c[3]);

bool insideCyclicFZ(double rod[3], int order);

bool insideDihedralFZ(double rod[3], int order);

bool insideCubicFZ(double rod[3], char ot[3]);

struct SampleRFZ_return SampleRFZ(int nsteps, int pgnum);

#endif
