//
//  so3-sampler.c
//
//
//  Created by Saransh Singh on 5/14/14.
//
//

// ###################################################################
// Copyright (c) 2014, Saransh Singh, Marc De Graef/Carnegie Mellon University
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification, are
// permitted provided that the following conditions are met:
//
//     - Redistributions of source code must retain the above copyright notice, this list
//        of conditions and the following disclaimer.
//     - Redistributions in binary form must reproduce the above copyright notice, this
//        list of conditions and the following disclaimer in the documentation and/or
//        other materials provided with the distribution.
//     - Neither the names of Marc De Graef, Carnegie Mellon University nor the names
//        of its contributors may be used to endorse or promote products derived from
//       this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
// USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// ###################################################################

//--------------------------------------------------------------------------
// CTEMsoft2013:so3 module
//--------------------------------------------------------------------------
//
// HEADER FILE: so3
//
//> @author Saransh Singh, Carnegie Mellon University
//
//> @brief
//> routines needed for uniform sampling of so3 via the cubochoric parameterization
//
//> @details This module is used as follows:
//>
//> SampleRFZ(pgnum, nsteps)
//>
//> where pgnum is the point group number (International Tables) and nsteps is the
//> required number of sampling steps along the cubic semi-edge.  The SampleRFZ
//> routine creates a linked list named FZlist which has FZcnt entries in it.  The
//> individual entries are the Rodrigues vector components and a link to the next
//> entry.  So, going through the list afterwards can be done as follows:
//>
//> FZtmp = FZlist                         ! point to the top of the list
//> do i = 1, FZcnt                        ! loop over all entries
//>   do something with FZtmp->rod         ! anything, really ...
//>   FZtmp = FZtmp->next                  ! point to the next entry
//> end do
//>
//
//> @date 5/14/14 SS 1.0 original, extracted from 5/12/14 version of CTEMsoft2013 main library
//--------------------------------------------------------------------------


#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <stdbool.h>
#include "so3-sampler.h"

// define function to find greater of two values

double max2(double a, double b){
    if (a >= b) {
        return a;
    }
    return b;
}

// define the equivalent of the maxval function in fortran (for a 3 element array)

double maxval(double a[3]){
    double res;
    res = a[0];
    for (int i = 0; i < 3; i++) {
        if (a[i] >= res) {
            res = a[i];
        }
    }
    return res;
}

// define function for fabsolute value of an array (for a 3 element array)

double *fabsarr(double a[3]){
    static double b[3];
    for(int i = 0; i < 3; i++) {
        if (a[i] < 0) {
            b[i] = -a[i];
        }
        else{
            b[i] = a[i];
        }
        
    }
    return b;
}


// The following two arrays are used to determine the FZtype (FZtarray) and primary rotation axis order (FZoarray)
// for each of the 32 crystallographic point group symmetries (in the order of the International Tables)
//
// 1 (C1), -1 (Ci), [triclinic]
// 2 (C2), m (Cs), 2/m (C2h), [monoclinic]
// 222 (D2), mm2 (C2v), mmm (D2h), [orthorhombic]
// 4 (C4), -4 (S4), 4/m (C4h), 422 (D4), 4mm (C4v), -42m (D2d), 4/mmm (D4h), [tetragonal]
// 3 (C3), -3 (C3i), 32 (D3), 3m (C3v), -3m (D3d), [trigonal]
// 6 (C6), -6 (C3h), 6/m (C6h), 622 (D6), 6mm (C6v), -6m2 (D3h), 6/mmm (D6h), [hexagonal]
// 23 (T), m3 (Th), 432 (O), -43m (Td), m-3m (Oh) [cubic]
//
// FZtype
// 0        no symmetry at all
// 1        cyclic symmetry
// 2        dihedral symmetry
// 3        tetrahedral symmetry
// 4        octahedral symmetry
//


struct LambertParameters InitLambertParameters(void){
    double dpi;
    struct LambertParameters LPs;
    dpi = 4.0*atan(1.0);
    
    LPs.Pi = dpi;
    LPs.iPi = 1.0/dpi;
    LPs.sPi = sqrt(dpi);
    LPs.srt = sqrt(3.0)/2.0;
    LPs.pref = sqrt(6.0/dpi);
    
    LPs.a = pow(dpi,(5.0/6.0))/pow(6.0,(1.0/6.0));
    LPs.ap = pow(dpi,(2.0/3.0));
    LPs.sc = LPs.a/LPs.ap;
    LPs.beta = 0.50*LPs.a;
    LPs.R1 = pow((3.0*dpi/4.0),(1.0/3.0));
    LPs.r2 = sqrt(2.0);
    LPs.r22 = sqrt(2.0)/2.0;
    LPs.pi12 = dpi/12.0;
    LPs.pi8 = tan(dpi/8.0);
    LPs.prek =  (LPs.R1*pow(2.0,0.25))/LPs.beta;
    LPs.r24 = sqrt(24.0);
    
    // we need to make sure that we do not consider Rodrigues vectors of infinite length,
    // which would correspond to 180 degree rotations; we'll put the max length so that
    // we can deal with 179.999 degree rotations, which should be good enough for most
    // practical cases
    
    LPs.rvmax2 = pow(tan(179.999*dpi/(2.0*180.0)),2);      // square of max rodrigues length
    
    LPs.BP[0] = 0.0;
    LPs.BP[1] = tan(dpi/4.0);
    LPs.BP[2] = tan(dpi/6.0);
    LPs.BP[3] = LPs.pi8;
    LPs.BP[4] = 0.0;
    LPs.BP[5] = tan(LPs.pi12);
    return LPs;
    
}

//-------------------------------------------------------------
//-------------------------------------------------------------
//
// Check the paper by Rosca, Morawiec, and De Graef
// for details on this cube-to-ball mapping, in particular
// the labeling of the pyramids that make up the cube.
//
// Citation:  "A new method of constructing a grid in
// the space of 3D rotations and its applications to
// texture analysis," D. Rosca, A. Morawiec, and M. De Graef,
// submitted to Modeling and Simulations in Materials Science
// and Engineering (April 2014).
//
//-------------------------------------------------------------
//-------------------------------------------------------------

//--------------------------------------------------------------------------
//
// FUNCTION: GetPyramid
//
//> @author Saransh Singh, Carnegie Mellon University
//
//> @brief determine to which pyramid a point in a cubic grid belongs
//
//> @details check the first figure in the cube-to-ball paper for the pyramid labels
//
//> @param xyz 3D coordinates to be considered (double precision)
//> @date 05/14/14    SS 1.0 original
//--------------------------------------------------------------------------


int GetPyramid(double xyz[3]){
    int res;
    
    if ((fabs(xyz[0]) <= xyz[2]) && (fabs(xyz[1]) <= xyz[2])) {
        res = 1;        // pyramid 1
        return res;
    }
    
    if ((fabs(xyz[0]) <= -xyz[2]) && (fabs(xyz[1]) <= -xyz[2])) {
        res = 2;        // pyramid 2
        return res;
    }
    
    if ((fabs(xyz[2]) <= xyz[0]) && (fabs(xyz[1]) <= xyz[0])) {
        res = 3;        // pyramid 3
        return res;
    }
    
    if ((fabs(xyz[2]) <= -xyz[0]) && (fabs(xyz[1]) <= -xyz[0])) {
        res = 4;        // pyramid 4
        return res;
    }
    
    if ((fabs(xyz[0]) <= xyz[1]) && (fabs(xyz[2]) <= xyz[1])) {
        res = 5;        // pyramid 5
        return res;
    }
    
    res = 6;            // pyramid 6
    
    return res;
    
}


//--------------------------------------------------------------------------
//
// FUNCTION: LambertCubeToBall
//
//> @author Saransh Singh, Carnegie Mellon University
//
//> @brief map from 3D cubic grid to 3D ball
//
//> @param lxyz 3D coordinates to be considered (double precision)
//> @param ierr error flag 0 = OK, 1 = outside of unit cube
//
//> @date 5/14/14    SS 1.0 original
//--------------------------------------------------------------------------


struct LambertCubeToBall_return LambertCubeToBall(double lxyz[3]){
    static double res[3];
    struct LambertCubeToBall_return ret;
    struct LambertParameters LPs;
    LPs = InitLambertParameters();
    double XYZ[3], sXYZ[3], T1, T2, c, s, q, LamXYZ[3], edge;
    int p, ierr;
    ierr = 0;
    edge = 0.5 * LPs.ap;
    if (maxval(fabsarr(lxyz)) > edge) {
        res[0] = 0.0;
        res[1] = 0.0;
        res[2] = 0.0;
        ierr = 1;
        ret.rod[0] = res[0];
        ret.rod[1] = res[1];
        ret.rod[2] = res[2];
        ret.ierr = ierr;
        return ret;
    }
    // determine which pyramid pair the point lies in and copy coordinates in correct order (see paper)
    p = GetPyramid(lxyz);
    switch(p){
        case 1:
            sXYZ[0] = lxyz[0];
            sXYZ[1] = lxyz[1];
            sXYZ[2] = lxyz[2];
            break;
        case 2:
            sXYZ[0] = lxyz[0];
            sXYZ[1] = lxyz[1];
            sXYZ[2] = lxyz[2];
            break;
        case 3:
            sXYZ[0] = lxyz[1];
            sXYZ[1] = lxyz[2];
            sXYZ[2] = lxyz[0];
            break;
        case 4:
            sXYZ[0] = lxyz[1];
            sXYZ[1] = lxyz[2];
            sXYZ[2] = lxyz[0];
            break;
        case 5:
            sXYZ[0] = lxyz[2];
            sXYZ[1] = lxyz[0];
            sXYZ[2] = lxyz[1];
            break;
        case 6:
            sXYZ[0] = lxyz[2];
            sXYZ[1] = lxyz[0];
            sXYZ[2] = lxyz[1];
            break;
    }
    // scale by grid parameter ratio sc
    XYZ[0] = LPs.sc * sXYZ[0];
    XYZ[1] = LPs.sc * sXYZ[1];
    XYZ[2] = LPs.sc * sXYZ[2];
    
    
    // transform to the sphere grid via the curved square, and intercept the zero point
    if (maxval(fabsarr(XYZ)) == 0.0){
        LamXYZ[0] = 0.0;
        LamXYZ[1] = 0.0;
        LamXYZ[2] = 0.0;
    }
    else {
        // intercept all the points along the z-axis
        if (XYZ[0]== 0.0 && XYZ[1] == 0.0){
            LamXYZ[0] = 0.0;
            LamXYZ[1] = 0.0;
            LamXYZ[2] = LPs.pref * XYZ[2];
        }
        else { // this is a general grid point
            if (fabs(XYZ[1]) <= fabs(XYZ[0])){
                c = cos(LPs.pi12 * XYZ[1]/XYZ[0]);
                s = sin(LPs.pi12 * XYZ[1]/XYZ[0]);
                q = LPs.prek * XYZ[0]/sqrt(LPs.r2 - c);
                T1 = (LPs.r2*c - 1.0) * q;
                T2 = LPs.r2 * s * q;
            }
            else {
                c = cos(LPs.pi12 * XYZ[0]/XYZ[1]);
                s = sin(LPs.pi12 * XYZ[0]/XYZ[1]);
                q = LPs.prek * XYZ[1] / sqrt(LPs.r2 - c);
                T1 = LPs.r2 * s * q;
                T2 = (LPs.r2*c - 1.0) * q;
            }
            //transform to sphere grid (inverse Lambert)
            //[note that there is no need to worry about dividing by zero, since XYZ[2] can not become zero]
            c = T1*T1 + T2*T2;
            s = LPs.Pi * c/(24.0 * XYZ[2]*XYZ[2]);
            c = (LPs.sPi * c) / (LPs.r24 * XYZ[2]);
            q = sqrt( 1.0 - s );
            LamXYZ[0] = T1 * q;
            LamXYZ[1] = T2 * q;
            LamXYZ[2] = LPs.pref * XYZ[2] - c;
            
        }
    }
    
    // reverse the coordinates back to the regular order according to the original pyramid number
    
    switch(p){
        case 1:
            res[0] = LamXYZ[0];
            res[1] = LamXYZ[1];
            res[2] = LamXYZ[2];
            break;
        case 2:
            res[0] = LamXYZ[0];
            res[1] = LamXYZ[1];
            res[2] = LamXYZ[2];
            break;
        case 3:
            res[0] = LamXYZ[2];
            res[1] = LamXYZ[0];
            res[2] = LamXYZ[1];
            break;
        case 4:
            res[0] = LamXYZ[2];
            res[1] = LamXYZ[0];
            res[2] = LamXYZ[1];
            break;
        case 5:
            res[0] = LamXYZ[1];
            res[1] = LamXYZ[2];
            res[2] = LamXYZ[0];
            break;
        case 6:
            res[0] = LamXYZ[1];
            res[1] = LamXYZ[2];
            res[2] = LamXYZ[0];
            break;
    }
    ret.rod[0] = res[0];
    ret.rod[1] = res[1];
    ret.rod[2] = res[2];
    ret.ierr = ierr;
    return ret;
}

//--------------------------------------------------------------------------
//
// FUNCTION: ho2ro
//
//> @author Saransh Singh, Carnegie Mellon University
//
//> @brief convert homochoric to Rodrigues
//
//> @param h homochoric coordinates (double precision)
//
//> @date 5/16/14   SS 1.0 original
//--------------------------------------------------------------------------

double *ho2ro(double h[3]){
    static double res[3];
    struct LambertParameters LPs;
    LPs = InitLambertParameters();
    int i;
    double hn[3];
    double hmag, s, hm;
    // fit parameters determined with Mathematica
    const double c[] = { -0.50000961491703210, -0.024866061488717310,
        -0.0045493817793628190, 0.00051186683663875260,
        -0.00165008273335755480, 0.00075933522033887180,
        -0.00020404225025668760 };
    // normalize h and store the magnitude
    hmag = h[0]*h[0] + h[1]*h[1] + h[2]*h[2];
    if (hmag == 0.0){
        res[0] = 0.0;
        res[1] = 0.0;
        res[2] = 0.0;
        return res;
    }
    else {
        hm = hmag;
        hn[0] = h[0] / sqrt(hmag);
        hn[1] = h[1] / sqrt(hmag);
        hn[2] = h[2] / sqrt(hmag);
        // convert the magnitude to the rotation angle
        s = c[0] * hmag;
        for ( i = 1; i<7 ; i++){
            hm *= hmag;
            s += c[i] * hm;
        }
        s = tan(acos(1.0 + s));
        res[0] = hn[0]*s;
        res[1] = hn[1]*s;
        res[2] = hn[2]*s;
        return res;
    }
    
}

//--------------------------------------------------------------------------
//
// FUNCTION: cu2ro
//
//> @author Saransh Singh, Carnegie Mellon University
//
//> @brief convert cubochoric to Rodrigues
//
//> @param c cubochoric coordinates  (double precision)
//
//> @note calling program MUST have initialized the Lambert parameters first!!!
//
//> @date 5/16/14   SS 1.0 original
//--------------------------------------------------------------------------

double *cu2ro(double c[3]){
    
    double *res;
    double cc[3];
    struct LambertCubeToBall_return tmp;
    struct LambertParameters LPs;
    LPs = InitLambertParameters();
    
    for (int i = 0; i < 3; i++){
        cc[i] = c[i];
    }
    // calling program must have initialized the Lambert parameters!!!!
    tmp = LambertCubeToBall(cc);
    
    // if ierr=1, then the input point does not lie inside the sampling cube.
    // the calling program should make sure that this never happens, but if
    // it does, we need to alert the user and abort the program right here...
    
    if(tmp.ierr == 1){
        printf("Fatal Error: the sampling point coordinates are outside sampling cube...\n");
        printf("%f %f %f\n",c[0],c[1],c[2]);
        printf("Sampling cube has semi edge length %f\n",0.5*LPs.ap);
        exit(0);
        
    }
    res = ho2ro(tmp.rod);
    return res;
    
}


//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
// next, we define a number of logical routines, that decide whether or not
// a point in Rodrigues representation lies inside the fundamental zone (FZ)
// for a given crystal symmetry. This follows the Morawiec@Field paper:
//
// A. Morawiec & D. P. Field (1996) Rodrigues parameterization for orientation
// and misorientation distributions, Philosophical Magazine A, 73:4, 1113-1130,
// DOI: 10.1080/01418619608243708
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------


//--------------------------------------------------------------------------
//
// FUNCTION: insideCyclicFZ
//
//> @author Saransh Singh, Carnegie Mellon University
//
//> @brief does rodrigues point lie inside cyclic FZ (for 2, 3, 4, and 6-fold)?
//
//> @param rod Rodrigues coordinates  (double precision)
//> @param order depending on main symmetry axis
//
//> @note calling program MUST have initialized the Lambert parameters first!!!
//
//> @date 5/27/14   SS 1.0 original
//--------------------------------------------------------------------------

bool insideCyclicFZ(double rod[3], int order){
    
    bool res;
    struct LambertParameters LPs;
    LPs = InitLambertParameters();
    res = (fabs(rod[2]) <= LPs.BP[order-1]);
    return res;
}

//--------------------------------------------------------------------------
//
// FUNCTION: insideDihedralFZ
//
//> @author Saransh Singh, Carnegie Mellon University
//
//> @brief does Rodrigues point lie inside dihedral FZ (for 2, 3, 4, and 6-fold)?
//
//> @param rod Rodrigues coordinates (double precision)
//> @param order depending on main symmetry axis
//
//> @note calling program MUST have initialized the Lambert parameters first!!!
//
//> @todo we ignore here the fact that, among others, the 3m point group can be oriented in two ways;
//> @todo this should be fixable in the future with an additional optional argument
//
//> @date 5/27/14   SS 1.0 original
//--------------------------------------------------------------------------

bool insideDihedralFZ(double rod[3], int order){
    
    bool res, c1, c2;
    const double r1 = 1.0;
    struct LambertParameters LPs;
    LPs = InitLambertParameters();
    // first, check the z-component vs. tan(pi/2n)  (same as insideCyclicFZ)
    c1 = (fabs(rod[2]) <= LPs.BP[order-1]);
    res = 0;
    
    // check the square boundary planes if c1=.1.
    if(c1 == 1){
        switch(order){
            case 2:
                c2 = (fabs(rod[0]) <= r1) && (fabs(rod[1]) <= r1);
                break;
            case 3:
                c2 = fabs(LPs.srt*rod[0] + 0.5*rod[1]) <= r1;
                c2 = c2 && fabs(LPs.srt*rod[0] - 0.5*rod[1]) <=r1;
                c2 = (c2 && (fabs(rod[1]) <=r1));
                break;
            case 4:
                c2 = (fabs(rod[0]) <= r1) && (fabs(rod[1]) <= r1);
                c2 = (c2 && ((LPs.r22*fabs(rod[0] + rod[1])) <= r1) && ((LPs.r22*fabs(rod[0] - rod[1])) <= r1));
                break;
            case 6:
                c2 = (fabs(0.5*rod[0] + LPs.srt*rod[1]) <= r1);
                c2 = c2 && ((fabs(LPs.srt*rod[0] + 0.5*rod[1])) <= r1);
                c2 = c2 && ((fabs(LPs.srt*rod[0] - 0.5*rod[1])) <= r1);
                c2 = c2 && ((fabs(0.5*rod[0] - LPs.srt*rod[1])) <= r1);
                c2 = c2 && ((fabs(rod[1])) <= r1);
                c2 = c2 && ((fabs(rod[0])) <= r1);
                break;
                
        }
        res = c2;
    }
    
    return res;
}

//--------------------------------------------------------------------------
//
// FUNCTION: insideCubicFZ
//
//> @author Saransh Singh, Carnegie Mellon University
//
//> @brief does rodrigues point lie inside cubic FZ (octahedral or tetrahedral)?
//
//> @param rod Rodrigues coordinates  (double precision)
//> @param ot 'oct' or 'tet', depending on symmetry
//
//> @note calling program MUST have initialized the Lambert parameters first!!!
//
//> @date 5/28/14   SS 1.0 original
//--------------------------------------------------------------------------

bool insideCubicFZ(double rod[3], char ot[3]){
    bool res, c1, c2;
    double rx, ry, rz;
    const double r1 = 1.0;
    struct LambertParameters LPs;
    LPs = InitLambertParameters();
    
    rx = rod[0];
    ry = rod[1];
    rz = rod[2];
    
    res = 0;
    
    // primary cube planes (only needed for octahedral case)
    if (strcmp(ot,"oct") == 0){
        c1 = (maxval(fabsarr(rod)) <= LPs.pi8);
    }
    else{
        c1 = 1;
    }
    // octahedral truncation planes, both for tetrahedral and octahedral point groups
    
    c2 = (fabs(rx + ry + rz) <= r1);
    c2 = c2 && (fabs(-rx + ry + rz) <= r1);
    c2 = c2 && (fabs(rx - ry + rz) <= r1);
    c2 = c2 && (fabs(rx + ry - rz) <= r1);
    
    if (c1 && c2){
        res = 1;
    }
    return res;
}

//--------------------------------------------------------------------------
//
// FUNCTION: SampleRFZ
//
//> @author Saransh Singh, Carnegie Mellon University
//
//> @brief Generate a uniform sampling of a Rodriguess FZ
//
//> @note This routine fills in a linked list FZlist of Rodrigues points that
//> are inside a specific fundamental zone determined by the sample point group;
//> this list can then be further dealt with in the calling program.
//
//> @param nsteps number of steps along semi-edge in cubochoric grid
//> @param pgnum point group number to determine the appropriate Rodrigues fundamental zone
//
//> @note calling program MUST have initialized the Lambert parameters first!!!
//
//> @date 5/28/14   MDG 1.0 original
//--------------------------------------------------------------------------

struct SampleRFZ_return SampleRFZ(int nsteps, int pgnum){
    
    const int FZtarray[] = { 0,0,1,1,1,2,2,2,1,1,1,2,2,2,2,1,1,2,2,2,1,1,1,2,2,2,2,3,3,4,3,4 };
    const int FZoarray[] = { 0,0,2,2,2,2,2,2,4,4,4,4,4,4,4,3,3,3,3,3,6,6,6,6,6,6,6,0,0,0,0,0 };
    double x, y, z, s, delta, rval;
    double *rod;
    int FZtype, FZorder;
    bool insideFZ;
    struct FZpointd *FZlist;    // pointer to root
    struct FZpointd *FZtmp;     // temporary pointer
    int FZcnt;                  // counts number of entries
    struct SampleRFZ_return res;

    // initialize all special constants
    struct LambertParameters LPs;
    LPs = InitLambertParameters();
    
    // cube semi-edge length
    s = 0.5*LPs.ap;
    
    // step size for sampling of grid; total number of samples = (2*nsteps+1)^3
    delta = s/(double)nsteps;
    
    
    // declare linked list pointers and counters
    
    struct FZpointd *FZtmp2;
    
    // set the counter to zero
    FZcnt = 0;
    
    // make sure the linked lists are empty
    // Commenting this part out
    //Could lead to errors as memory block of FZtmp2 was not assigned using malloc
    
    /*if (FZlist != NULL) {
        FZtmp = FZlist->next;
        FZtmp2 = FZlist;
        while (1) {
            free(FZtmp2);
            FZtmp2 = NULL;
            if (FZtmp == NULL) {
                break;
            }
            FZtmp2 = FZtmp;
            FZtmp = FZtmp->next;
        }
        FZlist = NULL;
    }*/
    
    
    // determine which function we should call for this point group symmetry
    FZtype = FZtarray[pgnum-1];
    FZorder = FZoarray[pgnum-1];
    
    // loop over the cube of volume pi^2; note that we do not want to include
    // the opposite edges/facets of the cube, to avoid double counting rotations
    // with a rotation angle of 180 degrees.
    
    x = -s;
    while (x < s) {
        y = -s;
        while (y < s) {
            z = -s;
            while (z < s) {
                // convert to Rodrigues representation
                double c[] = { x, y, z };
                rod = cu2ro(c);
                // is this point inside the selected Rodrigues FZ ?
                // limit the divergent Rodrigues space to 179.999 degrees, which
                // corresponds to a length of the rod vector of 113924.
                // (actually, LPs%rvmax2 is the square of this value, to avoid having to take square roots)
                // this only applies to FZtypes 0 and 1; the other FZs are always finite.
                switch (FZtype) {
                    case 0:
                        rval = rod[0]*rod[0] + rod[1]*rod[1] + rod[2]*rod[2];
                        if (rval < LPs.rvmax2) {
                            insideFZ = 1;
                        }
                        break;
                        
                    case 1:
                        rval = rod[0]*rod[0] + rod[1]*rod[1] + rod[2]*rod[2];
                        if (rval < LPs.rvmax2) {
                            insideFZ = insideCyclicFZ(rod, FZorder);
                        }
                        break;
                        
                    case 2:
                        insideFZ = insideDihedralFZ(rod, FZorder);
                        break;
                        
                    case 3:
                        insideFZ = insideCubicFZ(rod, "tet");
                        break;
                        
                    case 4:
                        insideFZ = insideCubicFZ(rod, "oct");
                        break;
                }
                
                // If insideFZ=TRUE, then add this point to the linked list FZlist and keep
                // track of how many points there are on this list
                
                if (insideFZ){
                    if (FZlist == NULL) {
                        FZlist = (struct FZpointd *) malloc( sizeof(struct FZpointd) );
                        FZtmp = FZlist;
                        FZtmp->next = NULL;
                        FZtmp->rod[0] = rod[0];
                        FZtmp->rod[1] = rod[1];
                        FZtmp->rod[2] = rod[2];
                    }
                    
                    else {
                        FZtmp->next = (struct FZpointd *) malloc( sizeof(struct FZpointd) );
                        FZtmp = FZtmp->next;
                        FZtmp->next = NULL;
                        FZtmp->rod[0] = rod[0];
                        FZtmp->rod[1] = rod[1];
                        FZtmp->rod[2] = rod[2];
                    }
                    
                    FZcnt++;
                    
                }
                
                z += delta;
            }
            
            y += delta;
            
        }
        
        x += delta;
        
    }
    res.FZpointer = FZlist;
    res.FZcnt = FZcnt;
    return res;
}
