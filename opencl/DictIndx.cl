/*
 ###################################################################
 ! Copyright (c) 2013-2014, Saransh Singh/Carnegie Mellon University
 ! All rights reserved.
 !
 ! Redistribution and use in source and binary forms, with or without modification, are
 ! permitted provided that the following conditions are met:
 !
 !     - Redistributions of source code must retain the above copyright notice, this list
 !        of conditions and the following disclaimer.
 !     - Redistributions in binary form must reproduce the above copyright notice, this
 !        list of conditions and the following disclaimer in the documentation and/or
 !        other materials provided with the distribution.
 !     - Neither the names of Marc De Graef, Carnegie Mellon University nor the names
 !        of its contributors may be used to endorse or promote products derived from
 !        this software without specific prior written permission.
 !
 ! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 ! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 ! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 ! ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 ! LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 ! DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 ! SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 ! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 ! OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
 ! USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 ! ###################################################################
 */

 #define PI 3.14159265359f
 #define e  2.71828182845f
 #define BLOCK_SIZE 16


float4 quatmult(float4,float4);
float factorial(int);
float Bessel0(float);
float Bessel1(float);
float Bessel2(float);
float Cp(float);
float CalcVMF(float4,float4,float,float4);

 /*
 //--------------------------------------------------------------------------
 //
 // FUNCTION: quatmult
 //
 //> @author Saransh Singh, Carnegie Mellon University
 //
 //> @brief function to multiply two quaternions i.e. quat1 x quat2
 //
 //> @param quat1 first quaternion
 //> @param quat2 second quaternion
 //> @date 01/06/15    SS 1.0 original
 //--------------------------------------------------------------------------
 */

 float4 quatmult(float4 quat1, float4 quat2)
{
    float4 res;
    res.x = quat1.x*quat2.x - quat1.y*quat2.y - quat1.z*quat2.z - quat1.w*quat2.w;
    res.y = quat1.y*quat2.x + quat1.x*quat2.y - quat1.w*quat2.z + quat1.z*quat2.w;
    res.z = quat1.z*quat2.x + quat1.w*quat2.y + quat1.x*quat2.z - quat1.y*quat2.w;
    res.w = quat1.w*quat2.x - quat1.z*quat2.y + quat1.y*quat2.z + quat1.x*quat2.w;
    
    return res;
}

 /*
 //--------------------------------------------------------------------------
 //
 // FUNCTION: factorial
 //
 //> @author Saransh Singh, Carnegie Mellon University
 //
 //> @brief function to compute factorial of a function
 //
 //> @param x input integer parameter
 //> @date 01/06/15    SS 1.0 original
 //--------------------------------------------------------------------------
 */

 float factorial(int x)
{
    int fact = 1;
    for (int i = 0; i <= x; ++i){
        fact *= i;
    }
    return fact;
}

 /*
 //--------------------------------------------------------------------------
 //
 // FUNCTION: Bessel0
 //
 //> @author Saransh Singh, Carnegie Mellon University
 //
 //> @brief function to compute modified Bessel function of first kind of order 0
 //
 //> @param x input parameter
 //> @date 01/06/15    SS 1.0 original
 //--------------------------------------------------------------------------
 */
float Bessel0(float x)
{
    float bessI0;
    float Y,AX,BX;
    
    float P1 = 1.0f;
    float P2 = 3.5156229f;
    float P3 = 3.0899424f;
    float P4 = 1.2067492f;
    float P5 = 0.2659732f;
    float P6 = 0.360768e-1f;
    float P7 = 0.45813e-2f;
    
    float Q1 = 0.39894228f;
    float Q2 = 0.1328592e-1f;
    float Q3 = 0.225319e-2f;
    float Q4 = 0.157565e-2f;
    float Q5 = 0.916281e-2f;
    float Q6 = -0.2057706e-1f;
    float Q7 = 0.2635537e-1f;
    float Q8 = -0.1647633e-1f;
    float Q9 = 0.392377e-2f;
    
    if (fabs(x) <= 3.75f){
        Y = (x/3.75f)*(x/3.75f);
        bessI0 = P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7)))));
    }
    else {
        AX = fabs(x);
        Y  = 3.75f/AX;
        BX = exp(AX)/sqrt(AX);
        AX = Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9)))))));
        bessI0 = AX*BX;
    }
    
    return bessI0;

}

 /*
 //--------------------------------------------------------------------------
 //
 // FUNCTION: Bessel1
 //
 //> @author Saransh Singh, Carnegie Mellon University
 //
 //> @brief function to compute modified Bessel function of first kind of order 1
 //
 //> @param x input parameter
 //> @date 01/06/15    SS 1.0 original
 //--------------------------------------------------------------------------
 */
float Bessel1(float x)
{
    float bessI1;
    float Y,AX,BX;
    
    float P1 = 0.5f;
    float P2 = 0.87890594f;
    float P3 = 0.51498869f;
    float P4 = 0.15084934f;
    float P5 = 0.2658733e-1f;
    float P6 = 0.301532e-2f;
    float P7 = 0.32411e-3f;
    
    float Q1 = 0.39894228f;
    float Q2 = -0.3988024e-1f;
    float Q3 = 0.362018e-2f;
    float Q4 = 0.163801e-2f;
    float Q5 = -0.1031555e-1f;
    float Q6 = 0.2282967e-1f;
    float Q7 = -0.2895312e-1f;
    float Q8 = 0.1787654e-1f;
    float Q9 = -0.420059e-2f;
    
    if (fabs(x) <= 3.75f){
        Y = (x/3.75f)*(x/3.75f);
        bessI1 = x*(P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))));
    }
    else {
        AX = fabs(x);
        Y  = 3.75f/AX;
        BX = exp(AX)/sqrt(AX);
        AX = Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9)))))));
        bessI1 = AX*BX;
    }
    
    return bessI1;
    
}
 /*
 //--------------------------------------------------------------------------
 //
 // FUNCTION: Bessel2
 //
 //> @author Saransh Singh, Carnegie Mellon University
 //
 //> @brief function to compute modified Bessel function of first kind of order 2
 //
 //> @param x input parameter
 //> @date 01/06/15    SS 1.0 original
 //--------------------------------------------------------------------------
 */
float Bessel2(float x)
{
    float res;
    if ( x >= 1.0e-5f){
        res = Bessel0(x) - (2/x)*Bessel1(x);
    }
    else {
        res = 0.0f;
    }
    return res;
}

 /*
 //--------------------------------------------------------------------------
 //
 // FUNCTION: Cp
 //
 //> @author Saransh Singh, Carnegie Mellon University
 //
 //> @brief function to generate Cp in the VMF distribution
 //
 //> @param kappa input value
 //> @date 01/06/15    SS 1.0 original
 //--------------------------------------------------------------------------
 */
 float Cp(float kappa)
{
    float res;
    res = kappa/((4*PI*PI)*Bessel1(kappa));
    return res;
}

 /*
 //--------------------------------------------------------------------------
 //
 // FUNCTION: CalcVMF
 //
 //> @author Saransh Singh, Carnegie Mellon University
 //
 //> @brief function to compute VMF density
 //
 //> @param quat unit quaternion
 //> @param mu mean value
 //> @param kappa concentartion parameter
 //> @date 01/06/15    SS 1.0 original
 //--------------------------------------------------------------------------
 */
 float CalcVMF(float4 quat, float4 mu, float kappa, float4 sym)
{
    float VMF;
    float4 meansym;
    float dp;
    meansym = quatmult(sym,mu);
    dp = meansym.x*quat.x + meansym.y*quat.y + meansym.z*quat.z + meansym.w*quat.w;
    VMF = Cp(kappa)*exp(kappa*dp);
    return VMF;
}

 /*
 !--------------------------------------------------------------------------
 !
 ! PROGRAM:ParamEstm
 !
 !> @author Saransh Singh, Carnegie Mellon University
 !
 !> @brief perform parameter estimation for the dictionary indexing
 !
 !> @param quaternion list of quaternions for which parameter estimation is done
 !> @param mean output variable containing the mean values i.e. mu
 !> @param kappa output variable containing the concentration parameter
 !> @param symmetry the list of symmetry operators
 !> @param numk number of quaternion from which mean is calculated
 !> @param numsym number of symmetry elements in the symmetry matrix
 !
 !> @date 01/06/15  SS 1.0 original
 !--------------------------------------------------------------------------
 */

__kernel void ParamEstm(__global float* quaternion, __global float* mean, __global float* ConcParam, __global float* symmetry, const int numk, const int numsym)
{
    int tx = get_global_id(0);
    int ty = get_global_id(1);
    
    int id = get_global_size(0)*ty + tx;
    
    float4 quatlist;
    
    float4 gamma = (float4)(0.0f,0.0f,0.0f,0.0f);
    float modgamma;
    float prefact;
    
    float4 mu = (float4)(1.0f,1.0f,1.0f,1.0f);
    float kappa = 1.0f;
    float4 sym;
    float rim;
    
    for (int l = 0; l < 30; ++l){
        gamma = (float4)(0.0f,0.0f,0.0f,0.0f);
        for (int i = 0; i < numk; ++i){
            quatlist = (float4)(quaternion[4*id*numk+i],quaternion[4*id*numk+i+1],quaternion[4*id*numk+i+2],quaternion[4*id*numk+i+3]);
            prefact = 0.0f;
            for (int j = 0; j < numsym; ++j){
                sym = (float4)(symmetry[4*j],symmetry[4*j+1],symmetry[4*j+2],symmetry[4*j+3]);
                prefact += CalcVMF(quatlist,mu,kappa,sym);
            }
            
            for (int k = 0; k < numsym; ++k){
                sym = (float4)(symmetry[4*k],symmetry[4*k+1],symmetry[4*k+2],symmetry[4*k+3]);
                rim = CalcVMF(quatlist,mu,kappa,sym)/prefact;
                gamma += rim*quatmult(sym,quatlist);
                
            }
        }
        modgamma = gamma.x*gamma.x + gamma.y*gamma.y + gamma.z*gamma.z + gamma.w*gamma.w;
        modgamma = sqrt(modgamma);
        mu = gamma/modgamma;
        kappa = Bessel1(modgamma/numk)/Bessel2(modgamma/numk);
    }
    
    mean[4*id] = mu.x;
    mean[4*id + 1] = mu.y;
    mean[4*id + 2] = mu.z;
    mean[4*id + 3] = mu.w;
    ConcParam[id] = kappa;
    
}
/*
!--------------------------------------------------------------------------
!
! PROGRAM:InnerProd
!
!> @author Saransh Singh, Carnegie Mellon University
!
!> @brief perform inner product calculations for normalized images as a matrix multiplication
!
!> @param exp experimental pattern chunk
!> @param dict dictionary pattern chunk
!> @param L size of one image in pixels
!> @param Ne number of experimental patterns in one chunk
!> @param Nd number of dictionary patterns in one chunk
!> @param result result of the dot product
!
!> @date 12/09/14  SS 1.0 original
!--------------------------------------------------------------------------
*/

//#define AS(i,j) As[i * 64 + j]
//#define BS(i,j) Bs[i * 64 + j]

// Notice that the Block Size is pre-set to 16. This needs to maybe become more flexible in the future.

__kernel void InnerProd(__global float* exp, __global float* dict, int Wexp, int Wdict, __global float* result)
{
    // Block index
    int bx = get_group_id(0);
    int by = get_group_id(1);
    
    // Thread index inside the block
    int tx = get_local_id(0);
    int ty = get_local_id(1);
    
    int aBegin = Wexp * BLOCK_SIZE * by;
    int aEnd = aBegin + Wexp - 1;
    int aStep = BLOCK_SIZE;
    
    int bBegin = BLOCK_SIZE * bx;
    //int bEnd = bBegin + Wdict * (get_num_groups(1));
    int bstep = BLOCK_SIZE * Wdict;
    float Csub = 0.0f;
    
    __local float As[BLOCK_SIZE][BLOCK_SIZE];
    __local float Bs[BLOCK_SIZE][BLOCK_SIZE];
    
    for (int a = aBegin, b = bBegin; /*b <= bEnd*/a <= aEnd; a += aStep, b += bstep){
        
        As[ty][tx] = exp[a + Wexp * ty + tx];
        
        Bs[ty][tx] = dict[b + Wdict * ty + tx];
        
        barrier(CLK_LOCAL_MEM_FENCE);
        
        #pragma unroll
        for (int k = 0; k < BLOCK_SIZE; ++k){
            Csub += As[ty][k] * Bs[k][tx];
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    //int c = Wdict * BLOCK_SIZE * by + BLOCK_SIZE * bx;
    result[get_global_id(1) * get_global_size(0) + get_global_id(0)] = Csub;
    //result[get_global_id(1) * get_global_size(0) + get_global_id(0)] = float(tx);
    //result[c + Wdict*ty + tx] = Csub;
    
}