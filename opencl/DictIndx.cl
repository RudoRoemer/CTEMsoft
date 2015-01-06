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

float4 quatmult(float4,float4);
float factorial(int);
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
    float tol = 1e-5f;
    float sum;
    int k = 1;
    float term = x/4;
    sum = term;
    while (term >= tol){
        term = (x/2)*pow((x*x/4),k)/(factorial(k)*factorial(k+2));
        sum += term;
        k++;
    }
    return sum;

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
    float tol = 1e-5f;
    float sum;
    int k = 1;
    float term = x*x/24;
    sum = term;
    while (term >= tol){
        term = (x*x/4)*pow((x*x/4),k)/(factorial(k)*factorial(k+3));
        sum += term;
        k++;
    }
    return sum;
    
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
    
    float4 mu = (float4)(0.0f,0.0f,0.0f,0.0f);
    float kappa = 0.0f;
    float4 sym;
    float rim;
    
    for (int l = 0; l < 100; ++l){
        gamma = (float4)(0.0f,0.0f,0.0f,0.0f);
        for (int i = 0; i < numk; ++i){
            quatlist = (float4)(quaternion[4*id],quaternion[4*id+1],quaternion[4*id+2],quaternion[4*id+3]);
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