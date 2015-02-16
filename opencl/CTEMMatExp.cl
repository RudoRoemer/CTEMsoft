// ###################################################################
// Copyright (c) 2014, Saransh Singh/Carnegie Mellon University
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
//        this software without specific prior written permission.
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
// CTEMsoft2013:CTEMMatExp.cl
//--------------------------------------------------------------------------
//
// FUNCTION: cmplxmult
//
//> @author Saransh Singh, Carnegie Mellon University
//
//> @brief function to calculate multiplication for two compllex numbers
//
//> @date 09/23/14  SS  1.0 Original
//--------------------------------------------------------------------------

float2 cmplxmult(float2, float2);

float2 cmplxmult(float2 a, float2 b){
    float2 res;
    res.x = a.x*b.x - a.y*b.y;
    res.y = a.x*b.y + a.y*b.x;
    
    return res;
}

//--------------------------------------------------------------------------
// CTEMsoft2013:CTEMMatExp.cl
//--------------------------------------------------------------------------
//
// KERNEL: CTEMMatExp
//
//> @author Saransh Singh, Carnegie Mellon University
//
//> @brief OpenCL kernel for calculating exponential of matrix using Taylor series when multiple exponentials are calculated simultaneously
//
//> @date 09/09/14  SS  1.0 Original
//> @date 02/16/15  SS  1.1 bug fixes and simplification using cmplxmult
//--------------------------------------------------------------------------

__kernel void MatExp(__global float2* cl_expA,__global float2* cl_A,__global float2* cl_AA,__global float2* cl_AAA,const int nn,__global float2* cl_coeff,__global float2* cl_T1,__global float2* cl_T2,__global float2* cl_wavefncoeff,__global float2* cl_wavefncoeffintd,const int ns)
{
    int tx,ty,id;
    tx = get_global_id(0);
    ty = get_global_id(1);
    id = get_global_size(0)*ty + tx;
    
    float2 sum = (float2)(0.0f,0.0f);
    
    // calculating A^2 and A^3
    
    for (int i = 0; i < nn; i++){
        for (int j = 0; j < nn; j++){
            for (int k = 0; k < nn; k++){
                sum += cmplxmult(cl_A[id*nn*nn + i*nn + k],cl_A[id*nn*nn + k*nn + j]);
            }
            cl_AA[id*nn*nn + i*nn + j] = sum;
            sum = (float2)(0.0f,0.0f);
        }
    }

    
    sum = (float2)(0.0f,0.0f);
    for (int i = 0; i < nn; i++){
        for (int j = 0; j < nn; j++){
            for (int k = 0; k < nn; k++){
                sum += cmplxmult(cl_AA[id*nn*nn + i*nn + k],cl_A[id*nn*nn + k*nn + j]);
            }
            cl_AAA[id*nn*nn + i*nn + j] = sum;
            sum = (float2)(0.0f,0.0f);
        }
    }
    
    // Calculating the three factors for the exponential
    
    for (int i = 0; i < nn; i++){
        for (int j = 0; j < nn; j++){
            if ( i == j){
                cl_expA[id*nn*nn + i*nn + j] = cl_AAA[id*nn*nn + i*nn + j] - cmplxmult(cl_coeff[0],cl_AA[id*nn*nn + i*nn + j]) + cmplxmult(cl_coeff[1],cl_A[id*nn*nn + i*nn + j]) - cl_coeff[2];
            
                cl_T1[id*nn*nn + i*nn + j] = cl_AAA[id*nn*nn + i*nn + j] - cmplxmult(cl_coeff[3],cl_AA[id*nn*nn + i*nn + j]) + cmplxmult(cl_coeff[4],cl_A[id*nn*nn + i*nn + j]) - cl_coeff[5];

            }
            else {
                
                cl_expA[id*nn*nn + i*nn + j] = cl_AAA[id*nn*nn + i*nn + j] - cmplxmult(cl_coeff[0],cl_AA[id*nn*nn + i*nn + j]) + cmplxmult(cl_coeff[1],cl_A[id*nn*nn + i*nn + j]);
                
                cl_T1[id*nn*nn + i*nn + j] = cl_AAA[id*nn*nn + i*nn + j] - cmplxmult(cl_coeff[3],cl_AA[id*nn*nn + i*nn + j]) + cmplxmult(cl_coeff[4],cl_A[id*nn*nn + i*nn + j]);
                
                
            }
        }
    }

    
    sum = (float2)(0.0f,0.0f);
    for (int i = 0; i < nn; i++){
        for (int j = 0; j < nn; j++){
            for (int k = 0; k < nn; k++){
                sum += cmplxmult(cl_T1[id*nn*nn + i*nn + k],cl_expA[id*nn*nn + k*nn + j]);
            }
            cl_T2[id*nn*nn + i*nn + j] = sum;
            sum = (float2)(0.0f,0.0f);
        }
    }
    
    for (int i = 0; i < nn; i++){
        for (int j = 0; j < nn; j++){
            if ( i == j){
                
                cl_T1[id*nn*nn + i*nn + j] = cl_AAA[id*nn*nn + i*nn + j] - cmplxmult(cl_coeff[6],cl_AA[id*nn*nn + i*nn + j]) + cmplxmult(cl_coeff[7],cl_A[id*nn*nn + i*nn + j]) - cl_coeff[8];
                
            }
            else {
                
                cl_T1[id*nn*nn + i*nn + j] = cl_AAA[id*nn*nn + i*nn + j] - cmplxmult(cl_coeff[6],cl_AA[id*nn*nn + i*nn + j]) + cmplxmult(cl_coeff[7],cl_A[id*nn*nn + i*nn + j]);
                
            }
        }
    }


    sum = (float2)(0.0f,0.0f);
    for (int i = 0; i < nn; i++){
        for (int j = 0; j < nn; j++){
            for (int k = 0; k < nn; k++){
                sum += cmplxmult(cl_T2[id*nn*nn + i*nn + k],cl_T1[id*nn*nn + k*nn + j]);
            }
            sum /= 362880;
            cl_expA[id*nn*nn + i*nn + j] = sum;
            sum = (float2)(0.0f,0.0f);
        }
    }
    
    // squaring operation as matrix was scaled in the host code
    
    for (int l = 0; l < ns; l++){
        sum = (float2)(0.0f,0.0f);
        for (int i = 0; i < nn; i++){
            for (int j = 0; j < nn; j++){
                for (int k = 0; k < nn; k++){
                    
                    sum += cmplxmult(cl_expA[id*nn*nn + i*nn + k],cl_expA[id*nn*nn + k*nn + j]);

                }
                cl_T1[id*nn*nn + i*nn + j] = sum;
                sum = (float2)(0.0f,0.0f);
            }
        }

        for (int i = 0; i < nn; i++){
            for (int j = 0; j < nn; j++){
                cl_expA[id*nn*nn + i*nn + j] = cl_T1[id*nn*nn + i*nn + j];
            }
        }
    }

    // cl_expA now has the exponential of the structure matrix. We now multiply
    // by the column vector [1 0 0 ...... 0] to the the fourier coefficients of
    // wavefunction at different depths and subsequently the depth integrated
    // intensity
    

    
    for (int i = 0; i < nn; i++){
        
        if ( i == 0) {
            cl_wavefncoeff[id*nn + i] = (float2)(1.0f,0.0f);
        }
        
        else {
            cl_wavefncoeff[id*nn + i] = (float2)(0.0f,0.0f);
        }
        
    }
  
    for (int i = 0; i < 25; i++){
        
        sum = (float2)(0.0f,0.0f);
        
        for (int l = 0; l < nn; l++){
            cl_wavefncoeffintd[id*nn + l] = cl_wavefncoeff[id*nn + l];
        }
        
        for (int j = 0; j < nn; j++){
            for (int k = 0; k < nn; k++){
                
                sum += cmplxmult(cl_expA[id*nn*nn + k*nn + j],cl_wavefncoeffintd[id*nn + k]);

            }
            
            cl_wavefncoeff[id*nn + j] = sum;
            sum = (float2)(0.0f,0.0f);
            
        }
        
    }
    
    // we have the fourier coefficients of the wavefunctions.
    
    

}



