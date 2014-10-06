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
float2 cmplxmult(float2,float2);

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
//> @brief OpenCL kernel for calculating exponential of matrix using Taylor series when matrix is spread across the whole device
//
//> @date 09/09/14  SS  1.0 Original
//--------------------------------------------------------------------------
__kernel void MatExp(__global float2* cl_expA,__global float2* cl_A,__global float2* cl_AA,__global float2* cl_AAA,const int nn,__global float2* cl_coeff,__global float2* cl_T1,__global float2* cl_T2,__global float2* cl_T3)
{
    int tx,ty;
    tx = get_global_id(0);
    ty = get_global_id(1);
    
    float2 value;
    value = (float2)(0.0f,0.0f);
    
    for(int k = 0; k < nn; ++k){
        value.x += cl_A[ty * nn + k].x*cl_A[k * nn + tx].x - cl_A[ty * nn + k].y*cl_A[k * nn + tx].y;
        value.y = cl_A[ty * nn + k].x*cl_A[k * nn + tx].y + cl_A[ty * nn + k].y*cl_A[k * nn + tx].x;
    }
    cl_AA[ty*nn+tx] = value;
    
    //barrier(CLK_GLOBAL_MEM_FENCE);
    
    value = (float2)(0.0f,0.0f);

    for(int k = 0; k < nn; ++k){
        value.x += cl_AA[ty * nn + k].x*cl_A[k * nn + tx].x - cl_AA[ty * nn + k].y*cl_A[k * nn + tx].y;
        value.y = cl_AA[ty * nn + k].x*cl_A[k * nn + tx].y + cl_AA[ty * nn + k].y*cl_A[k * nn + tx].x;
    }
    
    cl_AAA[ty*nn+tx] = value;
    
    //barrier(CLK_GLOBAL_MEM_FENCE);

    cl_T1[ty*nn + tx] = cl_AAA[ty*nn + tx] -cl_coeff[1]*cl_AA[ty*nn + tx] + cl_coeff[2]*cl_A[ty*nn + tx] - cl_coeff[3];
    cl_T2[ty*nn + tx] = cl_AAA[ty*nn + tx] -cl_coeff[4]*cl_AA[ty*nn + tx] + cl_coeff[5]*cl_A[ty*nn + tx] - cl_coeff[6];
    cl_T3[ty*nn + tx] = cl_AAA[ty*nn + tx] -cl_coeff[7]*cl_AA[ty*nn + tx] + cl_coeff[8]*cl_A[ty*nn + tx] - cl_coeff[9];
    
    //barrier(CLK_GLOBAL_MEM_FENCE);

    value = (float2)(0.0f,0.0f);
    
    for(int k = 0; k < nn; ++k){
        value.x += cl_T2[ty*nn + k].x*cl_T3[tx*nn + k].x - cl_T2[ty*nn + k].y*cl_T3[tx*nn + k].y;
        value.y += cl_T2[ty*nn + k].x*cl_T3[tx*nn + k].y - cl_T2[ty*nn + k].y*cl_T3[tx*nn + k].x;
    }
    
    cl_AA[ty*nn+tx] = value;
    
    //barrier(CLK_GLOBAL_MEM_FENCE);
    
    value = (float2)(0.0f,0.0f);

    for(int k = 0; k < nn; ++k){
        value.x += cl_T1[ty*nn + k].x*cl_AA[tx*nn + k].x - cl_T1[ty*nn + k].y*cl_AA[tx*nn + k].y;
        value.y += cl_T1[ty*nn + k].x*cl_AA[tx*nn + k].y - cl_T1[ty*nn + k].y*cl_AA[tx*nn + k].x;
    }

    cl_expA[ty*nn+tx] = value;


}

//--------------------------------------------------------------------------
// CTEMsoft2013:CTEMMatExp.cl
//--------------------------------------------------------------------------
//
// KERNEL: CTEMMatExpImg
//
//> @author Saransh Singh, Carnegie Mellon University
//
//> @brief OpenCL kernel for calculating exponential of matrix using Taylor series when multiple exponentials are calculated simultaneously
//
//> @date 09/09/14  SS  1.0 Original
//--------------------------------------------------------------------------

__kernel void MatExpImg(__global float2* cl_expA,__global float2* cl_A,__global float2* cl_AA,__global float2* cl_AAA,const int nn,__global float2* cl_coeff,__global float2* cl_T1,__global float2* cl_T2,__global float2* cl_T3,__global float2* cl_T1T2,__global float2* cl_T1T2T3,__global float2* cl_sqr,const int ns,const int threadcnt)
{
    int tx,ty,id;
    tx = get_global_id(0);
    ty = get_global_id(1);
    id = threadcnt*ty + tx;
    
    //float2 Aloc[nn][nn];
    //for (int i = 0; i < nn; i++){
    //    for (int j = 0; i < nn; j++){
    //        Aloc[i][j] = cl_A[id*nn*nn + i*nn + j];
    //    }
    //}
    
    float2 sum = (float2)(0.0f,0.0f);
    
    for (int i = 0; i < nn; i++){
        for (int j = 0; j < nn; j++){
            for (int k = 0; k < nn; k++){
                sum += cmplxmult(cl_A[id*nn*nn + i*nn + k],cl_A[id*nn*nn + k*nn + j]);
            }
            cl_AA[id*nn*nn + i*nn + j] = sum;
            sum = (float2)(0.0f,0.0f);
        }
    }
    
    //float2 AAloc[nn][nn];
    //for (int i = 0; i < nn; i++){
    //    for (int j = 0; i < nn; j++){
    //        AAloc[i][j] = cl_AA[id*nn*nn + i*nn + j];
    //    }
   // }

    
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
    
    //float2 AAAloc[nn][nn];
    //for (int i = 0; i < nn; i++){
    //    for (int j = 0; i < nn; j++){
    //        AAAloc[i][j] = cl_AAA[id*nn*nn + i*nn + j];
    //    }
    // }
    
    float2 tmp1,tmp2;
    for (int i = 0; i < nn; i++){
        for (int j = 0; j < nn; j++){
            if(i == j){
                tmp1 = cmplxmult(cl_coeff[0],cl_AA[id*nn*nn + i*nn + j]);
                tmp2 = cmplxmult(cl_coeff[1],cl_A[id*nn*nn + i*nn + j]);
                cl_T1[id*nn*nn + i*nn + j] = cl_AAA[id*nn*nn + i*nn + j] - tmp1 + tmp2 - cl_coeff[2];
            
                tmp1 = cmplxmult(cl_coeff[3],cl_AA[id*nn*nn + i*nn + j]);
                tmp2 = cmplxmult(cl_coeff[4],cl_A[id*nn*nn + i*nn + j]);
                cl_T2[id*nn*nn + i*nn + j] = cl_AAA[id*nn*nn + i*nn + j] - tmp1 + tmp2 - cl_coeff[5];
            
                tmp1 = cmplxmult(cl_coeff[6],cl_AA[id*nn*nn + i*nn + j]);
                tmp2 = cmplxmult(cl_coeff[7],cl_A[id*nn*nn + i*nn + j]);
                cl_T3[id*nn*nn + i*nn + j] = cl_AAA[id*nn*nn + i*nn + j] - tmp1 + tmp2 - cl_coeff[8];
            }
            else {
                tmp1 = cmplxmult(cl_coeff[0],cl_AA[id*nn*nn + i*nn + j]);
                tmp2 = cmplxmult(cl_coeff[1],cl_A[id*nn*nn + i*nn + j]);
                cl_T1[id*nn*nn + i*nn + j] = cl_AAA[id*nn*nn + i*nn + j] - tmp1 + tmp2;
                
                tmp1 = cmplxmult(cl_coeff[3],cl_AA[id*nn*nn + i*nn + j]);
                tmp2 = cmplxmult(cl_coeff[4],cl_A[id*nn*nn + i*nn + j]);
                cl_T2[id*nn*nn + i*nn + j] = cl_AAA[id*nn*nn + i*nn + j] - tmp1 + tmp2;
                
                tmp1 = cmplxmult(cl_coeff[6],cl_AA[id*nn*nn + i*nn + j]);
                tmp2 = cmplxmult(cl_coeff[7],cl_A[id*nn*nn + i*nn + j]);
                cl_T3[id*nn*nn + i*nn + j] = cl_AAA[id*nn*nn + i*nn + j] - tmp1 + tmp2;
            }
        }
    }


    sum = (float2)(0.0f,0.0f);
    for (int i = 0; i < nn; i++){
        for (int j = 0; j < nn; j++){
            for (int k = 0; k < nn; k++){
                sum += cmplxmult(cl_T1[id*nn*nn + i*nn + k],cl_T2[id*nn*nn + k*nn + j]);
            }
            cl_T1T2[id*nn*nn + i*nn + j] = sum;
            sum = (float2)(0.0f,0.0f);
        }
    }

    sum = (float2)(0.0f,0.0f);
    for (int i = 0; i < nn; i++){
        for (int j =0; j < nn; j++){
            for (int k = 0; k < nn ; k++){
                sum += cmplxmult(cl_T1T2[id*nn*nn + i*nn + k],cl_T3[id*nn*nn + k*nn + j]);
            }
            sum /= 362880.0f;
            cl_AA[id*nn*nn + i*nn + j] = sum;
            sum = (float2)(0.0f,0.0f);
        }
    }

    
    for (int l = 0; l < ns; l++){
        sum = (float2)(0.0f,0.0f);
        for (int i = 0; i < nn; i++){
            for (int j = 0; j < nn; j++){
                for (int k = 0; k < nn; k++){
                    sum += cmplxmult(cl_AA[id*nn*nn + i*nn + k],cl_AA[id*nn*nn + k*nn + j]);
                }
                cl_AAA[id*nn*nn + i*nn + j] = sum;
                sum = (float2)(0.0f,0.0f);
            }
        }
        for (int i = 0; i < nn; i++){
            for (int j = 0; j < nn; j++){
                cl_AA[id*nn*nn + i*nn + j] = cl_AAA[id*nn*nn + i*nn + j];
            }
        }
    }
    for (int i = 0; i < nn; i++){
        for (int j = 0; j < nn; j++){
            cl_expA[id*nn*nn + i*nn + j] = cl_AA[id*nn*nn + i*nn + j];
        }
    }
}



