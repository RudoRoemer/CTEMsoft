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

#define PI 3.14159265359f

//--------------------------------------------------------------------------
// EMsoft:Defectmodule.cl
//--------------------------------------------------------------------------

float2 cmplxmult(float2, float2);
float2 conjg(float2);
float2 cmplxexp(float);
float  sqrabs(float2);
float2 cmplxexpTable(float,__global float2*,const int);

//--------------------------------------------------------------------------
// FUNCTION: cmplxmult
//
//> @author Saransh Singh, Carnegie Mellon University
//
//> @brief function to calculate multiplication for two compllex numbers
//
//> @date 09/23/14  SS  1.0 Original
//--------------------------------------------------------------------------

float2 cmplxmult(float2 a, float2 b){
    float2 res;
    res.x = a.x*b.x - a.y*b.y;
    res.y = a.x*b.y + a.y*b.x;
    
    return res;
}

//--------------------------------------------------------------------------
//
// FUNCTION: conjg
//
//> @author Saransh Singh, Carnegie Mellon University
//
//> @brief function to calculate conjugate of a compllex number
//
//> @date 02/19/15  SS  1.0 Original
//--------------------------------------------------------------------------

float2 conjg(float2 a)
{
    float2 ret;
    ret = (float2)(a.x,-a.y);
    return ret;
}

//--------------------------------------------------------------------------
// FUNCTION: cmplxexp
//
//> @author Saransh Singh, Carnegie Mellon University
//
//> @brief function to calculate exponential of imaginary number
//
//> @date 06/15/15  SS  1.0 Original
//--------------------------------------------------------------------------

float2 cmplxexp(float a)
{
    float2 ret;
    ret = (float2)(cos(a),sin(a));
    return ret;
}

//--------------------------------------------------------------------------
// FUNCTION: cmplxexpTable
//
//> @author Saransh Singh, Carnegie Mellon University
//
//> @brief function to calculate exponential of imaginary number using LUT
//
//> @date 06/25/15  SS  1.0 Original
//--------------------------------------------------------------------------

float2 cmplxexpTable(float a,__global float2* table, const int n)
{
    float ap;
    ap = a - floor(a/2.0f/PI)*2.0f*PI;
    for (int i = 0; i < n-1; i++){
        if( (2.0f*PI*i/n - ap)*(2.0f*PI*(i+1)/n - ap) <= 0.0f){
            return (float2)(table[i].x + (table[i+1].x - table[i].x)*(ap*n - i*2.0f*PI)/2.0f/PI,table[i].y + (table[i+1].y - table[i].y)*(ap*n - i*2.0f*PI)/2.0f/PI);
        }
    }
    return (float2)(0.0f,0.0f);
}

//--------------------------------------------------------------------------
// FUNCTION: sqrabs
//
//> @author Saransh Singh, Carnegie Mellon University
//
//> @brief function to calculate absolute suared value of a complex number
//
//> @date 06/17/15  SS  1.0 Original
//--------------------------------------------------------------------------

float sqrabs(float2 a)
{
    float ret;
    ret = a.x*a.x + a.y*a.y;
    return ret;
}

/*
//--------------------------------------------------------------------------
// EMsoft:CTEMDefect.cl
//--------------------------------------------------------------------------
//
// KERNEL: CTEMDefect
//
//> @author Saransh Singh, Carnegie Mellon University
//
//> @brief OpenCL kernel for defect simulation in the CTEM geometry i.e. with
//  a parallel illumination beam. The algorithm closely follows the one 
//  outlined in the paper "Systematic row and zone axis STEM defect image simulations"
//  P.J. Philips, M.j. Mills, M. De Graef Phil. Mag., 2011. The exponential is calculated
//  using the optimized taylor's algorithm.
//
//  @param cl_A list of structure matrices
//  @param cl_AA,Cl_AAA,cl_expA,cl_T1,cl_T2,cl_wavefncoeff,cl_wavefncoeffintd intermediate vars
//  @param arrsize size of each array in cl_A
//  @param cl_coeff complex coefficients for the optimized taylor algorithm
//  @param sqrsize parameter for scalingand squaring algorithm
//  @param offset sum of squares of elements in arrsize till i
//  @param arrsizesum sum of elements in arrsize till i
//  @param defectindex index of a defect in each thickness of a pixel

//> @date 06/08/15  SS  1.0 Original
//--------------------------------------------------------------------------


__kernel void CTEMDefect(__global float2* cl_expA,
                         __global float2* cl_Aperfect,
                         __global float2* cl_Adefect,
                         __global float2* cl_AA,
                         __global float2* cl_AAA,
                         __global int* arrsize,
                         __global float2* cl_coeff,
                         __global float2* cl_T1,
                         __global float2* cl_T2,
                         __global float2* cl_wavefncoeff,
                         __global float2* cl_wavefncoeffintd,
                         __global int* sqrsize,
                         __global int* offset,
                         __global int* arrsizesum,
                         const int numdepth,
                         __global int2* defectindex)
		      
			 
{

    int tx,ty,id;
    tx = get_global_id(0);
    ty = get_global_id(1);
    id = get_global_size(0)*ty + tx;
    
    float2 sum = (float2)(0.0f,0.0f);

    // calculating A^2 and A^3

}*/

//--------------------------------------------------------------------------
// EMsoft:SEMDefect.cl
//--------------------------------------------------------------------------
//
// KERNEL: SEMDefect
//
//> @author Saransh Singh, Carnegie Mellon University
//
//> @brief OpenCL kernel for defect simulation where the list of precomputed
//  scattering matrices is used along with bilinear interpolation to calculate
//  the intensity

//> @date 06/08/15  SS  1.0 Original
//--------------------------------------------------------------------------

__kernel void PreCalcScatMat(   __global float2* cl_A,
                                __global float2* cl_expA,
                                __global float2* cl_AA,
                                __global float2* cl_AAA,
                                __global float2* cl_T1,
                                __global float2* cl_T2,
                                __global float2* cl_coeff,
                                const int nn,
                                const int ns )

{
    int tx, ty, id;
	tx = get_global_id(0);
	ty = get_global_id(1);
	id = get_global_size(0)*ty + tx;
    float scaling;
    scaling = pow(2.0f,ns);
    
    // calculating initial input matrix
    
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
    
    // squaring operation as matrix was scaled
    
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
    
    
}


//--------------------------------------------------------------------------
// EMsoft:CalcScatMatDefects.cl
//--------------------------------------------------------------------------
//
// KERNEL: CalcScatMatDefects
//
//> @author Saransh Singh, Carnegie Mellon University
//
//> @brief OpenCL kernel computing defect imaging using parallel beam illumination
//
//> @date 06/11/15  SS  1.0 Original
//--------------------------------------------------------------------------



__kernel void CalcScatMatDefects(	__global float2* cl_DynMat,
                                    __global float2* cl_ScalFact,
                                    __global float2* cl_A,
                                    __global float2* cl_expA,
                                    __global float2* cl_AA,
                                    __global float2* cl_AAA,
                                    __global float2* cl_T1,
                                    __global float2* cl_T2,
                                    __global float2* cl_wavefncoeff,
                                    __global float2* cl_wavefncoeffintd,
                                    __global float2* cl_coeff,
                                    __global float2* cl_gdotR,
                                    const int nn,
                                    const int ns,
                                    const int numdepth)


{
	int tx, ty, id;
	tx = get_global_id(0);
	ty = get_global_id(1);
	id = get_global_size(0)*ty + tx;
	int offset;
    float scaling;
	float2 gdotR,arg1,arg2,DynMatScald;
    scaling = pow(2.0f,ns);
    
    for (int i = 0; i < nn; i++){
        if ( i == 0) {
            cl_wavefncoeff[id*nn + i] = (float2)(1.0f,0.0f);
        }
        else {
            cl_wavefncoeff[id*nn + i] = (float2)(0.0f,0.0f);
        }
        
    }

	for (int nstep = 0; nstep < numdepth; nstep++){
	
// calculating initial input matrix
    
        offset = nstep*get_global_size(0)*get_global_size(1);
        gdotR = cl_gdotR[offset + id];
 
        for (int i = 0; i < nn; i++){
            for (int j = 0; j < nn; j++){
                arg1 = cmplxexp(-2.0f*PI*cl_ScalFact[i*nn+j].x*gdotR.x);
                arg2 = cmplxexp(-2.0f*PI*cl_ScalFact[i*nn+j].y*gdotR.y);
                arg1 = cmplxmult(arg1,arg2);
                DynMatScald = cl_DynMat[id*nn*nn + i*nn + j]/scaling;
		cl_A[id*nn*nn + i*nn + j] = cmplxmult(DynMatScald,arg1); 

            }
        }
	//barrier(CLK_GLOBAL_MEM_FENCE);					

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
    
// squaring operation as matrix was scaled
    
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
        
        sum = (float2)(0.0f,0.0f);
        for (int l = 0; l < nn; l++){
            cl_wavefncoeffintd[id*nn + l] = cl_wavefncoeff[id*nn + l];
        }
        for (int j = 0; j < nn; j++){
            for (int k = 0; k < nn; k++){
                sum += cmplxmult(cl_expA[id*nn*nn + j*nn + k],cl_wavefncoeffintd[id*nn + k]);
            }
            cl_wavefncoeff[id*nn + j] = sum;
            sum = (float2)(0.0f,0.0f);
        }

// we have the fourier coefficients of the wavefunction here
// just take the square modulus to get the intensity
        
	}
    barrier(CLK_GLOBAL_MEM_FENCE);
    
}
                                

//--------------------------------------------------------------------------
// EMsoft:CalcScatMatDefectsTable.cl
//--------------------------------------------------------------------------
//
// KERNEL: CalcScatMatDefectsTable
//
//> @author Saransh Singh, Carnegie Mellon University
//
//> @brief OpenCL kernel for defect imaging for parallel beam illumination using a precomputed list of sin and cosine.
//
//> @date 06/16/15  SS  1.0 Original
//--------------------------------------------------------------------------



__kernel void CalcScatMatDefectsTable(	__global float2* cl_DynMat,
                                 __global float2* cl_ScalFact,
                                 __global float2* cl_A,
                                 __global float2* cl_expA,
                                 __global float2* cl_AA,
                                 __global float2* cl_AAA,
                                 __global float2* cl_T1,
                                 __global float2* cl_T2,
                                 __global float2* cl_wavefncoeff,
                                 __global float2* cl_wavefncoeffintd,
                                 __global float2* cl_coeff,
                                 __global float2* cl_gdotR,
                                 __global float2* cl_table,
                                 const int nn,
                                 const int ns,
                                 const int numdepth,
                                 const int numentries)


{
	int tx, ty, id;
	tx = get_global_id(0);
	ty = get_global_id(1);
	id = get_global_size(0)*ty + tx;
	int offset;
    float scaling;
	float2 gdotR,arg1,arg2,DynMatScald;
    scaling = pow(2.0f,ns);
    
    for (int i = 0; i < nn; i++){
        if ( i == 0) {
            cl_wavefncoeff[id*nn + i] = (float2)(1.0f,0.0f);
        }
        else {
            cl_wavefncoeff[id*nn + i] = (float2)(0.0f,0.0f);
        }
        
    }
    
	for (int nstep = 0; nstep < numdepth; nstep++){
        
        // calculating initial input matrix
        
        offset = nstep*get_global_size(0)*get_global_size(1);
        gdotR = cl_gdotR[offset + id];
        
        for (int i = 0; i < nn; i++){
            for (int j = 0; j < nn; j++){
                arg1 = cmplxexpTable(-2.0f*PI*cl_ScalFact[i*nn+j].x*gdotR.x,cl_table,numentries);
                arg2 = cmplxexpTable(-2.0f*PI*cl_ScalFact[i*nn+j].y*gdotR.y,cl_table,numentries);
                arg1 = cmplxmult(arg1,arg2);
                DynMatScald = cl_DynMat[id*nn*nn + i*nn + j]/scaling;
                cl_A[id*nn*nn + i*nn + j] = cmplxmult(DynMatScald,arg1);
                
            }
        }
        //barrier(CLK_GLOBAL_MEM_FENCE);
        
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
        
        // squaring operation as matrix was scaled
        
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
        
        sum = (float2)(0.0f,0.0f);
        for (int l = 0; l < nn; l++){
            cl_wavefncoeffintd[id*nn + l] = cl_wavefncoeff[id*nn + l];
        }
        for (int j = 0; j < nn; j++){
            for (int k = 0; k < nn; k++){
                sum += cmplxmult(cl_expA[id*nn*nn + j*nn + k],cl_wavefncoeffintd[id*nn + k]);
            }
            cl_wavefncoeff[id*nn + j] = sum;
            sum = (float2)(0.0f,0.0f);
        }
        
        // we have the fourier coefficients of the wavefunction here
        // just take the square modulus to get the intensity
        
	}
    
}
                                           


