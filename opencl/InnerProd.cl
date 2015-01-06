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

#define BLOCK_SIZE 16
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