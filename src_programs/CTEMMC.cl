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
// CTEMsoft2013:CTEMMC.cl
//--------------------------------------------------------------------------
//
// PROGRAM: CTEMMC.cl
//
//> @author Saransh Singh, Carnegie Mellon University
//
//> @brief OpenCL kernel for Monte Carlo BSE simulation
//
//> @detail Monte Carlo Electron Trajectory Simulation for EBSD
//>	This version uses the modified Lambert projection to store
//>	the MC output data, so that we are not dependent on a
//>	particular detector geometry.  We store the energy and direction
//>	cosines of a BSE electron along with depth information in the
//>	Lambert projection array, which needs to be sufficiently fine in
//>	terms of sampling so that we can deal with a detector that's relatively
//>	far away.
//> This is the part which runs on the GPU. This is independent of the platform. The data/memory management is done using MainMC.c.
//
//> @date 11/**/12  PGC 1.0 IDL version
//> @date 12/04/12  MDG 1.1 conversion to Fortran-90
//> @date 12/06/12  MDG 1.2 conversion to OpenMP, with new random number generator
//> @date 12/06/12  MDG 1.3 added energy histogram sampling
//> @date 12/07/12  MDG 1.4 added energy vs. depth sampling
//> @date 03/11/13  MDG 2.0 replaced regular storage by modified Lambert projection
//> @date 07/23/13  MDG 3.0 complete rewrite
//> @date 09/25/13  MDG 3.1 modified output file format
//> @date 06/13/14  SS  1.0 OpenCL implementation
//--------------------------------------------------------------------------
#define RAND_MAX  2147483647.0f
#define PI        3.14159f

float *LambertSphereToSquare(float xyz[3]);
int rand_r(unsigned int);
//--------------------------------------------------------------------------
//
// FUNCTION: rand_r
//
//> @author Saransh Singh, Carnegie Mellon University
//
//> @brief function to generate random number
//
//> @details this function generates a random number based on the seed value supplied. the random number generated serves as seed for the subsequent random number generated and so on. this is the the same implementation as in ANSI C
//
//> @positive integer seed
//> @date 05/14/14    SS 1.0 original
//--------------------------------------------------------------------------


int rand_r (unsigned int seed)
{
    unsigned int next = seed;
    int result;
    
    next *= 1103515245;
    next += 12345;
    result = (unsigned int) (next / 65536) % 2048;
    
    next *= 1103515245;
    next += 12345;
    result <<= 10;
    result ^= (unsigned int) (next / 65536) % 1024;
    
    next *= 1103515245;
    next += 12345;
    result <<= 10;
    result ^= (unsigned int) (next / 65536) % 1024;
    
    seed = next;
    
    return result;
}

//--------------------------------------------------------------------------
//
// FUNCTION: LambertSphereToSquare
//
//> @author Saransh Singh, Carnegie Mellon University
//
//> @brief Lambert projection of a point of the unit sphere
//
//> @param xyz 3D coordinates to be considered
//> @date 05/14/14    SS 1.0 original
//--------------------------------------------------------------------------


float *LambertSphereToSquare(float xyz[3]){
    float eps = 1e-4f;
    float q, LPssPi2;
    LPssPi2 = 0.886226925452758f;
    float res[2];
    if (fabs(1.0f-(xyz[0]*xyz[0] + xyz[1]*xyz[1] + xyz[2]*xyz[2])) > eps){
        res[0] = 0.0f;
        res[1] = 0.0f;
    }
    
    else {
        if(fabs(xyz[2]) == 1.0f){
            res[0] = 0.0f;
            res[1] = 0.0f;
        }
        else {
            if (fabs(xyz[1]) <= fabs(xyz[0])){
                q = fabs(xyz[0])/xyz[0] * sqrt(2.0f*(1.0f - xyz[2]));
                res[0] = q * LPssPi2;
                res[1] = q * atan(xyz[1]/xyz[0])/LPssPi2;
            }
            else{
                q = fabs(xyz[1])/xyz[1] * sqrt(2.0f*(1.0f-xyz[2]));
                res[0] = q * atan(xyz[0]/xyz[1])/LPssPi2;
                res[1] = q * LPssPi2;
            }
            
        }
    }
    return res;
    
}

//--------------------------------------------------------------------------
//
// FUNCTION: MC
//
//> @author Saransh Singh, Carnegie Mellon University
//
//> @brief Kernel which does MC simulations for electron scattering
//
//> @details This is the OpenCL kernel which runs on the device and does the simulation. the output variables are Lamx and Lamy (this can only be a void function). The input variables include the energy of the beam, the material parameters and the parameter count and num_max. count is the no. of threads started in each dimension and num_max is the maximum no. of electrons to carry out the simulation for in one run of the kernel. the kernel also needs a prime number to seed the random number generator
//
//> @param Lamx, Lamy, energy of electron, no. of threads, material parameters, maximum electrons for which the simulation is done and a prime number seed
//> @date 05/14/14    SS 1.0 original
//--------------------------------------------------------------------------

__kernel void MC(__global float* Lamx, __global float* Lamy, const float E, const int count, const float z, const float rho, const float A, const int num_max, const int prime, const float sig, const float omega, __global float* depth, __global float* energy)
{
    int tx, ty;
    tx = get_global_id(0);
    ty = get_global_id(1);
    float dir_cos[3];
    float *res;
    
    int seed = count*ty + tx + prime;
	int id = count*ty + tx;
    float rand;
	int rand_seed;
    
    int num;
    num = num_max/(count*count);
    
    int counter1, counter2;
    
    float4 c_new, r_new;
    float E_new, alpha, de_ds, phi, psi, mfp,sig_eNA,step, dsq, dsqi, absc0z;
    
    float J;    // refer to Monte Carlo simulation for Electron Microscopy and Microanalysis, David C. Joy
    J = (9.76f*z + 58.5f*powr(z,-0.19f))*1E-3f;
    
    float4 r0 = (float4)(0.0f, 0.0f, 0.0f, 0.0f);
    float4 c0 = (float4)(cos(sig)*sin(omega), sin(sig)*sin(omega), cos(omega), 0.0f);
    float escape_depth;


// Setting all values to -10. Any value other than -10 will denote a backscattered electron with the x and y component of the Lambert Projection

	for (int i = 0; i < num; ++i){
		Lamx[num*id + i] = -10;
		Lamy[num*id + i] = -10;
        depth[num*id + i] = -10;
        energy[num*id + i] = 0;
	}

    
    
    for (int i = 0; i < num; ++i){
        rand_seed = rand_r(seed);
        seed = rand_seed;
        rand = rand_seed/RAND_MAX; //some random no. generator in gpu
        r0 = (float4)(0.0f, 0.0f, 0.0f, 0.0f);
        c0 = (float4)(cos(sig)*sin(omega), sin(sig)*sin(omega), cos(omega), 0.0f);
        E_new = E;
        //Lamx[num*id + i] = z;
        c_new = c0;
        alpha = (3.4E-3f)*powr(z,0.67f)/E_new;
        sig_eNA = (5.21f * 602.3f)*((z*z)/(E_new*E_new))*((4.0f*PI)/(alpha*(1+alpha)))*((E_new + 511.0f)*(E_new + 511.0f)/((E_new + 1024.0f)*(E_new + 1024.0f)));
        mfp = A/(rho*sig_eNA);
        step = -mfp * log(rand);
        r_new = (float4)(r0.x + step*c_new.x, r0.y + step*c_new.y, r0.z + step*c_new.z, 0.0f);
        r0 = r_new;
        counter1 = 0;   // This is used as a counter for the number of monte carlo steps to carry out for each electron. This is due to the lock step nature of the GPU code. We have arbitly set this to a 1000, though this is material dependent
        
        counter2 = 0;   // This counter is used to figure out if the electron has left the sample or not. Again, this is because of the lock step nature. All steps have to be executed on each thread irrespective of the fact that the electron may have actually left the sample
        
        while (counter1 < 1000){
// inline code rather than function call
// Taken from book Monte Carlo simulation for Electron Microscopy and Microanalysis, David C. Joy

            alpha = (3.4E-3f)*powr(z,0.67f)/E_new;
            sig_eNA = (5.21f * 602.3f)*((z*z)/(E_new*E_new))*((4*PI)/(alpha*(1+alpha)))*((E_new + 511.0f)*(E_new + 511.0f)/((E_new + 1024.0f)*(E_new + 1024.0f)));
            mfp = A/(rho*sig_eNA);
            rand_seed = rand_r(seed);
            seed = rand_seed;
            rand = rand_seed/RAND_MAX; //some random no. generator in gpu
            step = -mfp * log(rand);
// This is the Continuous Slowing Down approximation that we want to get rid of

            de_ds = -78500.0f*(z/(A*E_new)) * log((1.66f*(E_new + 0.85f*J))/J);
            rand_seed = rand_r(seed);
            seed = rand_seed;
            rand = rand_seed/RAND_MAX;
            phi = acos(1 - ((2*alpha*rand)/(1 + alpha - rand)));
            rand_seed = rand_r(seed);
            seed = rand_seed;
            rand = rand_seed/RAND_MAX;
            psi = 2*PI*rand;
            
// new direction cosines of the electrons after scattering event
            if ((c0.z >= 0.999f) || (c0.z <= -0.999f) ){
                absc0z = c0.z;
                absc0z = fabs(absc0z);
                c_new = (float4)(sin(phi) * cos(psi), sin(phi) * sin(psi), (c0.z/absc0z)*cos(phi), 0.0f);
            }
            else {
                dsq = sqrt(1-c0.z*c0.z);
                dsqi = 1/dsq;
                c_new = (float4)(sin(phi)*(c0.x*c0.z*cos(psi) - c0.y*sin(psi))*dsqi + c0.x*cos(phi), sin(phi) * (c0.y * c0.z * cos(psi) + c0.x * sin(psi)) * dsqi + c0.y * cos(phi), -sin(phi) * cos(psi) * dsq + c0.z * cos(phi), 0.0f);
            }
            escape_depth = r_new.z;
            r_new = (float4)(r0.x + step*c_new.x, r0.y + step*c_new.y, r0.z + step*c_new.z, 0.0f);
            r0 = r_new;
            c0 = c_new;
            E_new += step*rho*de_ds;
            if (r0.z <= 0 && counter2 == 0){
                dir_cos[0] = c0.x;
                dir_cos[1] = c0.y;
                dir_cos[2] = c0.z;
                res = LambertSphereToSquare(dir_cos);
                Lamx[num*id + i] = res[0];
                Lamy[num*id + i] = res[1];
                depth[num*id + i] = escape_depth;
                energy[num*id + i] = E_new;
                counter2 = 1;
            }
			
            counter1++ ;
        }
        
        
    }
}