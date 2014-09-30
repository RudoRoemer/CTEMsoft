//
//  MatrixExpo.c
//  
//
//  Created by Saransh Singh on 9/9/14.
//
//

#define NUM (40)
#define MAX_SOURCE_SIZE (0x100000)


#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <OpenCL/opencl.h>
#include <time.h>


cl_float2 cmplxmult(cl_float2 a, cl_float2 b){
    cl_float2 res;
    res.s[0] = a.s[0]*b.s[0] - a.s[1]*b.s[1];
    res.s[1] = a.s[0]*b.s[1] + a.s[1]*b.s[0];
    
    return res;
}


int main(int argc, char** argv){
    int err;                            // error code returned from api calls
    int ns = 1;
    unsigned int count = NUM;
    unsigned int threadcnt = 64;
    cl_ulong time_start, time_end;
    double total_time = 0;
    cl_float2 coeff[9],cvals[9],tmp;
    cl_float2 *A;
    cl_float2 *expAres;            // original data set given to device
    A = malloc(sizeof(cl_float2)*count*count* threadcnt * threadcnt);
    expAres = malloc(sizeof(cl_float2)*count*count* threadcnt * threadcnt);

    size_t global[2] ;                      // global domain size for our calculation
    size_t local[2] ;                       // local domain size for our calculation

    cl_device_id device_id;             // compute device id
    cl_context context;                 // compute context
    cl_command_queue commands;          // compute command queue
    cl_program program;                 // compute program
    cl_kernel kernel;                   // compute kernel
    cl_event event;

    cl_mem cl_expA;
    cl_mem cl_A;
    cl_mem cl_AA;
    cl_mem cl_AAA;
    cl_mem cl_T1,cl_T2,cl_T3,cl_T1T2,cl_T1T2T3,cl_sqr,cl_coeff;
    
    cvals[0].s[0] = -3.3335514852690488032942739163345055;
    cvals[0].s[1] = 0.0;
    
    cvals[1].s[0] = -3.0386480729366970892124687564926859;
    cvals[1].s[1] = -1.5868011957588383288038677051222921;
    
    cvals[2].s[0] = -3.0386480729366970892124687564926859;
    cvals[2].s[1] = 1.5868011957588383288038677051222921;
    
    cvals[3].s[0] = -2.1108398003026547374987047865183922;
    cvals[3].s[1] = -3.0899109287255009227777015426228801;
    
    cvals[4].s[0] = -2.1108398003026547374987047865183922;
    cvals[4].s[1] = 3.0899109287255009227777015426228801;
    
    cvals[5].s[0] = -0.38106984566311299903129424501333242;
    cvals[5].s[1] = -4.3846445331453979503692027283066828;
    
    cvals[6].s[0] = -0.38106984566311299903129424501333242;
    cvals[6].s[1] = 4.3846445331453979503692027283066828;
    
    cvals[7].s[0] = 2.6973334615369892273896047461916633;
    cvals[7].s[1] = -5.1841620626494141778340870727109629;
    
    cvals[8].s[0] = 2.6973334615369892273896047461916633;
    cvals[8].s[1] = 5.1841620626494141778340870727109629;
    
    
    coeff[0].s[0] = cvals[0].s[0]+cvals[1].s[0]+cvals[2].s[0];
    coeff[0].s[1] = cvals[0].s[1]+cvals[1].s[1]+cvals[2].s[1];

    coeff[1].s[0] = cmplxmult(cvals[0],cvals[1]).s[0] + cmplxmult(cvals[1],cvals[2]).s[0] + cmplxmult(cvals[0],cvals[2]).s[0];
    coeff[1].s[1] = cmplxmult(cvals[0],cvals[1]).s[1] + cmplxmult(cvals[1],cvals[2]).s[1] + cmplxmult(cvals[0],cvals[2]).s[1];
    
    tmp = cmplxmult(cvals[1],cvals[2]);
    coeff[2] = cmplxmult(cvals[0],tmp);
    
    
    coeff[3].s[0] = cvals[3].s[0]+cvals[4].s[0]+cvals[5].s[0];
    coeff[3].s[1] = cvals[3].s[1]+cvals[4].s[1]+cvals[5].s[1];
    
    coeff[4].s[0] = cmplxmult(cvals[3],cvals[4]).s[0] + cmplxmult(cvals[4],cvals[5]).s[0] + cmplxmult(cvals[3],cvals[5]).s[0];
    coeff[4].s[1] = cmplxmult(cvals[3],cvals[4]).s[1] + cmplxmult(cvals[4],cvals[5]).s[1] + cmplxmult(cvals[3],cvals[5]).s[1];

    tmp = cmplxmult(cvals[4],cvals[5]);
    coeff[5] = cmplxmult(cvals[3],tmp);
    
    coeff[6].s[0] = cvals[6].s[0]+cvals[7].s[0]+cvals[8].s[0];
    coeff[6].s[1] = cvals[6].s[1]+cvals[7].s[1]+cvals[8].s[1];
    
    coeff[7].s[0] = cmplxmult(cvals[6],cvals[7]).s[0] + cmplxmult(cvals[7],cvals[8]).s[0] + cmplxmult(cvals[6],cvals[8]).s[0];
    coeff[7].s[1] = cmplxmult(cvals[6],cvals[7]).s[1] + cmplxmult(cvals[7],cvals[8]).s[1] + cmplxmult(cvals[6],cvals[8]).s[1];

    tmp = cmplxmult(cvals[7],cvals[8]);
    coeff[8] = cmplxmult(cvals[6],tmp);
    
    FILE *fp;
	const char fileName[] = "CTEMMatExp.cl";
	size_t source_size;
	char *source_str;
    
	/* Load kernel source file */
	fp = fopen(fileName, "rb");
	if (!fp) {
		fprintf(stderr, "Failed to load kernel.\n");
		exit(1);
	}
	source_str = (char *)malloc(MAX_SOURCE_SIZE);
	source_size = fread(source_str, 1, MAX_SOURCE_SIZE, fp);
	fclose(fp);

    int nt;
    for(int i = 0; i < count*count* threadcnt * threadcnt; i++){                 // first input value
        A[i].s[0] = 1.0;//(rand() / (float)RAND_MAX)*0.01;
        A[i].s[1] = 1.0;//(rand() / (float)RAND_MAX)*0.01;
    }
    //ns = 1;
    ns = 7;
    float fact = pow(2.0,ns);

    for(int i = 0; i < count*count* threadcnt * threadcnt; i++){                 // first input value
        A[i].s[0] /= fact;//(rand() / (float)RAND_MAX)*0.01;
        A[i].s[1] /= fact;//(rand() / (float)RAND_MAX)*0.01;
    }
    err = clGetDeviceIDs(NULL, CL_DEVICE_TYPE_GPU, 1, &device_id, NULL);
    if (err != CL_SUCCESS)
    {
        printf("Error: Failed to create a device group!\n");
        return EXIT_FAILURE;
    }

    // Create a compute context
    //
    context = clCreateContext(NULL, 1, &device_id, NULL, NULL, &err);
    if (!context)
    {
        printf("Error: Failed to create a compute context!\n");
        return EXIT_FAILURE;
    }

    commands = clCreateCommandQueue(context, device_id, CL_QUEUE_PROFILING_ENABLE, &err);
    if (!commands)
    {
        printf("Error: Failed to create a command commands!\n");
        return EXIT_FAILURE;
    }

    program = clCreateProgramWithSource(context, 1, (const char **) & source_str, NULL, &err);
    if (!program)
    {
        printf("Error: Failed to create compute program!\n");
        return EXIT_FAILURE;
    }
    
    // Build the program executable
    //
    err = clBuildProgram(program, 1, &device_id, NULL, NULL, NULL);
    if (err != CL_SUCCESS)
    {
        size_t len;
        char buffer[source_size];
        
        printf("Error: Failed to build program executable!\n");
        clGetProgramBuildInfo(program, device_id, CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &len);
        printf("%s\n", buffer);
        exit(1);
    }
    
    // Create the compute kernel in the program we wish to run
    //
    kernel = clCreateKernel(program, "MatExpImg", &err);
    if (!kernel || err != CL_SUCCESS)
    {
        printf("Error: Failed to create compute kernel!\n");
        exit(1);
    }

    // Create the input and output arrays in device memory for our calculation
    //
    cl_expA = clCreateBuffer(context,  CL_MEM_WRITE_ONLY,  sizeof(cl_float2) * count * count * threadcnt * threadcnt, NULL, NULL);
    cl_A = clCreateBuffer(context,  CL_MEM_READ_WRITE,  sizeof(cl_float2) * count * count * threadcnt * threadcnt, NULL, NULL);
    cl_AA = clCreateBuffer(context,  CL_MEM_READ_WRITE,  sizeof(cl_float2) * count * count * threadcnt * threadcnt, NULL, NULL);
    cl_AAA = clCreateBuffer(context,  CL_MEM_READ_WRITE,  sizeof(cl_float2) * count * count * threadcnt * threadcnt, NULL, NULL);
    cl_T1 = clCreateBuffer(context,  CL_MEM_READ_WRITE,  sizeof(cl_float2) * count * count * threadcnt * threadcnt, NULL, NULL);
    cl_T2 = clCreateBuffer(context,  CL_MEM_READ_WRITE,  sizeof(cl_float2) * count * count * threadcnt * threadcnt, NULL, NULL);
    cl_T3 = clCreateBuffer(context,  CL_MEM_READ_WRITE,  sizeof(cl_float2) * count * count * threadcnt * threadcnt, NULL, NULL);
    cl_T1T2 = clCreateBuffer(context,  CL_MEM_READ_WRITE,  sizeof(cl_float2) * count * count * threadcnt * threadcnt, NULL, NULL);
    cl_T1T2T3 = clCreateBuffer(context,  CL_MEM_READ_WRITE,  sizeof(cl_float2) * count * count * threadcnt * threadcnt, NULL, NULL);
    cl_sqr = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(cl_float2) * count * count * threadcnt * threadcnt, NULL, NULL);
    cl_coeff = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(cl_float2) * 9, NULL, NULL);
    
    if (!cl_expA || !cl_A || !cl_AA || !cl_AAA || !cl_T1 || !cl_T2 || !cl_T3 || !cl_coeff || !cl_T1T2 || !cl_T1T2T3 || !cl_sqr)
    {
        printf("Error: Failed to allocate device memory!\n");
        exit(1);
    }
    
    
    // Write our data set into the input array in device memory
    //
    err = clEnqueueWriteBuffer(commands, cl_A, CL_TRUE, 0, sizeof(cl_float2) * count * count * threadcnt * threadcnt, A, 0, NULL, NULL);
    if (err != CL_SUCCESS)
    {
        printf("Error: Failed to write to source array!\n");
        exit(1);
    }
    
    err = clEnqueueWriteBuffer(commands, cl_coeff, CL_TRUE, 0, sizeof(cl_float2)*9, coeff, 0, NULL, NULL);
    if (err != CL_SUCCESS)
    {
        printf("Error: Failed to write to source array!\n");
        exit(1);
    }
    //err = clEnqueueReadBuffer( commands, cl_A, CL_TRUE, 0, sizeof(cl_float2) * count * count * threadcnt * threadcnt, expAres, 0, NULL, NULL );
    //if (err != CL_SUCCESS)
    // {
    //    printf("Error: Failed to read output array! %d\n", err);
    //    exit(1);
    //}
    //for (int i = 5000; i < 5010; i++) {
    //    printf("%f\t%f\n",expAres[i].s[0],expAres[i].s[1]);
    // }

    // Set the arguments to our compute kernel
    //
    err = 0;
    err   =  clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *)&cl_expA);
    err  |= clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *)&cl_A);
    err  |= clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *)&cl_AA);
    err  |= clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *)&cl_AAA);
    err  |= clSetKernelArg(kernel, 4, sizeof(unsigned int), (void *)&count);
    err  |= clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *)&cl_coeff);
    err  |= clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *)&cl_T1);
    err  |= clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *)&cl_T2);
    err  |= clSetKernelArg(kernel, 8, sizeof(cl_mem), (void *)&cl_T3);
    err  |= clSetKernelArg(kernel, 9, sizeof(cl_mem), (void *)&cl_T1T2);
    err  |= clSetKernelArg(kernel, 10, sizeof(cl_mem), (void *)&cl_T1T2T3);
    err  |= clSetKernelArg(kernel, 11, sizeof(cl_mem), (void *)&cl_sqr);
    err  |= clSetKernelArg(kernel, 12, sizeof(unsigned int), (void *)&ns);
    err  |= clSetKernelArg(kernel, 13, sizeof(unsigned int), (void *)&threadcnt);


    if (err != CL_SUCCESS)
    {
        printf("Error: Failed to set kernel arguments! %d\n", err);
        exit(1);
    }

    global[0] = threadcnt;
    global[1] = threadcnt;
    local[0] = 8;
    local[1] = 8;

    err = clEnqueueNDRangeKernel(commands, kernel, 2, NULL, global, local, 0, NULL, &event);
    if (err)
    {
        printf("Error: Failed to execute kernel!\n");
        return EXIT_FAILURE;
    }
    // Wait for the command commands to get serviced before reading back results
    //
    clFinish(commands);

    clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_START, sizeof(time_start), &time_start, NULL);
    clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_END, sizeof(time_end), &time_end, NULL);
    total_time = time_end - time_start;
    

    //printf("%f sec for GPU using cl profiler\n",total_time*1E-9 );
    
    err = clEnqueueReadBuffer( commands, cl_expA, CL_TRUE, 0, sizeof(cl_float2) * count * count * threadcnt * threadcnt, expAres, 0, NULL, NULL );
    if (err != CL_SUCCESS)
    {
        printf("Error: Failed to read output array! %d\n", err);
        exit(1);
    }
    int val = 0;
    for(int i = 0; i < count*count*threadcnt*threadcnt; i++){
        //printf("%0.10f\t%0.10f\n",expAres[i].s[0],expAres[i].s[1]);
        if ((expAres[i].s[0] == 0.0) || (expAres[i].s[1] == 0.0)) {
            val++;
        }
    }
    if (val > 0) {
        printf("Something wrong");
    }
    

    clReleaseMemObject(cl_expA);
    clReleaseMemObject(cl_A);
    clReleaseMemObject(cl_AA);
    clReleaseMemObject(cl_AAA);
    clReleaseMemObject(cl_T1);
    clReleaseMemObject(cl_T2);
    clReleaseMemObject(cl_T3);
    clReleaseMemObject(cl_coeff);
    clReleaseProgram(program);
    clReleaseKernel(kernel);
    clReleaseCommandQueue(commands);
    clReleaseContext(context);

}