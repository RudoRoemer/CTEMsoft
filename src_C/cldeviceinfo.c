//
//  clDeviceInfo.c
//  
//
//  Created by Saransh Singh on 4/25/14.
//
//

#include <fcntl.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <OpenCL/opencl.h>
#include <time.h>

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv)
{
    system("clear");
    int err;                            //error code returned from API call
    cl_uint num_platforms;
    cl_platform_id platforms;
    cl_ulong mem_param;
    cl_uint param;
    cl_uint NumDevices;
    size_t param1;
    size_t param2[3];
    char cldeviceName[256];
    
    int i = 0;
    
    err = clGetPlatformIDs(2, &platforms, &num_platforms);
    
    if (num_platforms == 0)
    {
        printf("No platforms found\n");
        exit(0);
    }
    
    printf("Number of platform = %d\n",num_platforms);
    
//////////////////////////////////////////////////////////////////////////////
    
    err = clGetDeviceIDs(NULL, CL_DEVICE_TYPE_GPU, 3, NULL, &NumDevices);
        if (err != CL_SUCCESS)
        {
                printf("Error: Failed to create a device group!\n");
        printf("No GPU detected\n");
                return EXIT_FAILURE;
        }
    printf("GPU detected on the machine. No. of GPU's = %d\n", NumDevices);
    
    cl_device_id device_id[NumDevices];
    err = clGetDeviceIDs(NULL, CL_DEVICE_TYPE_GPU, NumDevices, device_id, &NumDevices);
//////////////////////////////////////////////////////////////////////////////
    
    while (i < NumDevices){
        

        err = clGetDeviceInfo(device_id[i], CL_DEVICE_NAME, sizeof(cldeviceName), cldeviceName, NULL);
        if (err != CL_SUCCESS)
        {
            printf("Error: Failed to gather info about the address bit!\n");
            return EXIT_FAILURE;
        }
        printf("Name of device %d is %s\n", i+1 , cldeviceName);
        
/////////////////////////////////////////////////////////////////////////////

    err = clGetDeviceInfo(device_id[i], CL_DEVICE_ADDRESS_BITS, sizeof(param), &param, NULL);
    if (err != CL_SUCCESS)
        {
                printf("Error: Failed to gather info about the address bit!\n");
                return EXIT_FAILURE;
        }
    printf("Your compute device is a %d bit device\n", param);
    
/////////////////////////////////////////////////////////////////////////////
    
    err = clGetDeviceInfo(device_id[i], CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(param), &param, NULL);
    if (err != CL_SUCCESS)
        {
                printf("Error: Failed to find maximum no. of compute units!\n");
                return EXIT_FAILURE;
        }
    printf("No. of parallel compute core in the device is %d\n", param);
    
/////////////////////////////////////////////////////////////////////////////


    err = clGetDeviceInfo(device_id[i], CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(mem_param), &mem_param, NULL);
    if (err != CL_SUCCESS)
        {
                printf("Error: Failed to get global memory size!\n");
                return EXIT_FAILURE;
        }
    printf("Global memory size of device is %llu Mega bytes\n", mem_param/(1024*1024));
    
/////////////////////////////////////////////////////////////////////////////
    
    err = clGetDeviceInfo(device_id[i], CL_DEVICE_LOCAL_MEM_SIZE, sizeof(mem_param), &mem_param, NULL);
    if (err != CL_SUCCESS)
        {
                printf("Error: Failed to get local memory size!\n");
                return EXIT_FAILURE;
        }
    printf("Local memory size of device is %llu Kilo bytes\n", mem_param/(1024));

/////////////////////////////////////////////////////////////////////////////
    
    err = clGetDeviceInfo(device_id[i], CL_DEVICE_MAX_CLOCK_FREQUENCY, sizeof(param), &param, NULL);
    if (err != CL_SUCCESS)
        {
                printf("Error: Failed to get maximum clock speed !\n");
                return EXIT_FAILURE;
        }
    printf("Maximum clock speed of the device is %d MHz\n", param);
    
/////////////////////////////////////////////////////////////////////////////
    
    
    err = clGetDeviceInfo(device_id[i], CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(param1), &param1, NULL);
    if (err != CL_SUCCESS)
        {
                printf("Error: Failed to get maximum work group size!\n");
                return EXIT_FAILURE;
        }
    printf("Maximum work group size of the device is %zu\n", param1);
    
/////////////////////////////////////////////////////////////////////////////
    
    
    err = clGetDeviceInfo(device_id[i], CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS, sizeof(param), &param, NULL);
    if (err != CL_SUCCESS)
        {
                printf("Error: Failed to get maximum work item dimensions!\n");
                return EXIT_FAILURE;
        }
    printf("Maximum dimensions to specify local and global work items are %d\n", param);
    
/////////////////////////////////////////////////////////////////////////////
    
    
    err = clGetDeviceInfo(device_id[i], CL_DEVICE_MAX_WORK_ITEM_SIZES, sizeof(param2), param2, NULL);
    if (err != CL_SUCCESS)
        {
                printf("Error: Failed to get maximum work item size!\n");
                return EXIT_FAILURE;
        }
    printf("Maximum no. of work items that can be specified in each dimension are [%zu %zu %zu]\n", param2[0], param2[1], param2[2]);
    
/////////////////////////////////////////////////////////////////////////////
        i++;
        printf("////////////////////////////////////////////////////////////\n");
        printf("////////////////////////////////////////////////////////////\n\n");
    }
    return 0;
}

