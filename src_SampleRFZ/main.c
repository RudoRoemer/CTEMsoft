//
//  main.c
//  
//
//  Created by Saransh Singh on 6/17/14.
//
//

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

//--------------------------------------------------------------------------
//
// PROGRAM: so3sampler
//
//> @author Saransh Singh, Carnegie Mellon University
//
//> @brief Example program illustrating how to use the Rodrigues fundamental zone sampler
//
//> @date 5/30/14   SS 1.0 original
//--------------------------------------------------------------------------

int main(int argc, char **argv){
    
    int i, pgnum, nsteps;
    struct SampleRFZ_return data;
    struct FZpointd *FZtmp;
    
    printf("Enter the point group number:\t");
    scanf("%d",&pgnum);
    
    if (pgnum > 32 || pgnum < 1) {
        printf("There are only 32 point groups!\n Invalid point group entered.\n");
        exit(0);
    }
    
    printf("Enter the number of intervals along the cube semi-edge length:\t");
    scanf("%d",&nsteps);
    
    // get the linked list for the FZ for point group symmetry pgnum for 100 steps along the cubic semi-edge
    data = SampleRFZ(nsteps, pgnum);
    printf("Total number of orientations inside Rodrigues Fundamental Zone = %d\n", data.FZcnt);
    
    // now we have the linked list so we can do anything we want with it.
    // in this test program, we create a VTK file so that we can visualize the RFZ with ParaView
    
    FILE *fp;
    
    fp = fopen("Test_Saransh.vtk","w");
    fprintf(fp,"#vtk DataFile Version 2.0\n");
    fprintf(fp,"Uniform sampling of Rodrigues FZ\n");
    fprintf(fp,"ASCII\n");
    fprintf(fp,"DATASET POLYDATA\n");
    fprintf(fp,"POINTS %d float\n",data.FZcnt);
    
    // scan through the list
    FZtmp = data.FZpointer;
    for (int i = 0; i < data.FZcnt; i++) {
        fprintf(fp,"%0.4f\t %0.4f\t %0.4f\n", FZtmp->rod[0], FZtmp->rod[1], FZtmp->rod[2]);
        FZtmp = FZtmp->next;
    }
    
    // alternate way to print to file
    
    /*while (FZtmp != NULL) {
     fprintf(fp,"%0.4f \t %0.4f \t %0.4f\n", FZtmp->rod[0], FZtmp->rod[1], FZtmp->rod[2]);
     FZtmp = FZtmp->next;
     }*/
    
    fprintf(fp, " \n");
    fprintf(fp, "POINT_DATA %d\n",data.FZcnt);
    fprintf(fp, "SCALARS radii float 1\n");
    fprintf(fp, "LOOKUP_TABLE default\n");
    
    double a = 0.1;
    for (int i = 0; i < data.FZcnt; i++) {
        fprintf(fp, "%0.4f\n",a);
    }
    
    return 0;
    
}
