

//##########################################
    //File name: CalcP.c
    //Author: Harrison Agrusa
    //Date created:
    //Date last modified:
    //Description: calculate mutual potential between rubble piles
    //Input: command line text file with 4 columns: x,y,z, mass, color (to distinguish bodies) in any (self-consistent)units
    //Output: gravitational potential (scalar) (in form of txt file)
    //WARNING: This code uses whatever units that its given and does not multiply by G
    //        YOU MUST multiply the output answer by G in the correct units 
//##########################################

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <time.h>
//include <omp.h>

//A function to display an error message and then exit
void fatal(char *message) {
    char error_message[100];
    strcpy(error_message, "[!!] Fatal Error: ");
    strncat(error_message, message, 83); //83 = 100 - 17, max number of characters to append
    perror(error_message);
    exit(-1);
}

int safe_alloc(void * p) {
    if (!p) {
        printf("ALLOCATION ERROR\n");
        return -1;
    } else {
        return 0;
    }
}

int main(int argc, char ** argv) {
    char *fname = (char *) malloc(50); //text files cannot be more than 50 characters
    char *outname = (char *) malloc(50); //output file for energy (just going to be 1 number)
    strcpy(fname, argv[1]);
    strcpy(outname, argv[2]);
    FILE *file;
    FILE *outfile;

    file = fopen(fname, "r");

    int n = 0, i, j;
    int start_index = 0; //by default
    double dummy[5];

    if (argc != 3) {
        printf("Incorrect number of command line arguments.\n");
        printf("Usage: ./CalcPotential data.txt outname.txt\n");
        exit(1);
    }

    if (file == NULL) {
        fatal("Cannot open file");
    }
    //find number of particles
    while(!feof(file)) {
        if (fscanf(file, "%lf %lf %lf %lf %lf",
            dummy + 0, dummy + 1, dummy + 2, dummy + 3, dummy + 4) == 5) {
            n++;
        }
    }
    //printf("    Number of Particles:            %d\n", n);
    //printf("\nReading Input File (%s)......\n", fname);
    free(fname);
    rewind(file); //go back to begining on input file, start reading
    //initialize vectors
    double * mass;
    double * x;
    double * y;
    double * z;
    int * color;
    mass = (double *) malloc(n * sizeof(double));
    x = (double *) malloc(n * sizeof(double));
    y = (double *) malloc(n * sizeof(double));
    z = (double *) malloc(n * sizeof(double));
    color = (int *) malloc(n * sizeof(int));

    if (safe_alloc(mass) + safe_alloc(x) + safe_alloc(z) + safe_alloc(z)  + safe_alloc(color) < 0) {
        free(mass);
        free(x);
        free(y);
        free(z);
        free(color);
        return -1;
    }

    int error_count = n;
    for(i = 0; i < n; i++) {
        if (fscanf(file, "%lf %lf %lf %lf %d", mass + i,
            x + i,
            y + i,
            z + i,
            color + i) == 5) {
            error_count--;
        }
    }

    //printf("Successfully read particles with %d errors\n", error_count);
    fclose(file);
    // Test the nbody tree!
    //printf("%d\n", color[0]);
    //printf("\nCalculating Potential....\n");
    //clock_t start, end;
    double cpu_time_used;
    //start = clock();
    double current_potential;
    double r;
    outfile = fopen(outname, "w");

    //do computation
    for (i=0; i < n; i++) {
        //printf("%d\n", color[i]);
        current_potential = 0.0;
        for (j=0; j < n; j++) {
            if (i != j) {
                r = sqrt(   pow(x[i]-x[j], 2.0) +
                            pow(y[i]-y[j], 2.0) +
                            pow(z[i]-z[j], 2.0)
                        );
                current_potential += mass[i]*mass[j]/r;
            }
        }
        fprintf(outfile, "%lf \n", -current_potential);
    }

    //end = clock();
    //cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    //printf("\nFinished calculation %d particles. \n", n);
    //printf("\nTotal Computation Time: %f seconds\n", cpu_time_used);
    //printf("\nTotal Gravitational Potential: %f \n",E_grav);
    //write total potential to text file
    fclose(outfile);
    free(mass);
    free(x);
    free(y);
    free(z);
    free(color);
    return 0;
}
