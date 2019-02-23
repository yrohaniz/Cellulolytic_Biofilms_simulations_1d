/* This code runs the simulation in one dimension along a line of grid cells.
 * Here the number of attachment events is kept the same and determined by the
 * coarsest grid size. The only change is in the grid size and the starting seed.
 * The random number generation is performed in parallel and as the grid becomes
 * finer the spatial placement of attached bacteria is propagated through cells
 * to match that of the coarsest grid size.*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "CeBiSi.h"
#include <omp.h>
#include <time.h>
#include <sys/stat.h>
#include <string.h>

//****************************************************************************************
//DECLARATION OF CONSTANT GLOBAL AND GLOBAL VARIABLES
//****************************************************************************************
unsigned long NTHRDS; //The number of threads to be used for computation
double INIT_T, END_T; //Start and stop times of the simulation
unsigned long long SEED_init; //This variable is used to initialize the process of random number generation at the beginning of each simulation
unsigned long long SEED_fin; //This is the final value of seed, the simulation steps through the seed values in unit size increments
double INIT_M, INIT_C; //Initial biomass density and carbon concentration <--[0,1]
double INIT_PHI; //Initial attachment function value
unsigned long BMASS_INIT_LOC; //Initial biomass deposit location
double DOMAIN_L; //The length of the domain of simulation
double LOWER_BMASS; //The lower-bound of biomass density used for determining the time step size for a particular grid size
double UPPER_BMASS; //The upper-bound of biomass density used for limiting the time step size reduction due to time adaptive method
double T_PRTN; //Time increments used to record data
double BMASS_ATTACH_RATE; //This is the value of the density of attached biomass per unit time
unsigned long N_coarse; //Coarsest grid size
unsigned long N_fine; //Finest grid size, note that simulation uses powers of two as its grid sizes
unsigned long N; //Refinement of the simulation domain in x direction
#define NUM_NNS 2 //The number of nearest neighbors to each grid cell in the 1d case
const int LT=0, RT=1; //Left and right directions with respect to a grid cell
#define PADS 32 //This number is used to pad the arrays in parallel region for avoiding false sharing
long t_num_invs; //The number of the time intervals between start and stop times
double step; //The time step size
double dlt_x, dlt_x_2; //The cell size in 1d and its square
double inv_dlt_x; //The inverse of the cell size in 1d
double sqrtstep; //The square root of the time step used in the Euler-Maruyama recipe
//****************************************************************************************
//****************************************************************************************

static inline double min(double a, double b) {
    if (a > b)
        return b;
    return a;
}//This function is used to find the lesser of two floating point vars

int main() {

    //*************************************************************************************
    //VARIABLE DECLARATION FOR MAIN
    //*************************************************************************************
    FILE *f_out1, *f_out2; //Output files for biomass density and carbon concentration values respectively
    FILE *f_out3, *f_out4; //Output files for the spatially integrated biomass density and spatially integrated carbon concentration
    FILE *f_out5, *f_out6; //Output file for the simulation log and the attachment function phi
    unsigned long long **random_1, **random_2; //These are random variables that seed the random number generation of Marsaglia_polar_generator
    unsigned long long mult2_n, mult3_n, inc2_n, inc3_n, v;//The multipliers and increments used in lcg2 and lcg 3, v is seed counter
    unsigned long i, j, l, q, nthrds, id; //Loop counters, nthrds is number of threads, id for thread identification
    double g1, g2, **arr_g; //Declaration of random Gaussian variables and an array to store their values
    double **bmass, **bmass_crr, **carb, **carb_crr, flx[NUM_NNS], **tempb, **tempc, ti, tot_flx, **phi, **phi_crr, **tempphi;
    /*Pointers to biomass arrays (current and new), pointers to carbon arrays (current and new),
     * flux array for each cell (4 entries in 2d case and 2 entries in 1d case),
     * temporary pointer for biomass address, temporary pointer for carbon address,
     * time variable, sum of edge fluxes and the attachment function (current and new) and temporary
     * pointer to save the attachment function address.*/
    double **step_arr; //Pointer to array of time steps used for storing calculated time steps in the adaptive method
    double *rec_ti; // Pointer to data recording times
    double sptl_sum_bmass; //summation variable that holds the spatial integration of bmass dns
    double sptl_sum_carb; //summation variable that holds the spatial integration of carb concn
    double upper_step; //The upper-bound of the time step
    unsigned long rec_ti_size; // Size of the array pointed to by *rec_ti
    char str[100]; //String of 100 chars
    char src[100]; //String of 100 chars
    char sg[100]; //String of 100 chars
    struct stat st; // struct variable used for checking whether the folder where the files are saved, exists
    //*************************************************************************************
    //*************************************************************************************
    //openmp timing variables
    //double t1, t2;

    //Code run-time timing routine - start
    time_t curtime;
    time(&curtime);
    printf("Start @ = %s\n", ctime(&curtime));

    f_out5 = fopen ("Simulation_log", "w+"); //This file stores some parameters specific to a particular run
    params_r(); //Read the simulation parameters from the params file

    //Initialize the recording times here for all simulations with various refinements and seeds
    //*****************************************************************************************
    rec_ti_size = (unsigned long) round((END_T - INIT_T) / T_PRTN) + 1; // The size of the data recording times array; +1 to include stop time
    if ((rec_ti = malloc(rec_ti_size*sizeof(rec_ti))) == NULL) {
        fprintf(stderr,"malloc failed\n");
        exit(1);
    } //Initialize the array of data recording times
    for (j=0;j<rec_ti_size; j++){
        rec_ti[j] = ((double)j*T_PRTN) + INIT_T; //Multiply the recording time increment by the loop counter and add to start time
    }//Fill up data recording times array with the calculated values
    //*****************************************************************************************
    //*****************************************************************************************
    //The following for-loop starts the simulation and loops through grid sizes
    for (N=N_coarse; N<=N_fine; N*=2) {

        if (N % 2 != 0) {
            fprintf(stderr, "Make sure the grid size is an even number!\n");
            exit(1);
        }
        if (N < NTHRDS) {
            fprintf(stderr, "The number of allocated threads must be less than the grid size!\n");
            exit(1);
        }

        //*************************************************************************************
        //*************************************************************************************
        //MEMORY ALLOCATIONS AND FILE OPENS
        //*************************************************************************************
        //*************************************************************************************
        if ((bmass = malloc(N*sizeof(*bmass))) == NULL ||
            (bmass_crr = malloc(N*sizeof(*bmass))) == NULL) {
            fprintf(stderr, "malloc failed\n");
            exit(1);
        } //Initialize the new and current biomass arrays
        for(j=0; j<N; j++) {
            if ((bmass[j] = malloc(PADS*sizeof(*(bmass[j])))) == NULL ||
                (bmass_crr[j] = malloc(PADS*sizeof(*(bmass_crr[j])))) == NULL) {
                fprintf(stderr, "malloc failed\n");
                exit(1);
            }
        }//Pad the biomass arrays

        if ((carb =  malloc(N*sizeof(*carb))) == NULL ||
            (carb_crr = malloc(N*sizeof(*carb_crr))) == NULL) {
            fprintf(stderr, "malloc failed\n");
            exit(1);
        } //Initialize the new and current carbon arrays
        for(j=0; j<N; j++) {
            if ((carb[j] = malloc(PADS*sizeof(*(carb[j])))) == NULL ||
                (carb_crr[j] = malloc(PADS*sizeof(*(carb_crr[j])))) == NULL) {
                fprintf(stderr, "malloc failed\n");
                exit(1);
            }
        }//Pad the carbon arrays

        if ((arr_g = malloc(N*sizeof(*arr_g))) == NULL) {
            fprintf(stderr, "malloc failed\n");
            exit(1);
        } //Initialize the array of Gaussian random numbers
        for(j=0; j<N; j++) {
            if ((arr_g[j] = malloc(PADS*sizeof(*(arr_g[j])))) == NULL) {
                fprintf(stderr, "malloc failed\n");
                exit(1);
            }
        }//Pad the Gaussian random number array

        if ((phi = malloc(N*sizeof(*phi))) == NULL ||
            (phi_crr = malloc(N*sizeof(*phi_crr))) == NULL) {
            fprintf(stderr, "malloc failed\n");
            exit(1);
        } //Initialize the array of attachment function
        for(j=0; j<N; j++) {
            if ((phi[j] = malloc(PADS*sizeof(*(phi[j])))) == NULL ||
                (phi_crr[j] = malloc(PADS*sizeof(*(phi_crr[j])))) == NULL) {
                fprintf(stderr, "malloc failed\n");
                exit(1);
            }
        }//Pad the attachment function array

        if ((random_1 = malloc((N_coarse/2)*sizeof(*random_1))) == NULL ||
            (random_2 = malloc((N_coarse/2)*sizeof(*random_2))) == NULL) {
            fprintf(stderr,"malloc failed\n");
            exit(1);
        }//Initialize the array of random numbers used in the Marsaglia_polar_generator
        for(j=0; j<(N_coarse/2); j++) {
            if ((random_1[j] = malloc(PADS*sizeof(*(random_1[j])))) == NULL ||
                (random_2[j] = malloc(PADS*sizeof(*(random_2[j])))) == NULL) {
                fprintf(stderr, "malloc failed\n");
                exit(1);
            }
        }//Pad the random number arrays generated by lcg 2 and lcg 3

        if ((step_arr = malloc(N*sizeof(*step_arr))) == NULL) {
            fprintf(stderr, "malloc failed\n");
            exit(1);
        }//Initialize the array of adaptively calculated time steps
        for(j=0; j<N; j++) {
            if ((step_arr[j] = malloc(PADS*sizeof(*(step_arr[j])))) == NULL) {
                fprintf(stderr, "malloc failed\n");
                exit(1);
            }
        }//Pad the array of time steps

        //This for-loop loops through starting seeds determined by supplied parameters
        for (v=SEED_init; v<=SEED_fin; v++) {

            sprintf(src, "seed=%llu/", v);//print the given string to src
            if (stat(src, &st) == -1) {
                mkdir(src, 0700);
            }//Check whether the directory with the name stored in src exits and if not create and save it

            /*////////////////////////////////////*/
            //File creation and labeling opertions
            sprintf(str, "Biomass_density_Nx1=%lu", N); //Create the file with the given name and corresponding grid size
            strcpy(sg, src); //Copy the name of the directory stored in src to sg
            strcat(sg, str); // String concatenation i.e. append str to sg, str contains the name of the file
            f_out1 = fopen (sg, "w+"); //Open file
            sprintf(str, "Carbon_concentration_Nx1=%lu", N);
            strcpy(sg, src);
            strcat(sg, str);
            f_out2 = fopen (sg, "w+");
            sprintf(str, "Spatially_integrated_carbconcn_Nx1=%lu", N);
            strcpy(sg, src);
            strcat(sg, str);
            f_out3 = fopen (sg, "w+");
            sprintf(str, "Spatially_integrated_bmassdns_Nx1=%lu", N);
            strcpy(sg, src);
            strcat(sg, str);
            f_out4 = fopen (sg, "w+");
            sprintf(str, "phi_Nx1=%lu", N);
            strcpy(sg, src);
            strcat(sg, str);
            f_out6 = fopen (sg, "w+");

            //*************************************************************************************
            //*************************************************************************************
            //INITIAL DATA BLOCK
            //*************************************************************************************
            //*************************************************************************************
            ti = INIT_T; //Initialize time variable to start time

            nthrds = 0; //Initialize the local var nthrds to zero for each starting seed,
                        //later this is used to check whether the allocated number of threads
                        //matches that of requested each time the parallel region is constructed


            for (j=0; j<N; j++) {
                carb_crr[j][0] = INIT_C;
                carb[j][0] = INIT_C;
            } //Fill up both current and new carbon arrays to the initial carbon value

            for (j=0; j<N; j++) {
                phi_crr[j][0] = INIT_PHI;
                phi[j][0] = INIT_PHI;
            } //Fill up both current and new attachment function arrays to the provided initial value

            for (j=0; j<N; j++) {
                bmass_crr[j][0] = 0.0;
                bmass[j][0] = 0.0;
            } //The biomass arrays are initialized to zero
            bmass_crr[BMASS_INIT_LOC][0] = INIT_M; //Both the current and new biomass arrays are initialized with arbitrary choice of cells
            bmass[BMASS_INIT_LOC][0] = INIT_M;

            //Determine the upper limit of the time step size using the CFL criterion (LOWER_BMASS is acquired through experimentation for stability)
            upper_step = pow(DOMAIN_L / (double) N_coarse, 2) / D_M(LOWER_BMASS); //Calculate the upper limit of time step using the coarsest grid size
            for (j=0; j<N; j++) {
                step_arr[j][0] = upper_step;
            }//Fill up the time step array with the upper limit value calculated above
            t_num_invs = 0; //Initialize the number of time intervals to zero for each starting seed
            step = upper_step; //Initialize step var to the upper limit as well
            sqrtstep = sqrt(step); //Take the square root of step and store it
            dlt_x = DOMAIN_L / (double) N; //Calculate the cell size for the current grid size
            dlt_x_2 = pow(DOMAIN_L / (double) N, 2); //Calculate the square of cell size for the current grid size
            inv_dlt_x = (double) N / DOMAIN_L ; //Determine the inverse of the cell size for the current grid size

            /* ******************************************************************************** */
            //**************************************************************************************
            //**************************************************************************************
            //**************************************************************************************
            //**************************************************************************************

            time(&curtime);
            fprintf(f_out5, "*********start for seed=%llu , grid size=%lu*********\n", v, N);
            fprintf(f_out5, "Start time for this run: %s", ctime(&curtime));
            fprintf(f_out5, "Cell size: %f\n", dlt_x);
            fprintf(f_out5, "Upper limit of time step: %f\n", upper_step);
            fflush(f_out5); //Flush the memory to log file

            //initialize and prepare the prngs that are used for random number generation in the parallel region by leap-frogging
            parseed(v, v, random_1, random_2, &mult2_n, &mult3_n, &inc2_n, &inc3_n);

            //This for-loop loops through the rec_ti array, picks a new recoding time
            //and performs the runtime saving of the data at the selected time.
            for (q=0; q<rec_ti_size; q++) {
                /*Start a while loop and increase the time variable using the time step value
                 * and calculate biomass and carbon for the new time. Do so up until the stop
                 * time is achieved.*/
                //t1 = omp_get_wtime(); //openmp routine to save time at the beginning of the while loop

                sptl_sum_bmass = 0.0; //Initialize the spatial integration value of bmass dns to 0.0 for this iteration
                sptl_sum_carb = 0.0; //Initialize the spatial integration value of carb concn to 0.0 for this iteration

                while (ti < rec_ti[q]) {
                    tempb = bmass;    //The next five lines swap the addresses of the current biomass(carbon and phi) values and the new-
                    bmass = bmass_crr;//-ones and fill up the space allocated for the current values with the newly calculated values
                    bmass_crr = tempb;
                    tempc = carb;
                    carb = carb_crr;
                    carb_crr = tempc;
                    tempphi = phi;
                    phi = phi_crr;
                    phi_crr = tempphi;

                    t_num_invs++; //Increment the number of time steps

                    ti += step; //Increment the time step
//Entering the parallel region
                    omp_set_num_threads((int) NTHRDS); //Allocate the requested number of threads
#pragma omp parallel default(shared) private(i, j, l, id, flx, tot_flx, g1, g2)
                    {
                        nthrds = (unsigned long) omp_get_num_threads(); //Assign the number of requested threads to nthrds

#pragma omp single
                        {
                            if(nthrds != (int) NTHRDS) {
                                fprintf(stderr, "Parallel region failed!\n");
                                exit(1);
                            }
                        }//If the allocated number of threads does not match the requested one exit and throw error

                        //parallel calculation of the normal prngs using the Marsaglia Polar generator
#pragma omp for schedule(auto)
                        for (i=0; i<(N_coarse/2); i++) { //loop through half the coarsest grid size since MP_PRNG generates two prns at a time
                            MP_PRNG(&g1, &g2, random_1[i][0], random_2[i][0]);
                            for (j=0; j<(N/N_coarse); j++) {
                                arr_g[2*(N/N_coarse)*i+j][0] = g1;
                                arr_g[2*(N/N_coarse)*i+j+(N/N_coarse)][0] = g2;
                            }//For higher refinement the Gaussian prns are repeated in adjacent cells to match their spatial distribution in the coarsest grid
                        }

                        id = (unsigned long) omp_get_thread_num(); //Initialize id to the thread numbers

                        //Updating prngs used for normal prngs by the leap-frogging method
                        for (l=id;l<(N_coarse/2);l+=nthrds) {
                            random_1[l][0] = (mult2_n * random_1[l][0] + inc2_n) % MODUL;
                            random_2[l][0] = (mult3_n * random_2[l][0] + inc3_n) % MODUL;
                        }

                        //Parallel for-loop for calculating the biomass and carbon using the Euler method
#pragma omp for schedule(auto)
                        for (i=0; i<N; i++) {
                            if (i == 0) { //If in leftmost cell, there is no left flux
                                flx[LT] = 0.0;
                                flx[RT] = flux(bmass_crr[i][0], bmass_crr[i+1][0]);
                                tot_flx = flx[LT] + flx[RT];
                                phi[i][0] = Euler_phi(phi_crr[i][0], bmass_crr[i][0], carb_crr[i][0], arr_g[i][0]);//Update phi using Euler-Maruyama scheme
                                bmass[i][0] = Euler_biomass(bmass_crr[i][0], tot_flx, carb_crr[i][0], phi_crr[i][0]);//Updage bmass using Euler scheme
                                carb[i][0] = Euler_carb(carb_crr[i][0], bmass_crr[i][0]);//Update carb using Euler scheme
                            }
                            else if (i == (N-1)) { //If in rightmost cell, there is no right flux
                                flx[LT] = flux(bmass_crr[i][0], bmass_crr[i-1][0]);
                                flx[RT] = 0.0;
                                tot_flx = flx[LT] + flx[RT];
                                phi[i][0] = Euler_phi(phi_crr[i][0], bmass_crr[i][0], carb_crr[i][0], arr_g[i][0]);
                                bmass[i][0] = Euler_biomass(bmass_crr[i][0], tot_flx, carb_crr[i][0], phi_crr[i][0]);
                                carb[i][0] = Euler_carb(carb_crr[i][0], bmass_crr[i][0]);
                            }
                            else { //All other cells have right and left fluxes
                                flx[LT] = flux(bmass_crr[i][0], bmass_crr[i-1][0]);
                                flx[RT] = flux(bmass_crr[i][0], bmass_crr[i+1][0]);
                                tot_flx = flx[LT] + flx[RT];
                                phi[i][0] = Euler_phi(phi_crr[i][0], bmass_crr[i][0], carb_crr[i][0], arr_g[i][0]);
                                bmass[i][0] = Euler_biomass(bmass_crr[i][0], tot_flx, carb_crr[i][0], phi_crr[i][0]);
                                carb[i][0] = Euler_carb(carb_crr[i][0], bmass_crr[i][0]);
                            }
                        }

#pragma omp single
                        {
                            step = upper_step; //Set time step to upper value in order to find its minimum possible value below (adaptive approach)
                        }

#pragma omp for schedule(auto) reduction(min:step)
                        for (i=0; i<N; i++) { //Go through all cells and apply CFL formula
                            if (bmass[i][0] > 0.0 && bmass[i][0] < UPPER_BMASS) {
                                step_arr[i][0] = dlt_x_2 / (D_M(bmass[i][0]));
                                step = min(step,step_arr[i][0]); //Find the minimum value of the time step using the CFL condition
                            }
                        }

#pragma omp single
                        {
                            sqrtstep = sqrt(step);//Update the square root of the time step using the adaptively obtained time step
                        }
                    }
                }

                //t2 = omp_get_wtime(); //openmp routine to save time at the end of the while loop

                //*******************************************************************************************
                //PRINTING STOP-TIMES' DATA TO FILES
                //*******************************************************************************************
                fprintf(f_out1, "%f ", ti); //Save the current time and the biomass dns values and the carbon conn. values-
                fprintf(f_out2, "%f ", ti); //-at this time to the output files.
                fprintf(f_out6, "%f ", ti); // File for saving the phi values
                for (j=0; j<N; j++) {
                    fprintf(f_out1, "%f ", bmass[j][0]);
                    fprintf(f_out2, "%f ", carb[j][0]);
                    fprintf(f_out6, "%f ", phi[j][0]);
                    sptl_sum_bmass += bmass[j][0]*dlt_x; // Multiply the bmass dns at this position by the cell size
                    sptl_sum_carb += carb[j][0]*dlt_x; // Multiply the carb concn at this position by the cell size
                }
                fprintf(f_out4, "%f %15.14f\n", ti, sptl_sum_bmass); //Save the spatial integration of biomass at this time to file
                fprintf(f_out3, "%f %15.14f\n", ti, sptl_sum_carb); //Save the spatial integration of carbon at this time to file
                fprintf(f_out1, "\n");
                fprintf(f_out2, "\n");
                fprintf(f_out6, "\n");
                //*******************************************************************************************
                //*******************************************************************************************
            }
            fprintf(f_out5, "# of time steps: %ld\n", t_num_invs); //Save simulation outputs to log file
            time(&curtime);
            fprintf(f_out5, "The number of allocated threads for this run: %lu\n", nthrds);
            fprintf(f_out5, "Stop time for this run: %s", ctime(&curtime));
            fprintf(f_out5, "*********end for seed=%llu , grid size=%lu*********\n\n", v, N);
            fflush(f_out5); //Flush memory to log file

            fclose(f_out1);
            fclose(f_out2);
            fclose(f_out3);
            fclose(f_out4);
            fclose(f_out6);
        }

        for(j=0; j<N; j++) { //Freeing space in memory at the end of current grid size simulation
            free(bmass_crr[j]);
            free(carb_crr[j]);
            free(bmass[j]);
            free(carb[j]);
            free(phi[j]);
            free(phi_crr[j]);
            free(arr_g[j]);
            free(step_arr[j]);
        }
        for(j=0; j<(N_coarse/2); j++) {
            free(random_1[j]);
            free(random_2[j]);
        }
        free(random_1);
        free(random_2);
        free(bmass_crr);
        free(carb_crr);
        free(bmass);
        free(carb);
        free(phi);
        free(phi_crr);
        free(arr_g);
        free(step_arr);
    }
    free(rec_ti); //Free memory for recording times at the end of simulations for all grid sizes and seeds
    fclose(f_out5); //Simulation log file closed after all is done

    //Code run-time timing routine - stop
    time(&curtime);
    printf("Stop @ = %s", ctime(&curtime));
    //printf("The time calculated by openmp time routine: %.8f\n", t2 - t1);

    return 0;
}
