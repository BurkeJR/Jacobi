/*
John Burke ~ Comp 233 A ~ Jacobi With Picture MPI

Runs an approximation of the Jacobi method

Starts with an array with a pre-set number of rows and cols
Sets the boundries to pre-set numbers
    - the boundries are the perimeter of the array
The interior cells are all given a pre-set value

For each iteration, the interior cells are set equal to the average of the 4 cells it touches
    - above, below, to the right and left
    - the boundry cells are not changed

The array continues to go through iterations until:
    1. the difference between the values of the cells in consecutive iterations converges to
        a value less than the given epsilon
    2. the number of iterations reaches the given maximum number of iterations

The convergence value can be printed at the every 1000 iterations to the terminal

After the final iteration, the values in the array are printed as a pixel map
The red is hot (100 degrees is the highest) and the blue is cold (0 degrees is the lowest)

The final version of the array can be printed to the terminal

The time it takes to run is calculated and is printed to a txt file

The epsilon value and the max iterations are both received from the command line.

code taken from https://www.mcs.anl.gov/research/projects/mpi/tutorial/mpiexmpl/src/jacobi/C/main.html

must be compiled with "-lm" added to the end for the sqrt() function
compile: mpicc jacobi.c -o jacobi -lm
run: mpirun --oversubscribe -np 4 ./jacobi <epsilon> <max iterations>
*/

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "mpi.h"
#include <time.h>

//the MAX_R and MAX_C can be changed to give different dimensions to the array

//the dimensions of the grid
#define MAX_R 400    //num rows
#define MAX_C 800    //num cols

//if defined, then print the iterations
//comment out if you don't want to print the iterations
//#define PRINT_ITER 

//if defined, then print the final version of the array
//#define PRINT_FINAL_ARRAY


int main( argc, argv )
int argc;
char **argv;
{
    /*
    start_time is where the timing begins
    stop_time is where the timing ends
    seconds_elapsed is the difference between start and stop_time, which gives the amount of time
        it took to run the program
    only master tracks the time

    my_rank is the rank of the process running the program
    comm_size is the number of processes in the COMM_WORLD

    iteration_count counts the number of iterations that have been done over the array

    r_first is where the processes first row is in the array
    r_last is where the processes last row is in the array
    these are different depending on the process due to boundry rows

    r and c are counters for rows and columns
    source is a counter for receiving information from other processes

    red, green, and blue are used for the colors printed to the ppm file
    red is (temperature/100) * 255, green is always 0, and blue is calculated from the red value

    north, south, east, west are pre-set values for the boundries of the array

    interior is the pre-set value for the interior cells of the array

    status is used for the status parameter in MPI_Receive

    diff_norm is used to store the sum of the differences between the new values and the old values of
        the cells in the array 

    converge_point is the value that all of the interior cells are converging too
        it is compared with epsilon to determine when to stop iterating

    epsilon is the value where the iterations stop if the cells converge to a value less than it

    max_iterations is the maximum number of iterations before the program stops

    fptr is used to print to the timing file and the ppm file
    only master prints to the file
    */

    double seconds_elapsed, start_time, stop_time;
    int my_rank, comm_size;
    long int iteration_count;
    int r_first, r_last, r, c, source;
    int red, green = 0, blue;
    int north = 100, east = 100, south = 0, west = 0;
    double interior = -1.0;
    MPI_Status status;
    double diff_norm, converge_point;
    FILE *fptr;

    //get the epsilon and max_interations from the command line
    char *endPtr;   //used to set the end of the string
    double epsilon = strtod(*(argv + 1), &endPtr);
    long int max_iterations = strtol(*(argv + 2), &endPtr, 10);

    //set up MPI
    MPI_Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );
    MPI_Comm_size( MPI_COMM_WORLD, &comm_size );

    //the master prints greeting and description and starts the timer
    if(my_rank == 0){
        printf("John Burke ~ Comp 233 A ~ Jacobi With Picture MPI\n");
        printf("Runs an approximation of the Jacobi method on %d processes.\n", comm_size);
        printf("Size is %d rows and %d cols.", MAX_R, MAX_C);
        printf("The boundry and interior values are set as follows:\n");
        printf("North: %d\tSouth: %d\tEast: %d\tWest: %d", north, south, east, west);
        printf("\t\tInterior: %.2lf\n", interior);
        printf("The given epsilon is %lf, and the given max iterations is %li\n\n", epsilon, max_iterations);

        //start timer
        start_time = clock();
    }


    //creates an array for each process to have a chunk of the data
    //2 is being added for the 2 ghost rows

    double **x_local = (double **)malloc((MAX_R/comm_size + 2) * sizeof(double *));
    for(r = 0; r < MAX_R/comm_size + 2; r++){
        x_local[r] = (double *)malloc(MAX_C * sizeof(double));
    }

    double **x_new = (double **)malloc((MAX_R/comm_size + 2) * sizeof(double *));
    for(r = 0; r < MAX_R/comm_size + 2; r++){
        x_new[r] = (double *)malloc(MAX_C * sizeof(double));
    }

    //used to swap the arrays
    double **temp;


    /* xlocal[][0] is lower ghostpoints, xlocal[][MAX_N+2] is upper */

    /* Note that top and bottom processes have one less row of interior
       points */
    r_first = 1;
    r_last  = MAX_R/comm_size;

    //increase r_first for master to make up for the one less row
    if (my_rank == 0) {
        r_first++;
    }

    //decrease r_last for last process to make up for the one less row
    if (my_rank == comm_size - 1) {
        r_last--;
    }

    //set the interior, east, and west values
    for (r = 1; r <= MAX_R/comm_size; r++) {
        for (c = 0; c < MAX_C; c++) {
            //if on the west (left side), set to west value
            if(c == 0){
                x_local[r][c] = west;
                x_new[r][c] = west;
            //if on the east (right side), set to east value
            }else if(c == MAX_C  -1){
                x_local[r][c] = east;
                x_new[r][c] = east;
            //if not east or west side (so in the middle), then set to interior value
            }else {
                x_local[r][c] = interior;
                x_new[r][c] = interior;
            }
        }
    }

    //set the north and south boundries
    for (c = 0; c < MAX_C; c++) {
        x_local[r_first-1][c] = north;
        x_local[r_last+1][c] = south;
    
        x_new[r_first-1][c] = north;
        x_new[r_last+1][c] = south;
    }


    //start at 0 iterations
    iteration_count = 0;

    /*
    will keep doing iterations until it reaches the max number of iterations or if the difference
    value converges to be less than epsilon
    */
    do {
        /* 
        Note the use of xlocal[i] for &xlocal[i][0] 
        sends only the rows that are needed by the other processes
        does not send the entire array
        */

        /* Send up unless I'm at the top, then receive from below */
        if (my_rank < comm_size - 1) {
            MPI_Send( x_local[MAX_R/comm_size], MAX_C, MPI_DOUBLE, my_rank + 1, 0, MPI_COMM_WORLD );
        }

        if (my_rank > 0){
            MPI_Recv( x_local[0], MAX_C, MPI_DOUBLE, my_rank - 1, 0, MPI_COMM_WORLD, &status );
        }

        /* Send down unless I'm at the bottom */
        if (my_rank > 0) {
            MPI_Send( x_local[1], MAX_C, MPI_DOUBLE, my_rank - 1, 1, MPI_COMM_WORLD );
        }

        if (my_rank < comm_size - 1) {
            MPI_Recv( x_local[MAX_R/comm_size+1], MAX_C, MPI_DOUBLE, my_rank + 1, 1, MPI_COMM_WORLD, &status );
        }

        /* Compute new values (but not on boundary) */
        iteration_count ++;
        diff_norm = 0.0;
        for (r = r_first; r <= r_last; r++) {
            for (c = 1; c < MAX_C - 1; c++) {
                //calculate the new value for each cell
                x_new[r][c] = (x_local[r][c+1] + x_local[r][c-1] + x_local[r+1][c] + x_local[r-1][c]) / 4.0;
                /*
                find the difference between the new value and the old one
                square it (to remove negatives) and add to the total
                the square root will be taken later witch gives up the absolute value
                */
                diff_norm += (x_new[r][c] - x_local[r][c]) * (x_new[r][c] - x_local[r][c]);
            }
        }

        //swap the arrays so the x_current has the new values
        temp = x_new;
        x_new = x_local;
        x_local = temp;
        temp = NULL;

        /*
        gets the sum of all of the squared differences between the original and the new value
        not using reduce because we need to send the result to all of the processes
        */
        MPI_Allreduce( &diff_norm, &converge_point, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

        //square root because diff_norm was squared
        converge_point = sqrt( converge_point );

        //print iterations if defined
        #ifdef PRINT_ITER
        //master prints the results for every 1000 iterations
        if (my_rank == 0) {
            if(iteration_count % 1000 == 0){
                printf( "At iteration %li, diff is %e\n", iteration_count, converge_point );
            }
        }
        #endif
    
    } while (converge_point > epsilon && iteration_count < max_iterations);

    //master prints final result
    if(my_rank == 0){
        printf("Converged at %lf in %li iterations.\n", converge_point, iteration_count);
    }


    //master stops timer
    if(my_rank == 0){
        //stop timer
        stop_time = clock();
    }
    

    /*
    master gets the final version of the array from each process
    finds red and blue values and prints to the ppm file

    if PRINT_FINAL_ARRAY is defined, master will print the final iteration of the array
    to the terminal
    */
    if(my_rank == 0){
        //temp array to hold the arrays received from other processes
        temp = (double **)malloc((MAX_R/comm_size + 2) * sizeof(double *));
        for(r = 0; r < MAX_R/comm_size + 2; r++){
            temp[r] = (double *)malloc(MAX_C * sizeof(double));
        }

        //open file
        fptr = fopen("./jacobiPicture.ppm", "wb");

        //print P6 to indicate binary version
        fprintf(fptr, "P6\n");

        //print header 
        fprintf(fptr, "#John Burke ~ Comp 233 A ~ Jacobi With Picture\n"
            "#Jacobi final state on %d x %d plate; epsilon 1e-2 was used\n", MAX_R, MAX_C);

        //print size and color value
        fprintf(fptr, "%d %d\n255\n", MAX_C, MAX_R);


        //start printing color values of the cells in the array

        //print top rows of table (the interior rows of 0 and the end boundry row)
        //is r_first - 1 in order to include the boundry
        for (r = r_first - 1; r <= r_last; r++) {
            for (c=0; c < MAX_C; c++) {
                //print to terminal if defined
                #ifdef PRINT_FINAL_ARRAY
                printf("%.1lf\t", x_local[r][c]);
                #endif

                //red = (temperature/100) * 255
                red = (x_local[r][c] / 100) * 255;
                blue = 255 - red;

                /*
                write the colors to a binary file
                static unsigned char method taken from
                https://rosettacode.org/wiki/Bitmap/Write_a_PPM_file#C
                */
                static unsigned char color[3];
                color[0] = red;
                color[1] = green;
                color[2] = blue;
                fwrite(color, sizeof(color), 1, fptr);
            }
            #ifdef PRINT_FINAL_ARRAY
            printf("\n");
            #endif
        }

        #ifdef PRINT_FINAL_ARRAY
        //extra line between process printouts for readability 
        printf("\n");
        #endif

        //print the table rows of the other processes, starting with the last process
        for(source = 1; source < comm_size; source++){
            int temp_first = 1, temp_last = MAX_R/comm_size, index;

            //receive the arrays from the other processes
            for(index = 0; index < temp_last + 2; index++){
                MPI_Recv(temp[index], MAX_C, MPI_DOUBLE, source, 0, MPI_COMM_WORLD, &status);
            }
            
            //print the rows 
            for (r = temp_first; r <= temp_last; r++) {
                for (c=0; c < MAX_C; c++) {
                    #ifdef PRINT_FINAL_ARRAY
                    printf("%.1lf\t", temp[r][c]);
                    #endif

                    //red = (temperature/100) * 255
                    red = (temp[r][c] / 100) * 255;
                    blue = 255 - red;
                    
                    /*
                    write the colors to a binary file
                    static unsigned char method taken from
                    https://rosettacode.org/wiki/Bitmap/Write_a_PPM_file#C
                    */
                    static unsigned char color[3];
                    color[0] = red;
                    color[1] = green;
                    color[2] = blue;
                    fwrite(color, sizeof(color), 1, fptr);
                }
                #ifdef PRINT_FINAL_ARRAY
                printf("\n");
                #endif
            }
            #ifdef PRINT_FINAL_ARRAY
            printf("\n");
            #endif
        }       

        //close file
        fclose(fptr);
        fptr = NULL;

    }
    else{      //if not master, then send your array to master 
        int index;        
        //send the rows 
        for(index = 0; index < MAX_R/comm_size+2; index++){
            MPI_Send(x_local[index], MAX_C, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        }
    }  

    //free the arrays

    //free x_local
    for(r = 0; r < MAX_R/comm_size + 2; r++){
        free(x_local[r]);
    }
    free(x_local);
    x_local = NULL;

    //free x_new
    for(r = 0; r < MAX_R/comm_size + 2; r++){
        free(x_new[r]);
    }
    free(x_new);
    x_new = NULL;

    //only master malloc'd temp, so only master needs to free temp's rows
    if(my_rank == 0){
        for(r = 0; r < MAX_R/comm_size + 2; r++){
            free(temp[r]);
        }
    }
    free(temp);
    temp = NULL;

    //master finds time and prints to terminal and to file
    if(my_rank == 0){

        //calculate how many seconds the program took to run
        seconds_elapsed = (stop_time - start_time) / 1000000;

        //open file
        //append because we may be adding multiple time values because we are running the program
        //with several different amounts of processes
        fptr = fopen("./jacobiTimes.txt", "a");

        //print time to file
        fprintf(fptr, "MPI %d processes time: %.3lf\n", comm_size, seconds_elapsed);

        //close file
        fclose(fptr);
        fptr = NULL;

        //print time to console
        printf("Finished in %.3lf seconds.", seconds_elapsed);
    }

    //done
    if(my_rank == 0){
        printf("\n<< Normal Termination >>\n");
    }
    MPI_Finalize( );
    return 0;
}