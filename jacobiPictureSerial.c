/*
John Burke ~ Comp 233 A ~ Jacobi With Picture Serial

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

The convergence value can be printed at the end of every 1000 iterations

After the final iteration, the values in the array are printed as a pixel map
The red is hot (100 degrees is the highest) and the blue is cold (0 degrees is the lowest)

The time it takes to run is calculated and is printed to a txt file

The epsilon value and the max iterations are both received from the command line.

base code taken from 
https://www.mcs.anl.gov/research/projects/mpi/tutorial/mpiexmpl/src/jacobi/C/main.html

"timer.h" taken from 
Pacheco, Peter. An Introduction to Parallel Programming. Elsevier Science &amp; Technology, 2018. 

must be compiled with "-lm" added to the end for the sqrt() function
compile: gcc jacobiPictureSerial.c -o jacobiPictureSerial -lm
run: ./jacobiPictureSerial <epsilon> <max iterations>
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "timer.h"

//the MAX_R and MAX_C can be changed to give different dimensions to the array

//the dimensions of the grid
#define MAX_R 400    //num rows
#define MAX_C 800    //num cols

//if defined, then print the iterations
//comment out if you don't want to print the iterations
//#define PRINT_ITER 


int main( argc, argv )
int argc;
char **argv;
{
    /*
    start_time is where the timing begins
    stop_time is where the timing ends
    seconds_elapsed is the difference between start and stop_time, which gives the amount of time
        it took to run the program

    iteration_count counts the number of iterations that have been done over the array

    r_first is the first row of the array that is not a boundry
    r_last is the last row of the array that is not a boundry

    r and c are counters for rows and columns

    red, green, and blue are used for the colors printed to the ppm file
    red is the (temperature/100) * 255, green is always 0, and blue is calculated from the red value

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
    */
    double seconds_elapsed, start_time, stop_time;
    long int iteration_count;
    int r_first = 1, r_last = MAX_R - 1, r, c;
    int red, green = 0, blue;
    int north = 100, east = 100, south = 0, west = 0;
    double interior = -1.0;
    double diff_norm, converge_point;
    FILE *fptr;
    

    //get the epsilon and max_interations from the command line
    char *endPtr;   //used to set the end of the string
    double epsilon = strtod(*(argv + 1), &endPtr);
    long int max_iterations = strtol(*(argv + 2), &endPtr, 10);

    
    //prints greeting and description
    printf("John Burke ~ Comp 233 A ~ Jacobi With Picture Serial\n");
    printf("Runs an approximation of the Jacobi method with %d rows and %d cols.\n", MAX_R, MAX_C);
    printf("The boundry and interior values are set as follows:\n");
    printf("North: %d\tSouth: %d\tEast: %d\tWest: %d", north, south, east, west);
    printf("\t\tInterior: %.2lf\n", interior);
    printf("The given epsilon is %lf, and the given max iterations is %li\n\n", epsilon, max_iterations);

    //start timer
    GET_TIME(start_time);
    

    //create the arrays

    //the values of the current array
    double **x_current = (double **)malloc(MAX_R * sizeof(double *));
    for(r = 0; r < MAX_R; r++){
        x_current[r] = (double *)malloc(MAX_C * sizeof(double));
    }

    //the new values after the jacobi calculation
    double **x_new = (double **)malloc(MAX_R * sizeof(double *));
    for(r = 0; r < MAX_R; r++){
        x_new[r] = (double *)malloc(MAX_C * sizeof(double));
    }

    //used to swap the arrays
    double **temp;


    //set the interior, east, and west values
    //from r_first to r_last because we don't want to do the top & bottom rows
    for (r = r_first; r < r_last; r++) {   
        for (c = 0; c < MAX_C; c++) {
            //if on the west (left side), set to west value
            if(c == 0){
                x_current[r][c] = west;
                x_new[r][c] = west;
            //if on the east (right side), set to east value
            }else if(c == MAX_C  -1){
                x_current[r][c] = east;
                x_new[r][c] = east;
            //if not east or west side (so in the middle), then set to interior value
            }else {
                x_current[r][c] = interior;
                x_new[r][c] = interior;
            }
        }
    }

    //set the north and south boundries
    for (c = 0; c < MAX_C; c++) {
        x_current[r_first - 1][c] = north;
        x_current[r_last][c] = south;

        x_new[r_first - 1][c] = north;
        x_new[r_last][c] = south;
    }

    //start at 0 iterations
    iteration_count = 0;

    /*
    will keep doing iterations until it reaches the max number of iterations or if the difference
    value converges to be less than epsilon
    */
    do {
        /* Compute new values (but not on boundary) */
        iteration_count++;
        diff_norm = 0.0;
        for (r = r_first; r < r_last; r++) {
            for (c = 1; c < MAX_C - 1; c++) {
                //calculate the new value for each cell
                x_new[r][c] = (x_current[r][c+1] + x_current[r][c-1] + x_current[r+1][c] + x_current[r-1][c]) / 4.0;
                /*
                find the difference between the new value and the old one
                square it (to remove negatives) and add to the total
                the square root will be taken later witch gives up the absolute value
                */
                diff_norm += (x_new[r][c] - x_current[r][c]) * (x_new[r][c] - x_current[r][c]);
            }
        }

        //swap the arrays so the x_current has the new values
        temp = x_new;
        x_new = x_current;
        x_current = temp;
        temp = NULL;
        

        //square root because diff_norm was squared
        converge_point = sqrt( diff_norm );

        //prints the iterations if PRINT_ITER is defined
        #ifdef PRINT_ITER
        //prints the results for every 1000 iterations 
        if(iteration_count % 1000 == 0){
            printf( "At iteration %li, diff is %e\n", iteration_count, converge_point );
        }
        #endif
    
    } while (converge_point > epsilon && iteration_count < max_iterations);

    //print results
    printf("Converged at %lf in %li iterations.\n", converge_point, iteration_count);

    //stop timer
    GET_TIME(stop_time);
    

    //the final version of the array is printed to a ppm file

    fptr = fopen("./jacobiPicture.ppm", "wb");

    //print P6 to indicate binary version
    fprintf(fptr, "P6\n");

    //print header 
    fprintf(fptr, "#John Burke ~ Comp 233 A ~ Jacobi With Picture\n"
        "#Jacobi final state on %d x %d plate; epsilon 1e-2 was used\n", MAX_R, MAX_C);

    //print size and color value
    fprintf(fptr, "%d %d\n255\n", MAX_C, MAX_R);


    //print RGB values for each cell in the array
    for(r = 0; r < MAX_R; r++){
        for(c = 0; c < MAX_C; c++){
            //red = (temperature/100) * 255
            red = (x_current[r][c] / 100) * 255;
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
    }
    
    //close file
    fclose(fptr);
    fptr = NULL;
    

    //free arrays

    //free x_current
    for(r = 0; r < MAX_R; r++){
        free(x_current[r]);
    }
    free(x_current);
    x_current = NULL;

    //free x_new
    for(r = 0; r < MAX_R; r++){
        free(x_new[r]);
    }
    free(x_new);
    x_new = NULL;

    //free temp
    free(temp);
    temp = NULL;


    //calculate time, print to file, print to terminal
    seconds_elapsed = stop_time - start_time;

    //open file
    fptr = fopen("./jacobiTimes.txt", "w");

    //print time to file
    fprintf(fptr, "Serial Time: %.3lf\n\n", seconds_elapsed);

    //close file
    fclose(fptr);
    fptr = NULL;
    
    //print time to terminal
    printf("Finished in %.3lf seconds.", seconds_elapsed);


    //done
    printf("\n<< Normal Termination >>\n");
    return 0;
}