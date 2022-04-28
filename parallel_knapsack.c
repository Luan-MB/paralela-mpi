#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <sys/time.h>

#define ROOT 0
#define STD_TAG 0

double timestamp() {
    struct timeval tp;
    gettimeofday(&tp, NULL);
    return((double)(tp.tv_sec*1000.0 + tp.tv_usec/1000.0));
}

// A utility function that returns
// maximum of two integers
int max(int a, int b) {return (a > b) ? a : b;}

// Driver program to test above function
int main(int argc, char **argv)
{
	int num_procs;
    int *upper_row, *lower_row;
    MPI_Status status;
    int rank;
    FILE *arq = fopen(argv[1], "r");
    double parallel_time, sequential_time;

    // initialize MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int sub_size, division_rest;
    int max_weight;
    int n, W;

    int *val, *wt;

    if (rank == ROOT) {
        sequential_time = timestamp();
        fscanf(arq, "%d %d", &n, &W);
    }

    MPI_Bcast(&n, 1, MPI_INT, ROOT, MPI_COMM_WORLD);

    val = (int*) malloc(n * sizeof(int));
	wt = (int*) malloc(n * sizeof(int));

    if (rank == ROOT) {

        max_weight = 0;

        // read all profits and weights
        for (int i = 0; i < n; ++i) {
            fscanf(arq, "%d %d", &(val[i]), &(wt[i]));
            // find the heaviest item weight
            if (max_weight < wt[i])
                max_weight = wt[i];
        }

        upper_row = calloc (W + 1, sizeof(int));
        lower_row = calloc (W + 1, sizeof(int));

        sub_size = (W + 1) / num_procs;
        division_rest = (W + 1) % num_procs;
    }

    // broadcast important variables
    MPI_Bcast(&sub_size, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(val, n, MPI_INT, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(wt, n, MPI_INT, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&max_weight, 1, MPI_INT, ROOT, MPI_COMM_WORLD);

    int *aux_row;

    if (rank != ROOT)
        aux_row = malloc(sub_size * rank * sizeof(int));

    int counts[num_procs];
    int offsets[num_procs];

    if (rank == ROOT) {

        // calculate the last process array size (can be the same as sub_size)
        int last_proc_size = sub_size + division_rest;
        // send the size to last process
        MPI_Send(&last_proc_size, 1, MPI_INT, num_procs - 1, STD_TAG, MPI_COMM_WORLD);

        // set the counts and offsets for each process
        for (int i = 0; i < num_procs - 1; ++i) {
            counts[i] = sub_size;
            offsets[i] = i * sub_size;
        }

        counts[num_procs - 1] = last_proc_size;
        offsets[num_procs - 1] = (num_procs - 1) * sub_size;
    }

    if (rank == num_procs - 1) {
        // the last process receives its size
        MPI_Recv(&sub_size, 1, MPI_INT, ROOT, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    }

    // allocate the sub arrays for each process
    int *sub_upper_row = malloc(sub_size * sizeof(int));
    int *sub_lower_row = malloc(sub_size * sizeof(int));

    // scatter the original arrays beetween all processes
    MPI_Scatterv(upper_row, counts, offsets, MPI_INT, sub_upper_row, sub_size, MPI_INT, ROOT, MPI_COMM_WORLD);
    MPI_Scatterv(lower_row, counts, offsets, MPI_INT, sub_lower_row, sub_size, MPI_INT, ROOT, MPI_COMM_WORLD);

    if (rank == ROOT) {
        parallel_time = timestamp();
    }

    for(int i = 0; i < n; i++) {
        // if line i is even
        if ((i % 2) == 0) {
            // send last calculated array to next processes
            if (rank != num_procs - 1) {
                for (int k = (rank + 1); k < num_procs; ++k) {
                    MPI_Send(sub_upper_row, sub_size, MPI_INT, k, STD_TAG, MPI_COMM_WORLD);
                }
            }
            // receive last calculated array from previous processes
            if (rank != ROOT) {
                for (int k = 0; k < rank; ++k) {
                    MPI_Recv(&aux_row[k * sub_size], sub_size, MPI_INT, k, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                }
            }

            for(int j = 0; j < sub_size; j++) {
                if(wt[i] <= (j + rank * sub_size)) { // could put item in knapsack
                    int previous_value = sub_upper_row[j];
                    int replace_items;
                    int index = j - wt[i];
                    // if index is inside the process sub array
                    if (index >= 0)
                        replace_items = val[i] + sub_upper_row[index];
                    // if index is out of process sub array bounds
                    // get value from aux array
                    else
                        replace_items = val[i] + aux_row[rank * sub_size + index];

                    // is it better to keep what we already got,
                    // or is it better to swap whatever we have in the bag that weights up to `j`
                    // and put item `i`?
                    sub_lower_row[j] = max(previous_value, replace_items);
                }
                else {
                    // can't put item `i`
                    sub_lower_row[j] = sub_upper_row[j];
                }
            }
        // if line i is odd
        } else {
            // send a part of the last calculated array to next process
            if (rank != num_procs - 1) {
                for (int k = (rank + 1); k < num_procs; ++k) {
                    MPI_Send(sub_lower_row, sub_size, MPI_INT, k, STD_TAG, MPI_COMM_WORLD);
                }
            }
            // receive a part of the last calculated array from the previous process
            if (rank != ROOT) {
                for (int k = 0; k < rank; ++k) {
                    MPI_Recv(&aux_row[k * sub_size], sub_size, MPI_INT, k, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                }
            }

            for(int j = 0; j < sub_size; j++) {
                if(wt[i] <= (j + rank * sub_size)) { // could put item in knapsack
                    int previous_value = sub_lower_row[j];
                    int replace_items;
                    int index = j - wt[i];
                    // if index is inside the process sub array
                    if (index >= 0)
                        replace_items = val[i] + sub_lower_row[index];
                    // if index is out of process sub array bounds
                    // get value from aux array
                    else
                        replace_items = val[i] + aux_row[rank * sub_size + index];

                    // is it better to keep what we already got,
                    // or is it better to swap whatever we have in the bag that weights up to `j`
                    // and put item `i`?
                    sub_upper_row[j] = max(previous_value, replace_items);
                }
                else {
                    // can't put item `i`
                    sub_upper_row[j] = sub_lower_row[j];
                }
            }
        }

        // wait until all processes reach barrier
        // to make sure all arrays were correctly calculated
        MPI_Barrier(MPI_COMM_WORLD);
    }

    if (rank == ROOT)
        parallel_time = timestamp() - parallel_time;

    // gather all sub arrays in the orignal array
    MPI_Gatherv(sub_upper_row, sub_size, MPI_INT, upper_row, counts, offsets, MPI_INT, ROOT, MPI_COMM_WORLD);
    MPI_Gatherv(sub_lower_row, sub_size, MPI_INT, lower_row, counts, offsets, MPI_INT, ROOT, MPI_COMM_WORLD);

    if (rank == ROOT)
        sequential_time = timestamp() - sequential_time;

    if (rank == ROOT) {

        int retval;

        if ((n % 2) == 0)
            retval = upper_row[W]; 
        else
            retval = lower_row[W];

        printf("%d     %g     %g\n", retval, parallel_time, sequential_time);

        free(lower_row);
        free(upper_row);
    }

        free(val);
        free(wt);

    MPI_Finalize();
    fclose(arq);

    return 0;
} 