#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define ROOT 0
#define STD_TAG 0

// A utility function that returns
// maximum of two integers
int max(int a, int b) {return (a > b) ? a : b;}

// Driver program to test above function
int main(int argc, char **argv)
{
	int num_procs;
    int *upper_row, *lower_row;
    MPI_Status status;
    FILE *arq = fopen("input.in", "r");

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int sub_size, division_rest;
    int n, W;

    int *val, *wt;

    if (rank == ROOT) {
        fscanf(arq, "%d %d", &n, &W);
    }

    MPI_Bcast(&n, 1, MPI_INT, ROOT, MPI_COMM_WORLD);

    val = (int*) malloc(n * sizeof(int));
	wt = (int*) malloc(n * sizeof(int));

    if (rank == ROOT) {
        for (int i = 0; i < n; ++i) {
            fscanf(arq, "%d %d", &(val[i]), &(wt[i])); 
        }

        upper_row = calloc (W + 1, sizeof(int));
        lower_row = calloc (W + 1, sizeof(int));

        sub_size = (W + 1) / num_procs;
        division_rest = (W + 1) % num_procs;
    }

    MPI_Bcast(&sub_size, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(val, n, MPI_INT, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(wt, n, MPI_INT, ROOT, MPI_COMM_WORLD);

    int *aux_row;

    if (rank != ROOT)
        aux_row = malloc(sub_size * sizeof(int));

    int counts[num_procs];
    int offsets[num_procs];

    if (rank == ROOT) {

        int last_proc_size = sub_size + division_rest;
        MPI_Send(&last_proc_size, 1, MPI_INT, num_procs - 1, STD_TAG, MPI_COMM_WORLD);

        for (int i = 0; i < num_procs - 1; ++i) {
            counts[i] = sub_size;
            offsets[i] = i * sub_size;
        }

        counts[num_procs - 1] = last_proc_size;
        offsets[num_procs - 1] = (num_procs - 1) * sub_size;
    }

    if (rank == num_procs - 1) {
        MPI_Recv(&sub_size, 1, MPI_INT, ROOT, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    }

    int *sub_upper_row = malloc(sub_size * sizeof(int));
    int *sub_lower_row = malloc(sub_size * sizeof(int));

    MPI_Scatterv(upper_row, counts, offsets, MPI_INT, sub_upper_row, sub_size, MPI_INT, ROOT, MPI_COMM_WORLD);
    MPI_Scatterv(lower_row, counts, offsets, MPI_INT, sub_lower_row, sub_size, MPI_INT, ROOT, MPI_COMM_WORLD);

    for(int i = 0; i < n; i++) {

        if ((i % 2) == 0) {
            if (rank != num_procs - 1) {
                MPI_Send(&sub_upper_row[sub_size - 100], 100, MPI_INT, rank + 1, STD_TAG, MPI_COMM_WORLD);
            }
            if (rank != ROOT) {
                MPI_Recv(aux_row, 100, MPI_INT, rank - 1, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            }

            for(int j = 0; j < sub_size; j++) {
                if(wt[i] <= (j + rank * sub_size)) { // could put item in knapsack
                    int previous_value = sub_upper_row[j];
                    int replace_items;
                    int index = j - wt[i];
                    if (index >= 0)
                        replace_items = val[i] + sub_upper_row[index];
                    else
                        replace_items = val[i] + aux_row[100 + index];

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
        
        } else {
            if (rank != num_procs - 1) {
                MPI_Send(&sub_lower_row[sub_size - 100], 100, MPI_INT, rank + 1, STD_TAG, MPI_COMM_WORLD);
            }
            if (rank != ROOT) {
                MPI_Recv(aux_row, 100, MPI_INT, rank - 1, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            }
            
            for(int j = 0; j < sub_size; j++) {
                if(wt[i] <= (j + rank * sub_size)) { // could put item in knapsack
                    int previous_value = sub_lower_row[j];
                    int replace_items;
                    int index = j - wt[i];
                    if (index >= 0)
                        replace_items = val[i] + sub_lower_row[index];
                    else
                        replace_items = val[i] + aux_row[100 + index];

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

        MPI_Barrier(MPI_COMM_WORLD);
    }

    MPI_Gatherv(sub_upper_row, sub_size, MPI_INT, upper_row, counts, offsets, MPI_INT, ROOT, MPI_COMM_WORLD);
    MPI_Gatherv(sub_lower_row, sub_size, MPI_INT, lower_row, counts, offsets, MPI_INT, ROOT, MPI_COMM_WORLD);

    if (rank == ROOT) {
    
        int retval;

        if ((n % 2) == 0)
            retval = upper_row[W]; 
        else
            retval = lower_row[W];

        printf("%d\n", retval);
    
        free(lower_row);
        free(upper_row);
    }

        free(val);
        free(wt);

    MPI_Finalize();
    fclose(arq);

    return 0;
}