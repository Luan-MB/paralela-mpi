#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define ROOT 0
#define STD_TAG 0

// A utility function that returns
// maximum of two integers
int max(int a, int b) {return (a > b) ? a : b;}

/* int** get_matrix(int rows, int columns) {   
    int **mat;
    int i;
    
    // for each line
    mat = (int**) calloc(rows, sizeof (int*));
    
    mat[0] = (int*) calloc(rows * columns, sizeof (int));

    // set up pointers to rows
    for (i = 1; i < rows; i++)
        mat[i] = mat[0] + i * columns;

    return mat;
}

void free_matrix(int** mat) {
    free(mat[0]);
    free(mat);
} */
/* 
int new_knapsack(int W, int wt[], int val[], int n) {

    int *V = calloc(W + 1, sizeof(int));

     for (int i = 1; i < n + 1; i++) {
        for (int w = W; w >= 0; w--) {
  
            if (wt[i - 1] <= w)
                // finding the maximum value
                V[w] = max(V[w],
                            V[w - wt[i - 1]] + val[i - 1]);
        }
    }
    return V[W];
} */

int knapsack(int W, int wt[], int val[], int n)
{
    // Matrix-based solution
    int *upper_row = calloc (W + 1, sizeof(int));
    int *lower_row = calloc (W + 1, sizeof(int));

    // V Stores, for each (1 + i, j), the best profit for a knapscak
    // of capacity `j` considering every item k such that (0 <= k < i)
    int i, j;

    // evaluate item `i`
    for(i = 0; i < n; i++) {
        for(j = 1; j <= W; j++) {
            if ((i % 2) == 0) {
                if(wt[i] <= j) { // could put item in knapsack
                    int previous_value = upper_row[j];
                    int replace_items = val[i] + upper_row[j - wt[i]];

                    // is it better to keep what we already got,
                    // or is it better to swap whatever we have in the bag that weights up to `j`
                    // and put item `i`?
                    lower_row[j]= max(previous_value, replace_items);
                }
                else {
                    // can't put item `i`
                    lower_row[j] = upper_row[j];
                }
            } else {
                if(wt[i] <= j) { // could put item in knapsack
                    int previous_value = lower_row[j];
                    int replace_items = val[i] + lower_row[j - wt[i]];

                    // is it better to keep what we already got,
                    // or is it better to swap whatever we have in the bag that weights up to `j`
                    // and put item `i`?
                    upper_row[j]= max(previous_value, replace_items);
                }
                else {
                    // can't put item `i`
                    upper_row[j] = lower_row[j];
                }
            }
        }
    }

    int retval;

    if ((n % 2) == 0)
        retval = lower_row[W]; 
    else
        retval = upper_row[W];
    
    free(lower_row);
    free(upper_row);
    
    return retval;
}

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

    int sub_size;
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
    }

    MPI_Bcast(&sub_size, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(val, n, MPI_INT, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(wt, n, MPI_INT, ROOT, MPI_COMM_WORLD);

    int *sub_upper_row = malloc(sub_size * sizeof(int));
    int *sub_lower_row = malloc(sub_size * sizeof(int));

    MPI_Scatter(upper_row, sub_size, MPI_INT, sub_upper_row, sub_size, MPI_INT, ROOT, MPI_COMM_WORLD);
    MPI_Scatter(lower_row, sub_size, MPI_INT, sub_lower_row, sub_size, MPI_INT, ROOT, MPI_COMM_WORLD);

    int *aux_row = malloc(sub_size * sizeof(int));

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

    MPI_Gather(sub_upper_row, sub_size, MPI_INT, upper_row, sub_size, MPI_INT, ROOT, MPI_COMM_WORLD);
    MPI_Gather(sub_lower_row, sub_size, MPI_INT, lower_row, sub_size, MPI_INT, ROOT, MPI_COMM_WORLD);

    if (rank == ROOT) {
    
        int retval;

        if ((n % 2) == 0)
            retval = lower_row[W]; 
        else
            retval = upper_row[W];

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