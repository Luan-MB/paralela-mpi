#include <stdio.h>
#include <stdlib.h>
#include <omp.h>


// A utility function that returns
// maximum of two integers
int max(int a, int b) {return (a > b) ? a : b;}

int knapsack(int MAXIMUM_CAPACITY, int wt[], int val[], int n)
{

    // start both arrays with 0s
    int *upper_row = calloc(MAXIMUM_CAPACITY + 1, sizeof(int));
    int *lower_row = calloc(MAXIMUM_CAPACITY + 1, sizeof(int));

    // evaluate item `i`
    for(int i = 0; i < n; i++) {
        if ((i % 2) == 0) {
            
            for(int j = 0; j <= MAXIMUM_CAPACITY; j++) {
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
            }

        } else {

            for(int j = 0; j <= MAXIMUM_CAPACITY; j++) {
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

    int result;

    if ((n % 2) == 0) {
        result =  upper_row[MAXIMUM_CAPACITY];
    } else {
        result =  lower_row[MAXIMUM_CAPACITY];
    }

    free(upper_row);
    free(lower_row);

    return result;
}

// Driver program to test above function
int main()
{
	int n, W;

	scanf("%d %d", &n, &W);
	int *val = (int*) calloc(n, sizeof(int));
	int *wt = (int*) calloc(n, sizeof(int));

	int i;
	for (i = 0; i < n; ++i) {
		scanf("%d %d", &(val[i]), &(wt[i])); 
	}
    
    printf("%d\n", knapsack(W, wt, val, n));

    free(val);
    free(wt);
    return 0;
}
