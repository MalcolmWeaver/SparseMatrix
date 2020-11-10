#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>
#include <string.h>
#include <math.h> /* must llink with -lm to have access to sqrt */


/*
typedef struct csr{ // recursive doulbing, then realloc to correct size
    double* values; // malloc [nnzs]
    unsigned int* col_idxs; // malloc nnzs
    unsigned int* row_idxs; // malloc num rows + 1
} csr;

typedef struct csc{
    double* values;
    unsigned int* row_idxs;
    unsigned int* col_idxs;
} csc;
*/
typedef struct value_with_col{
    double value;
    unsigned int col_idx;
} val_col;


typedef struct sparsematrix{
    unsigned int nnzs; /* first three are how many are initialized */
    unsigned int nrows; /* index of row that is currently being initialized ie last row */
    unsigned int ncols;
    unsigned int maxnnzs; /* how many are allocated */
    unsigned int maxprefixsumsz; /* 1 more than num of rows */
    unsigned int * row_prefix_sums;
    val_col* val_col_array;
} sparsematrix; /* 36 bytes */





sparsematrix* create_empty(unsigned int maxnnzs, int maxprefixsumsz);
void * delete(sparsematrix * matrix);
void print(sparsematrix * matrix);

void insert_value(sparsematrix * matrix, double new_value, unsigned int col_num);
void next_row_to_build(sparsematrix * matrix);

void increase_size_values_and_col_idxs(sparsematrix* matrix);
void increase_size_row_prefix_sums(sparsematrix* matrix);

void analyze(sparsematrix * matrix); /*The CSR format saves on memory only when NNZ < (m (n − 1) − 1) / 2*/


sparsematrix * scalar_multiply(sparsematrix * matrix, double scalar_multiplier); /* Assume overflow of doubles does not occur */
sparsematrix * scalar_divide(sparsematrix * matrix, double scalar_divisor); /* error check for 0. also, could be same func as mul */

/* Matrix "Math" */
/* not sure if the following 3 are actually a thing (broadcasting?) */
sparsematrix * scalar_add(sparsematrix * matrix, double scalar_summand); /* add to a given location */
sparsematrix * scalar_matrixminus(sparsematrix * matrix, double scalar_substrahend); /* could be same func as add */
sparsematrix * scalar_minusmatrix(sparsematrix * matrix, double scalar_minuend); /* could be add (multiply by -1) minuend */

/* Other Funcs */
void * empty(sparsematrix * matrix);
sparsematrix* copy(sparsematrix* matrix);


sparsemtrxzip* convert_to_zip(sparsematrix* matrix);
sparsematrix* convert_from_zip(sparsemtrxzip* zip);
void sort_by_col_idx(sparsematrix* matrix); /* could also remove zeros */



/*Read/write*/

int write_matrix(const char * filename, sparsematrix * matrix); // format determined by ext
int read_matrix(const char * filename, sparsematrix * matrix);

/* not required */

void sort_by_col_idx(sparsematrix* matrix); /* could also remove zeros */


/* Other Funcs */

void remove_zeros(sparsematrix* matrix);

sparsematrix * transpose(sparsematrix * matrix);
/* THIS IS THE MOST IMPORTANT FUNCTION. IT WILL BE USED TO MULTIPLY. WE WILL STORE A SEPERATE TRANSPOSED MATRIX. */

// check dims for next 3
sparsematrix * matrix_multiply(sparsematrix * matrix_factor1, sparsematrix * matrix_factor2);
sparsematrix * matrix_add(sparsematrix * matrix_summand1, sparsematrix * matrix_summand2); /* assert same sizes */
sparsematrix * matrix_subtract(sparsematrix * matrix_minuend, sparsematrix * matrix_substrahend); // could be comp of add  and mul -1





/* Test Funcs */

int Test_read_and_write(char* input_file, char* output_file, char* correct_read, char* correct_write);
int Test_create_matrix(char* stdin_inputs_file, sparsematrix* correct_create); // not sure of input format/type
int Test_Scalar_Math(sparsematrix* matrix, double scalar, sparsematrix** correct_results);
int Test_Matrix_Math(sparsematrix* matrix1, sparsematrix* matrix2, sparsematrix** correct_results);
int Test_print_analyze(sparsematrix* matrix, sparsematrix** correct_results);

void test_print_matrix(sparsematrix* matrix){
    matrix->ncols = 5;
    matrix->nrows = 4;
    matrix->nnzs = 7;
    matrix->val_col_array[0].value = 1;
    matrix->val_col_array[1].value = 2;
    matrix->val_col_array[2].value = 3;
    matrix->val_col_array[3].value = 4;
    matrix->val_col_array[4].value = 5;
    matrix->val_col_array[5].value = 6;
    matrix->val_col_array[6].value = 7;

    matrix->val_col_array[0].col_idx[0] = 0;
    matrix->val_col_array[1].col_idx[1] = 3;
    matrix->val_col_array[2].col_idx[2] = 1;
    matrix->val_col_array[3].col_idx[3] = 4;
    matrix->val_col_array[4].col_idx[4] = 2;
    matrix->val_col_array[5].col_idx[5] = 0;
    matrix->val_col_array[6].col_idx[6] = 3;

    matrix->row_prefix_sums[0] = 0;
    matrix->row_prefix_sums[1] = 2;
    matrix->row_prefix_sums[2] = 4;
    matrix->row_prefix_sums[3] = 5;
    matrix->row_prefix_sums[4] = 7;
}


// Build from scratch

sparsematrix* create_empty(unsigned int maxnnzs, int maxprefixsumsz){
    sparsematrix* new = (sparsematrix *) malloc(sizeof(*new));
    new->nnzs = 0;
    new->nrows = 0; /*this could be set to 0 to avoid printing out shit*/
    new->ncols = 0;
    new->maxnnzs = maxnnzs;
    new->val_col_array = (val_col *) calloc(new->maxnnzs, sizeof(*(new->val_col_array)));
    new->val_col_array.col_idx = (unsigned int *) calloc(new->maxnnzs, sizeof(*(new->val_col_array.col_idx)));
    new->maxprefixsumsz = maxprefixsumsz;
    new->row_prefix_sums = (unsigned int*) calloc(1+new->maxprefixsumsz, sizeof(*(new->row_prefix_sums)));
    /* here we malloc new itself, and within new are 3 mallocs to free */
    return new;
}

void * delete(sparsematrix * matrix){ // free pointer
    free(matrix->val_col_array);
    free(matrix->val_col_array.col_idx);
    free(matrix->row_prefix_sums);
    free(matrix);
}

void print(sparsematrix * matrix){ /* add a test to not print zeros */
    int i, j;
    printf("CSR FORMAT: \n\n");
    printf("Number of Rows: %u, Number of Columns: %u\n", matrix->nrows, matrix->ncols);
    for(i = 0; i < matrix->nrows; ++i){
        /* take a slice from csr.row_idx[i] to csr.row_idx[i+1] */
        assert(matrix->row_prefix_sums[i] <= matrix->row_prefix_sums[i+1]); /* correct prefix sum layout */
        for(j = (matrix->row_prefix_sums)[i]; j < (matrix->row_prefix_sums)[i+1]; ++j){
            printf("%f ", matrix->val_col_array[j].value);
        }
        printf("\n");
    }
}

void insert_value(sparsematrix * matrix, double new_value, unsigned int col_num){ /* column is given as zero indexed. This function merely adds to the end of the list and incs the last prefix sum */
    /*must start at the row where you want the first value */
    if(!matrix->nrows){
        next_row_to_build(matrix);
    }
    if((matrix->nrows + 1) > matrix->maxprefixsumsz){ /* row size of one must have 2 ps, row size of 2 must have at least 3 etc. */
        increase_size_row_prefix_sums(matrix);
    }

    if(col_num >= matrix->ncols){ /* can add a column further out than currently initialiuzed */
        matrix->ncols = (col_num + 1);
    }

    unsigned int i, col_repeat;
    char is_repeat = 0;

    for(i = matrix->row_prefix_sums[matrix->nrows - 1]; i < matrix->row_prefix_sums[matrix->nrows]; ++i){
        if(matrix->val_col_array[i].col_idxi] == col_num){
            col_repeat = i;
            is_repeat = 1;
            printf("Column repeat caught. Replacing.\n");
        }
    }
    if(is_repeat){
            matrix->val_col_array[col_repeat].value = new_value;
    }
    else{
    ++(matrix->nnzs);
    if(matrix->nnzs > matrix->maxnnzs){
        increase_size_values_and_col_idxs(matrix);
    }
    else{
        assert(matrix->row_prefix_sums[matrix->nrows - 1] == 0 || matrix->row_prefix_sums[matrix->nrows - 1]); /* ensure that there is at least one idx before current: prefix sum array must start with a 0 */
    }
    matrix->val_col_array[matrix->nnzs - 1].value = new_value; /* inserts at the back of list. it could ba valuable to sort if trying to add//multiply */
    matrix->val_col_array[matrix->nnzs - 1].col_idx = col_num;
    ++(matrix->row_prefix_sums[matrix->nrows]);

}
}
/* this function adds a nnz to arrnnzs, the nnz's col to arrcols, and 1 to the last in the cumulitive_row_counts_array
It updates ncols if it is the largest(do you allow an empty col?) col_num so far. Compare col_num and matrix-> numcols */

void next_row_to_build(sparsematrix * matrix)/* this adds a value (initialized to the previous or 0) to the new cumulitive_row_counts_array
(new mem allocation) */
{
    if(matrix->nrows + 1 >= matrix->maxprefixsumsz){
        increase_size_row_prefix_sums(matrix);
    }
    assert(matrix->row_prefix_sums[matrix->nrows+1] == 0);
    matrix->row_prefix_sums[matrix->nrows + 1] = matrix -> row_prefix_sums[matrix->nrows];
    matrix->nrows += 1;
}

void increase_size_values_and_col_idxs(sparsematrix* matrix)
{   unsigned int old_sz = matrix->maxnnzs;
    matrix->maxnnzs = 2*old_sz;
    matrix->val_col_array = (val_col*) realloc(matrix->val_col_array, matrix->maxnnzs * (sizeof(matrix->val_col_array[0].value)));
    memset((matrix->val_col_array + old_sz*sizeof(matrix->val_col_array[0].value)), 0, old_sz*sizeof(matrix->val_col_array[0].value));
    printf("Reallocated and set new memory to zero for values array.\n");
    matrix->val_col_array.col_idx = (unsigned int *) realloc(matrix->val_col_array.col_idx, matrix->maxnnzs * (sizeof(matrix->val_col_array[0].col_idx)));
    memset((matrix->val_col_array.col_idx + old_sz*sizeof(matrix->val_col_array[0].col_idx)), 0, old_sz*sizeof(matrix->val_col_array[0].col_idx));
    printf("Reallocated and set new memory to zero for column index array.\n");
}
/* uses calloc not malloc */

void increase_size_row_prefix_sums(sparsematrix* matrix)
{
    unsigned int old_sz = matrix->maxprefixsumsz;
    matrix->maxprefixsumsz = 2*old_sz;
    matrix->row_prefix_sums = (unsigned int *) realloc(matrix->row_prefix_sums, matrix->maxprefixsumsz * (sizeof(matrix->row_prefix_sums[0])));
    memset((matrix->row_prefix_sums + old_sz*sizeof(matrix->row_prefix_sums[0])), 0, old_sz*sizeof(matrix->row_prefix_sums[0]));
    printf("Reallocated and set new memory to zero for row prefix sum array.\n");

}
void analyze(sparsematrix * matrix) /*The CSR format saves on memory only when NNZ < (m (n − 1) − 1) / 2*/
{
    printf("\nAnalysis of Matrix:\n\n");
    if(!matrix->nrows || !matrix->ncols || !matrix->nnzs){
        printf("Matrix has zero elements.\n");
    }
    else{/* we will need to call transpose and sort/clean for this to work ideally */
        double density, row_avg_nnzs, col_avg_nnzs;
        double min = matrix-> values[0], max = matrix->val_col_array[0].value, sum = 0, variation, stdev, mean, total_nnzs = 0;
        unsigned int i, j;
        for(i = 0; i < matrix->nrows; ++i){
            for(j = matrix->row_prefix_sums[i]; j < matrix->row_prefix_sums[i + 1]; ++j){
                if(matrix->val_col_array[j].value < min){min = matrix->val_col_array[j].value;}
                else if(matrix->val_col_array[j].value > max){max = matrix->val_col_array[j].value;}
                sum += matrix->val_col_array[j].value;
                ++total_nnzs;
            }
        }

        assert(total_nnzs == matrix -> row_prefix_sums[matrix->nrows]);
        assert(total_nnzs);
        mean = sum /total_nnzs;

        for(i = 0; i < matrix->nrows; ++i){
            for(j = matrix->row_prefix_sums[i]; j < matrix->row_prefix_sums[i + 1]; ++j){
                variation += (mean - matrix->val_col_array[j].value) * (mean - matrix->val_col_array[j].value);
            }
        }
        density = total_nnzs / (matrix->nrows*matrix->ncols);
        row_avg_nnzs = total_nnzs / (matrix->nrows);
        col_avg_nnzs = total_nnzs / (matrix->ncols);
        stdev = sqrt(variation / total_nnzs);
        printf("Sparse Matrix Density: %f\n Row Density: %f\n Column Density: %f\n General Mean: %f\n Standard Deviation: %f\n Min: %f\n Max: %f\n", density, row_avg_nnzs, col_avg_nnzs, mean, stdev, min, max);
        if(total_nnzs < (matrix->nrows*(matrix->ncols - 1) - 1) / 2.0){
            printf("Using CSR Saves space.\n");
        }
        else{printf("Using CSR DOES NOT save space.\n");}
    }
}

sparsematrix * scalar_multiply(sparsematrix * matrix, double scalar_multiplier){
    /*assuming no overflow */
    for(int i = 0; i < matrix->nnzs; ++i){
        matrix->val_col_array[i].value *= scalar_multiplier;
    }
    return matrix;
} /*check for overflow (bits 52-63)? */

sparsematrix * scalar_divide(sparsematrix * matrix, double scalar_divisor){
    if(scalar_divisor == 0){
        printf("Invalid divisor. Returning NULL.");
        return NULL;
    }
    else{
        for(int i = 0; i < matrix->nnzs; ++i){
            matrix->val_col_array[i].value /= scalar_divisor;
        }
        return matrix;
    }
} /* error check for 0. also, could be same func as mul */

/* broadcast to nnzs, check overflow? could also check if sum is less than the bigger of the two */
sparsematrix * scalar_add(sparsematrix * matrix, double scalar_summand){
    /*assuming no overflow */
    for(int i = 0; i < matrix->nnzs; ++i){
        matrix->val_col_array[i].value += scalar_summand;
    }
    return matrix;

}
sparsematrix * scalar_matrixminus(sparsematrix * matrix, double scalar_substrahend){
    /*assuming no overflow */
    for(int i = 0; i < matrix->nnzs; ++i){
        matrix->val_col_array[i].value -= scalar_substrahend;
    }
    return matrix;

}
sparsematrix * scalar_minusmatrix(sparsematrix * matrix, double scalar_minuend){
    /*assuming no overflow */
    for(int i = 0; i < matrix->nnzs; ++i){
        matrix->val_col_array[i].value = scalar_minuend - matrix->val_col_array[i].value;
    }
    return matrix;

}

sparsematrix * copy(sparsematrix * matrix){
    int i;
    sparsematrix* new = (sparsematrix *) malloc(sizeof(*new));
    new->nnzs = matrix->nnzs;
    new->nrows = matrix->nrows; /*this could be set to 0 to avoid printing out shit*/
    new->ncols = matrix->ncols;
    new->maxnnzs = matrix->maxnnzs;
    new->val_col_array = (val_col *) calloc(new->maxnnzs, sizeof(*(new->val_col_array)));
    new->val_col_array.col_idx = (unsigned int *) calloc(new->maxnnzs, sizeof(*(new->val_col_array.col_idx)));
    for(i = 0; i < new->nnzs; ++i){
        new->val_col_array[i].value = matrix->val_col_array[i].value;
        new->val_col_array[i].col_idx = matrix->val_col_array[i].col_idx;
    }
    new->maxprefixsumsz = matrix->maxprefixsumsz;
    new->row_prefix_sums = (unsigned int*) calloc(1+new->maxprefixsumsz, sizeof(*(new->row_prefix_sums)));
    for(i = 0; i < new->maxprefixsumsz; ++i){
        new->row_prefix_sums[i] = matrix->row_prefix_sums[i];
    }
    /* here we malloc new itself, and within new are 3 mallocs to free */
    return new;
}

void * empty(sparsematrix * matrix){
    memset(matrix->val_col_array, 0, matrix->nnzs*sizeof(matrix->val_col_array[0].value));
    memset(matrix->val_col_array.col_idx, 0, matrix->nnzs*sizeof(matrix->val_col_array[].col_idx0]));
    memset(matrix->row_prefix_sums, 0, matrix->maxprefixsumsz*sizeof(matrix->row_prefix_sums[0]));
    printf("Memory for Matrix set to 0s.\n");
}


int main(){
    sparsematrix* matrix = create_empty(7, 5);
    test_print_matrix(matrix);
    matrix = scalar_divide(matrix, 11);
    matrix = scalar_minusmatrix(matrix, 100);
    print(matrix);
    sparsematrix * new = copy(matrix);
    empty(new);
    print(new);
    analyze(matrix);
    delete(matrix);

/*
    matrix = create_empty(0,0);
    print(matrix);
    analyze(matrix);
    delete(matrix);
*/

    /*sparsematrix* m3x3 = create_empty(3, 4);
    insert_value(m3x3, 12.34, 0);
    next_row_to_build(m3x3);
    insert_value(m3x3, 1.2345, 1);
    next_row_to_build(m3x3);
    insert_value(m3x3, 678.123243, 2);
    insert_value(m3x3, 987654.321, 1);
    insert_value(m3x3, 3.00, 3);
    print(m3x3);
    insert_value(m3x3, 4, 4);
    next_row_to_build(m3x3);
    insert_value(m3x3, 3, 5);
    next_row_to_build(m3x3);
    insert_value(m3x3, 300, 3);

    print(m3x3);
    delete(m3x3);
    */
    return 0;
}
