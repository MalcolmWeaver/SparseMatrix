#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>

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
typedef struct sparsematrix{
    unsigned int nnzs; /* first three are how many are initialized */
    unsigned int nrows; /* index of row that is currently being initialized ie last row */
    unsigned int ncols;
    unsigned int maxnnzs; /* how many are allocated */
    unsigned int maxprefixsumsz; /* 1 more than num of rows */
    double* values;
    unsigned int * col_idxs;
    unsigned int * row_prefix_sums;
} sparsematrix; /* 36 bytes */

void next_row_to_build(sparsematrix * matrix);
/*Read/write*/

int write_matrix(const char * filename, sparsematrix * matrix); // format determined by ext
int read_matrix(const char * filename, sparsematrix * matrix);


/* Matrix "Math" */

sparsematrix * scalar_multiply(sparsematrix * matrix, double scalar_multiplier); /* error check for overflow? */
sparsematrix * scalar_divide(sparsematrix * matrix, double scalar_divider); /* error check for 0. also, could be same func as mul */
/* not sure if the following 3 are actually a thing (broadcasting?) */
sparsematrix * scalar_add(sparsematrix * matrix, double scalar_summand); /* add to a given location */
sparsematrix * scalar_matrixminus(sparsematrix * matrix, double scalar_substrahend); /* could be same func as add */
sparsematrix * scalar_minusmatrix(sparsematrix * matrix, double scalar_minuend); /* could be add (multiply by -1) minuend */

sparsematrix * empty(sparsematrix * matrix);

sparsematrix * transpose(sparsematrix * matrix);
/* THIS IS THE MOST IMPORTANT FUNCTION. IT WILL BE USED TO MULTIPLY. WE WILL STORE A SEPERATE TRANSPOSED MATRIX. */

// check dims for next 3
sparsematrix * matrix_multiply(sparsematrix * matrix_factor1, sparsematrix * matrix_factor2);
sparsematrix * matrix_add(sparsematrix * matrix_summand1, sparsematrix * matrix_summand2); /* assert same sizes */
sparsematrix * matrix_subtract(sparsematrix * matrix_minuend, sparsematrix * matrix_substrahend); // could be comp of add  and mul -1

/* Other Funcs */


void analyze(sparsematrix * matrix); /*The CSR format saves on memory only when NNZ < (m (n − 1) − 1) / 2*/

/* not required */

void sort_by_col_idx(sparsematrix* matrix); /* could also remove zeros */

/* Test Funcs */

int Test_read_and_write(char* input_file, char* output_file, char* correct_read, char* correct_write);
int Test_create_matrix(char* stdin_inputs_file, sparsematrix* correct_create); // not sure of input format/type
int Test_Scalar_Math(sparsematrix* matrix, double scalar, sparsematrix** correct_results);
int Test_Matrix_Math(sparsematrix* matrix1, sparsematrix* matrix2, sparsematrix** correct_results);
int Test_print_analyze(sparsematrix* matrix, sparsematrix** correct_results);

void test_print_matrix(sparsematrix* matrix){
    assert(matrix->ncols == 5);
    assert(matrix->nrows == 4);
    assert(matrix->nnzs == 7);
    matrix->values[0] = 1;
    matrix->values[1] = 2;
    matrix->values[2] = 3;
    matrix->values[3] = 4;
    matrix->values[4] = 5;
    matrix->values[5] = 6;
    matrix->values[6] = 7;

    matrix->col_idxs[0] = 0;
    matrix->col_idxs[1] = 3;
    matrix->col_idxs[2] = 1;
    matrix->col_idxs[3] = 4;
    matrix->col_idxs[4] = 2;
    matrix->col_idxs[5] = 0;
    matrix->col_idxs[6] = 3;

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
    new->values = (double *) calloc(new->maxnnzs, sizeof(*(new->values)));
    new->col_idxs = (unsigned int *) calloc(new->maxnnzs, sizeof(*(new->col_idxs)));
    new->maxprefixsumsz = maxprefixsumsz;
    new->row_prefix_sums = (unsigned int*) calloc(1+new->maxprefixsumsz, sizeof(*(new->row_prefix_sums)));
    /* here we malloc new itself, and within new are 3 mallocs to free */
    return new;
}

void * delete(sparsematrix * matrix){ // free pointer
    free(matrix->values);
    free(matrix->col_idxs);
    free(matrix->row_prefix_sums);
    free(matrix);
}

void print(sparsematrix * matrix){ /* add a test to not print zeros */
    int i, j;
    printf("CSR FORMAT: \n\n");
    for(i = 0; i < matrix->nrows; ++i){
        /* take a slice from csr.row_idx[i] to csr.row_idx[i+1] */
        assert(matrix->row_prefix_sums[i] <= matrix->row_prefix_sums[i+1]); /* correct prefix sum layout */
        for(j = (matrix->row_prefix_sums)[i]; j < (matrix->row_prefix_sums)[i+1]; ++j){
            printf("R: %u C: %u Value: %f   ", i+1, (matrix->col_idxs)[j] + 1, matrix->values[j]);
        }
        printf("\n");
    }
}





void increase_size_values(sparsematrix* matrix){printf("Not capable of increasing VALUES yet\n");} /* use calloc not malloc */
void increase_size_col_idxs(sparsematrix* matrix){printf("Not capable of increasing COL_IDXS yet\n");}
void increase_size_row_prefix_sums(sparsematrix* matrix){printf("Not capable of increasing ROW_PREFIX_SUMS yet\n");}



int insert_value(sparsematrix * matrix, double new_value, unsigned int col_num){ /* column is given as zero indexed. This function merely adds to the end of the list and incs the last prefix sum */
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

    ++(matrix->nnzs);
    if(matrix->nnzs > matrix->maxnnzs){
        increase_size_values(matrix);
        increase_size_col_idxs(matrix);
    }
    else{
        assert(matrix->row_prefix_sums[matrix->nrows - 1] == 0 || matrix->row_prefix_sums[matrix->nrows - 1]); /* ensure that there is at least one idx before current: prefix sum array must start with a 0 */
    }
    matrix->values[matrix->nnzs - 1] = new_value; /* inserts at the back of list. it could ba valuable to sort if trying to add//multiply */
    matrix->col_idxs[matrix->nnzs - 1] = col_num;
    ++(matrix->row_prefix_sums[matrix->nrows]);
    return 0;
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
    matrix->row_prefix_sums[matrix->nrows] = matrix -> row_prefix_sums[matrix->nrows + 1];
    matrix->nrows += 1;
}


int main(){
    /*sparsematrix* matrix = create_empty(7, 4, 5);
    test_print_matrix(matrix);
    print(matrix);
    delete(matrix);
    matrix = create_empty(0,0,0);
    print(matrix);
    delete(matrix); */

    sparsematrix* m3x3 = create_empty(3, 4);
    // unsigned int current_row = 0;
    insert_value(m3x3, 12.34, 0);
    insert_value(m3x3, 12.3, 1);
    // current_row = next_row_to_build(m3x3, 1);
    // insert_value(m3x3, 12.3, 1, current_row);
    // current_row = next_row_to_build(m3x3, 2);
    // insert_value(m3x3, 12.3, 0, current_row);


    print(m3x3);
    delete(m3x3);
    return 0;
}
