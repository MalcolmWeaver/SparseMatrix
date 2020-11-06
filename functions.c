#include <stdio.h>
#include <stdlib.h>


typedef struct csr{ // recursive doulbing, then realloc to correct size
    double* values; // malloc [num_nonzeros]
    unsigned int* col_idxs; // malloc num_nonzeros
    unsigned int* row_idxs; // malloc num rows + 1
} csr;

typedef struct csc{
    double* values;
    unsigned int* row_idxs;
    unsigned int* col_idxs;
} csc;

typedef struct sparsematrix{
    unsigned int num_nonzeros;
    unsigned int num_rows;
    unsigned int num_cols;
    csr csr;
    csc csc;
} sparsematrix;


// Read/write

int write_matrix(const char * filename, sparsematrix * matrix); // format determined by ext
int read_matrix(const char * filename, sparsematrix * matrix);


int add_value_to_row(sparsematrix * matrix, double new_value, unsigned int col_num, unsigned int row_num);
/* this function adds a nnz to arrnnzs, the nnz's col to arrcols, and 1 to the last in the cumulitive_row_counts_array
It updates num_cols if it is the largest(do you allow an empty col?) col_num so far. Compare col_num and matrix-> numcols

insert to sll sorted by col. */
int next_row_to_add(sparsematrix * matrix, unsigned int row_num);
/* this adds a value (initialized to the previous or 0) to the new cumulitive_row_counts_array
(new mem allocation) */


/* Matrix "Math" */

sparsematrix * scalar_multiply(sparsematrix * matrix, double scalar_multiplier); /* error check for overflow? */
sparsematrix * scalar_divide(sparsematrix * matrix, double scalar_divider); /* error check for 0. also, could be same func as mul */
/* not sure if the following 3 are actually a thing (broadcasting?) */
sparsematrix * scalar_add(sparsematrix * matrix, double scalar_summand);
sparsematrix * scalar_matrixminus(sparsematrix * matrix, double scalar_substrahend); /* could be same func as add */
sparsematrix * scalar_minusmatrix(sparsematrix * matrix, double scalar_minuend); /* could be add (multiply by -1) minuend */

sparsematrix * empty(sparsematrix * matrix);
void * delete(sparsematrix * matrix); // free pointer
sparsematrix * transpose(sparsematrix * matrix);

// check dims for next 3
sparsematrix * matrix_multiply(sparsematrix * matrix_factor1, sparsematrix * matrix_factor2);
sparsematrix * matrix_add(sparsematrix * matrix_summand1, sparsematrix * matrix_summand2);
sparsematrix * matrix_subtract(sparsematrix * matrix_minuend, sparsematrix * matrix_substrahend); // could be comp of add  and mul -1

/* Other Funcs */


void analyze(sparsematrix * matrix);

/* not required */
void free_matrix(sparsematrix*);

/* Test Funcs */

int Test_read_and_write(char* input_file, char* output_file, char* correct_read, char* correct_write);
int Test_create_matrix(char* stdin_inputs_file, sparsematrix* correct_create); // not sure of input format/type
int Test_Scalar_Math(sparsematrix* matrix, double scalar, sparsematrix** correct_results);
int Test_Matrix_Math(sparsematrix* matrix1, sparsematrix* matrix2, sparsematrix** correct_results);
int Test_print_analyze(sparsematrix* matrix, sparsematrix** correct_results);

// Build from scratch

sparsematrix* create_empty(unsigned int num_nonzeros, unsigned int num_rows, unsigned int num_cols){
    sparsematrix* new = (sparsematrix *) malloc(sizeof(*new));
    new->num_nonzeros = num_nonzeros;
    new->num_rows = num_rows;
    new->num_cols = num_cols;
    new->csr.values = (double *) malloc(sizeof(*(new->csr.values))*new->num_nonzeros);
    new->csr.col_idxs = (unsigned int *) malloc(sizeof(*(new->csr.col_idxs))*new->num_nonzeros);
    new->csr.row_idxs = (unsigned int*) malloc(sizeof(*(new->csr.row_idxs))*(1+new->num_rows));

    new->csc.values = (double *) malloc(sizeof(*(new->csc.values))*new->num_nonzeros);
    new->csc.col_idxs = (unsigned int *) malloc(sizeof(*(new->csc.col_idxs))*new->num_nonzeros);
    new->csc.row_idxs = (unsigned int*) malloc(sizeof(*(new->csc.row_idxs))*(1+new->num_cols));
    /* here we malloc new itself, and within new are 6 mallocs to free */
}

void print(sparsematrix * matrix){
    int i, j;
    printf("CSR FORMAT: \n\n");
    for(i = 0; i < matrix->num_rows; ++i){
        /* take a slice from csr.row_idx[i] to csr.row_idx[i+1] */
        for(j = (matrix->csr.row_idxs)[i]; j < (matrix->csr.row_idxs)[i+1]; ++j){
            printf("R: %u C: %u Value: %f", i+1, (matrix->csr.col_idxs)[j], matrix->csr.values[j]);
        }
    }
}

int main(){

    sparsematrix* matrix = create_empty(3, 4, 4);
    return 0;
}
