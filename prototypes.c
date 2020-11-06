#include <stdio.h>
#include <stdlib.h>

typedef struct sparsematrix{
    unsigned int num_nonzeros;
    unsigned int num_rows;
    unsigned int num_cols;
    // perhaps have 3 modes of storage: csr, csc, coo for matrix multiplication and transposition

} sparsematrix;


// Read/write

int write_matrix_to_ascii(const char * filename, sparsematrix * matrix); // format determined by ext
int read_ascii_to_matrix(const char * filename, sparsematrix * matrix);

// Build from scratch

sparsematrix * create_empty();

int add_value_to_row(sparsematrix * matrix, double new_value, unsigned int col_num, unsigned int row_num);
// I am thinking, we are using csr at this point.
// We may or may not end up using row_num: we might make it so that you cannot go back.
// This is because it is easiest if the array of nnzs grows in alignment with (n*row + col)
// this function adds a nnz to arrnnzs, the nnz's col to arrcols, and 1 to the last in the cumulitive_row_counts_array
// It updates num_cols if it is the largest(do you allow an empty col?) col_num so far. Compare col_num and matrix-> numcols
int next_row_to_add(sparsematrix * matrix, unsigned int row_num);
// this adds a value (initialized to the previous or 0) to the new cumulitive_row_counts_array
// (new mem allocation)


// Matrix "Math"

sparsematrix * scalar_multiply(sparsematrix * matrix, double scalar_multiplier); // error check for overflow?
sparsematrix * scalar_divide(sparsematrix * matrix, double scalar_divider); // error check for 0. also, could be same func as mul
// not sure if the following 3 are actually a thing (broadcasting?)
sparsematrix * scalar_add(sparsematrix * matrix, double scalar_summand);
sparsematrix * scalar_matrixminus(sparsematrix * matrix, double scalar_substrahend); // could be same func as add
sparsematrix * scalar_minusmatrix(sparsematrix * matrix, double scalar_minuend); // could be add (multiply by -1) minuend

sparsematrix * empty(sparsematrix * matrix);
void * delete(sparsematrix * matrix); // free pointer
sparsematrix * transpose(sparsematrix * matrix);

// check dims for next 3
sparsematrix * matrix_add(sparsematrix * matrix_factor1, sparsematrix * matrix_factor2);
sparsematrix * matrix_add(sparsematrix * matrix_summand1, sparsematrix * matrix_summand2);
sparsematrix * matrix_add(sparsematrix * matrix_minuend, sparsematrix * matrix_substrahend); // could be comp of add  and mul -1

//Other Funcs

void print(sparsematrix * matrix);
void analyze(sparsematrix * matrix);

//not required
void init_coo_from_csr(sparsematrix * matrix);
void init_csc_from_coo(sparsematrix * matrix); // makes transposing real easy

// Test Funcs

int Test_read_and_write(char* input_file, char* output_file, char* correct_read, char* correct_write);
int Test_create_matrix(char* stdin_inputs_file, sparsematrix* correct_create); // not sure of input format/type
int Test_Scalar_Math(sparsematrix* matrix, double scalar, sparsematrix** correct_results);
int Test_Matrix_Math(sparsematrix* matrix1, sparsematrix* matrix2, sparsematrix** correct_results);
int Test_print_analyze(sparsematrix* matrix, sparsematrix** correct_results);
