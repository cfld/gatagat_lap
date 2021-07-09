#pragma GCC diagnostic ignored "-Wunused-result"

#include <stdlib.h>
#include <iostream>
#include <chrono>

#include "lapjv.h"

// #define VERBOSE

int_t n_rows;
int_t n_cols;
int_t n_nnz;

int_t* x;
int_t* y;
int_t* indptr;
int_t* indices;
cost_t* data;

using namespace std::chrono;

void load_data(std::string inpath) {
    FILE *ptr;
    ptr = fopen(inpath.c_str(), "rb");
    
    fread(&n_rows, sizeof(int_t), 1, ptr);
    fread(&n_cols, sizeof(int_t), 1, ptr);
    fread(&n_nnz,  sizeof(int_t), 1, ptr);
    
    indptr  = (int_t*)  malloc((n_rows + 1) * sizeof(int_t));
    indices = (int_t*)  malloc(n_nnz        * sizeof(int_t));
    data    = (cost_t*) malloc(n_nnz        * sizeof(cost_t));
    
    fread(indptr,  sizeof(int_t),  n_rows + 1 , ptr);  // send directy to the memory since thats what the thing is.
    fread(indices, sizeof(int_t),  n_nnz      , ptr);
    fread(data,    sizeof(cost_t), n_nnz      , ptr);

#ifdef VERBOSE
    printf("n_rows : %d\n", n_rows);
    printf("n_cols : %d\n", n_cols);
    printf("n_nnz  : %d\n", n_nnz);
    
    for(int i = 0; i < 10; i++)
      std::cout << indptr[i] << " ";
    std::cout << std::endl;
    
    for(int i = 0; i < 10; i++)
      std::cout << indices[i] << " ";
    std::cout << std::endl;
    
    for(int i = 0; i < 10; i++)
      std::cout << data[i] << " ";
    std::cout << std::endl;
#endif
    
}

int main(int n_args, char** argument_array) {
    load_data(argument_array[1]);
    
    int_t* x = (int_t*)malloc(n_rows * sizeof(int_t));
    int_t* y = (int_t*)malloc(n_rows * sizeof(int_t));
    
    for(int i = 0; i < n_rows; i++) x[i] = 0;
    for(int i = 0; i < n_rows; i++) y[i] = 0;
    
    fp_t fp_version = FP_1;
    cost_t large    = 9999999;
    
    cost_t max_cost = 0;
    for(int offset = 0; offset < n_nnz; offset++) {
      if(data[offset] > max_cost) 
        max_cost = data[offset];
    }
    
    cost_t* idata = (cost_t*)malloc(n_nnz * sizeof(cost_t));
    for(int offset = 0; offset < n_nnz; offset++) {
      idata[offset] = max_cost - data[offset];
    }
    
    auto t1 = high_resolution_clock::now();
    lapmod_internal(n_rows, idata, indptr, indices, x, y, fp_version, large);
    auto elapsed = high_resolution_clock::now() - t1;
    long long ms = duration_cast<microseconds>(elapsed).count();
    
    cost_t cost = 0;
    for(int i = 0; i < n_rows; i++) {
      for(int offset = indptr[i]; offset < indptr[i + 1]; offset++) {
        int_t j = indices[offset];
        if(j == x[i]) {
          cost += data[offset];
        }
      }
    }
    
    printf("cost = %f | ms = %lld \n", cost, ms);
}