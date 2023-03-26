//
//
//

#ifndef AREPO_CPP_FUNCTIONS_H
#define AREPO_CPP_FUNCTIONS_H

#include <lapacke.h>
#ifdef __cplusplus
extern "C"{
#endif

// functions defined in residual_distribution_solver.cpp, which may be used in other C files
void cpp_hello_world();
void compute_residuals(tessellation*);
int boundary_triangle_check_responsibility_thistask(tessellation* , int);
int boundary_triangle_compare(tessellation*T, int, int);
void rd_test_func(tessellation*);
void apply_FluxRD_list(void);
int FluxRD_list_data_compare(const void *, const void *);
void apply_DualArea_list(void);
int DualArea_list_data_compare(const void *, const void *);
lapack_int mat_inv(double *A, unsigned n);
//void write_residual(char* fname,int Ndt_thistask, int* thistask_triangles,double* Residual_List);
// functions defined in C files, which may be used in the CPP file residual_distribution_solver.cpp
//void write_delaunay_triangulation(tessellation *, char *, int, int);



#ifdef __cplusplus
}
#endif






#endif  // AREPO_CPP_FUNCTIONS_H
