//
//
//

#ifndef AREPO_CPP_FUNCTIONS_H
#define AREPO_CPP_FUNCTIONS_H


#ifdef __cplusplus
extern "C"{
#endif

// functions defined in residual_distribution_solver.cpp, which may be used in other C files
void cpp_hello_world();
void compute_residuals(tessellation*);
int boundary_triangle_check_responsibility_thistask(tessellation* , int);
int boundary_triangle_compare(tessellation*T, int, int);
void rd_test_func(tessellation*);
void triangle_vertex_do_time_extrapolation(struct state_primitive *delta, struct state_primitive *st,struct grad_data *grad, double dt_Extrapolation);
void triangle_vertex_add_extrapolation(struct state_primitive *delta, struct state_primitive *st);
void apply_FluxRD_list(void);
int FluxRD_list_data_compare(const void *, const void *);
void apply_DualArea_list(void);
int DualArea_list_data_compare(const void *, const void *);

// functions defined in C files, which may be used in the CPP file residual_distribution_solver.cpp
// void mpi_printf(const char *fmt, ...);
//void write_delaunay_triangulation(tessellation *, char *, int, int);



#ifdef __cplusplus
}
#endif






#endif  // AREPO_CPP_FUNCTIONS_H
