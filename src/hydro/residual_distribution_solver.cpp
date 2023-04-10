
#include <lapacke.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <iostream>

#include "../main/allvars.h"
#include "../main/cpp_functions.h"
#include "../main/proto.h"
#include "../mesh/mesh.h"
#include "../mesh/voronoi/voronoi.h"

using namespace std;

// static struct flux_list_data
//{
//  int task, index;
//  double dM, dP[3];
//#ifdef MHD
//  double dB[3];
//#endif /* #ifdef MHD */
//
//#ifndef ISOTHERM_EQS
//  double dEnergy;
//#endif /* #ifndef ISOTHERM_EQS */
//#ifdef MAXSCALARS
//  double dConservedScalars[MAXSCALARS];
//#endif /* #ifdef MAXSCALARS */
//} * FluxList;

// extern "C" void mpi_printf(const char *fmt, ...);
void cpp_hello_world() { cout << "This is a test function in C++!" << endl; }

#ifdef RESIDUAL_DISTRIBUTION

static struct FluxRD_list_data
{
  int task, index;
  double dMass_Dual;
  double dMomentum_Dual[3];
  double dEnergy_Dual;

} * FluxRD_list;

static struct DualArea_list_data
{
  int task, index;
  double DualArea;
} * DualArea_list;

static int N_DualArea_export, Max_N_FluxRD_export, N_FluxRD_export;
struct triangle_normals *tri_normals_list;
extern struct primexch *PrimExch;
extern struct grad_data *GradExch;

/*! \brief Compute residuals of triangles/tetrahedra and distribute them.
 *  This function is used for Residual Distribution hydrodynamics method, which is equivalent to
 *  compute_interface_fluxes() in the original AREPO's Finite Volume approach.
 *
 *  \param[in] T Pointer to tessellation.
 *
 *  \return void
 */
void compute_residuals(tessellation *T)
{
#ifdef NOHYDRO
  return;
#endif /* #ifdef NOHYDRO */
  TIMER_START(CPU_RESIDUAL_DISTRIBUTION);
  point *DP = T->DP;
  tetra *DT = T->DT;
  int i, j = 0, k = 0, p;
  int Ndp = T->Ndp;
  int Ndt = T->Ndt;
  char *DT_label;
  DT_label      = (char *)mymalloc_movable(&DT_label, "DT_label",
                                      Ndt * sizeof(char)); /* array of labels of the triangles, l = local, b = boundary, o = other*/
  int Ndt_local = 0, Ndt_boundary = 0, Ndt_other = 0, Ndt_boundary_thistask = 0;

  /* classification of Delaunay triangles*/
  for(i = 0; i < Ndt; i++)
    {
      int pmin = imin_array(DT[i].p, DIMS + 1);
      int pmax = imax_array(DT[i].p, DIMS + 1);

      if(pmin >= 0 && pmax <= NumGas - 1) /* all vertices are local (excluding local ghost)*/
        {
          DT_label[i] = 'l';
          Ndt_local += 1;
        }
      else if(pmin >= 0 && pmin <= NumGas - 1) /*at least one vertex is local (excluding local ghost), but not all of them are local*/
        {
          DT_label[i] = 'b';
          Ndt_boundary += 1;
          if(boundary_triangle_check_responsibility_thistask(T, i) == ThisTask)
            {
              DT_label[i] = 't';
              Ndt_boundary_thistask += 1;
            }
        }
      else
        {
          DT_label[i] = 'o';
          Ndt_other += 1;
        }
    }

  int *local_triangles    = (int *)mymalloc_movable(&local_triangles, "local_triangles", Ndt_local * sizeof(int));
  int *boundary_triangles = (int *)mymalloc_movable(&boundary_triangles, "boundary_triangles", Ndt_boundary_thistask * sizeof(int));

  for(i = 0; i < Ndt; i++)
    {
      if(DT_label[i] == 'l')
        {
          local_triangles[j] = i;
          j += 1;
        }
      else if(DT_label[i] == 't')
        {
          boundary_triangles[k] = i;
          k += 1;
        }
    }
  /* check uniqueness of boundary triangles of this task*/
  int Ndt_boundary_thistask_repeated = 0;
  for(i = 0; i < Ndt_boundary_thistask; i++)
    {
      for(j = 0; j < i; j++)
        {
          if(boundary_triangle_compare(T, boundary_triangles[i], boundary_triangles[j]))
            {
              DT_label[boundary_triangles[i]] = 'r';
              boundary_triangles[i]           = -1;
              Ndt_boundary_thistask_repeated += 1;
            }
        }
    }

  int Ndt_thistask        = Ndt_local + Ndt_boundary_thistask - Ndt_boundary_thistask_repeated;
  int *thistask_triangles = (int *)mymalloc_movable(&thistask_triangles, "thistask_triangles", Ndt_thistask * sizeof(int));
  for(i = 0; i < Ndt_local; i++)
    {
      thistask_triangles[i] = local_triangles[i];
    }
  j = Ndt_local;
  for(i = 0; i < Ndt_boundary_thistask; i++)
    {
      if(boundary_triangles[i] >= 0)
        {
          thistask_triangles[j] = boundary_triangles[i];
          j++;
        }
    }

  /* compute residual: calculate normal vectors and triangular area; assign dual area to vertices
   tri_normals_list is a shorter list compared to DT because it only keeps necessary triangles for
   this task (Ndt_thistask). tri_normal_list[i]  <-->  DT[thistask_triangles[i]]
   */
  tri_normals_list = (struct triangle_normals *)mymalloc_movable(&tri_normals_list, "tri_normals_list",
                                                                 Ndt_thistask * sizeof(struct triangle_normals));

  for(i = 0; i < Ndt_thistask; i++)
    {
#ifdef TWODIMS
      triangle_get_normals_area(T, thistask_triangles[i], &tri_normals_list[i]);
#endif
    }

  for(i = 0; i < NumGas; i++)
    {
      SphP[i].DualArea = 0.0;
    }

  N_DualArea_export = 0;
  for(i = 0; i < Ndt_thistask; i++)
    {
      for(j = 0; j < DIMS + 1; j++)
        {
          if(DP[DT[thistask_triangles[i]].p[j]].task == ThisTask)
            {
              int SphP_index = DP[DT[thistask_triangles[i]].p[j]].index;
              //              if(SphP_index < 0)
              //                continue;   not necessary here, since the triangles we selected should not contain external points
              if(SphP_index > NumGas)
                SphP_index -= NumGas;

              SphP[SphP_index].DualArea += tri_normals_list[i].area / (DIMS + 1);
            }
          else
            {
              N_DualArea_export += 1;
            }
        }
    }

  DualArea_list = (struct DualArea_list_data *)mymalloc_movable(&DualArea_list, "DualArea_list",
                                                                N_DualArea_export * sizeof(struct DualArea_list_data));
  k             = 0;
  for(i = Ndt_local; i < Ndt_thistask; i++)
    {
      for(j = 0; j < DIMS + 1; j++)
        {
          if(DP[DT[thistask_triangles[i]].p[j]].task != ThisTask)
            {
              DualArea_list[k].task     = DP[DT[thistask_triangles[i]].p[j]].task;
              DualArea_list[k].index    = DP[DT[thistask_triangles[i]].p[j]].originalindex;
              DualArea_list[k].DualArea = tri_normals_list[i].area / (DIMS + 1);
              k += 1;
            }
        }
    }

  apply_DualArea_list();
  //  MPI_Barrier(MPI_COMM_WORLD);
  //  char triangulation_name[1024];
  //  snprintf(triangulation_name, 100, "%s/triangulation_dual_%03d", All.OutputDir, 0);
  //  write_delaunay_triangulation(T, triangulation_name, 0, NTask - 1);

  // debug
  //  double *Residual_List;
  //  Residual_List = (double *)mymalloc_movable(&Residual_List, "Residual_List", Ndt_thistask *4* sizeof(double));

  //
  //  face_dt = (((integertime)1) << timeBin) * All.Timebase_interval
  //

  Max_N_FluxRD_export = N_DualArea_export;
  N_FluxRD_export     = 0;
  FluxRD_list =
      (struct FluxRD_list_data *)mymalloc_movable(&FluxRD_list, "FluxRD_list", Max_N_FluxRD_export * sizeof(struct FluxRD_list_data));

  // main loop through triangles this task responsible for
  for(i = 0; i < Ndt_thistask; i++)
    {

      //get triangle timebin/timestep and skip inactive triangles
      int timebin_vertices[DIMS + 1];
      for(j = 0; j < DIMS + 1; j++)
        {
          if(DP[DT[thistask_triangles[i]].p[j]].task == ThisTask)
            {
              int SphP_index = DP[DT[thistask_triangles[i]].p[j]].index;
              if(SphP_index > NumGas)
                SphP_index -= NumGas;

              timebin_vertices[j] = P[SphP_index].TimeBinHydro;
            }
          else
            {
              int PrimExch_index  = DP[DT[thistask_triangles[i]].p[j]].index;
              timebin_vertices[j] = PrimExch[PrimExch_index].TimeBinHydro;
            }
        }

      bool is_active = false;
      int timebin_this_triangle = timebin_vertices[0];
      
      for(j = 0; j < DIMS + 1; j++)
        {
          if(TimeBinSynchronized[timebin_vertices[j]])
            is_active = true;

          if(timebin_vertices[j] < timebin_this_triangle)
            timebin_this_triangle = timebin_vertices[j];
        }

      if(is_active == false)
        continue;

      double triangle_dt = (((integertime)1) << timebin_this_triangle) * All.Timebase_interval;

      triangle_dt *= 0.5;  //RK2 half timestep

      // compute residual: set up initial states
      double U_fluid[DIMS + 1][DIMS + 2];  // specific conserved fluid variables
      double C_sound[DIMS + 1];
      double Pressure[DIMS + 1];


      for(j = 0; j < DIMS + 1; j++)
      {
        if(DP[DT[thistask_triangles[i]].p[j]].task == ThisTask)
        {
          int SphP_index = DP[DT[thistask_triangles[i]].p[j]].index;
          if(SphP_index > NumGas)
            SphP_index -= NumGas;

          //extrapolation to current time
          struct grad_data *grad = &SphP[SphP_index].Grad;

          struct state_primitive vertex_state;
          struct state_primitive delta_time;
          vertex_state.rho = SphP[SphP_index].Density;
          vertex_state.press = SphP[SphP_index].Pressure;
          vertex_state.velx = P[SphP_index].Vel[0];
          vertex_state.vely = P[SphP_index].Vel[1];
          vertex_state.velz = P[SphP_index].Vel[2];

          double dt_Extrapolation = All.Time - SphP[SphP_index].TimeLastPrimUpdate;

          triangle_vertex_do_time_extrapolation(&delta_time,&vertex_state,grad,dt_Extrapolation);
          triangle_vertex_add_extrapolation(&delta_time,&vertex_state);


#ifdef TWODIMS
          U_fluid[j][0] = vertex_state.rho;
          U_fluid[j][1] = U_fluid[j][0]*vertex_state.velx;
          U_fluid[j][2] = U_fluid[j][0]*vertex_state.vely;
          Pressure[j] = vertex_state.press;
          double kinetic_energy = 0.0;
          kinetic_energy += pow(vertex_state.velx,2)+pow(vertex_state.vely,2);
          kinetic_energy *= 0.5 * U_fluid[j][0];
          U_fluid[j][DIMS + 1] = kinetic_energy + Pressure[j] / (GAMMA_MINUS1);
#endif
          if(U_fluid[j][DIMS + 1] <= 0)
          {
            printf("sph energy <= 0 error %d %d   %f  %f  %f  %f\n", ThisTask, thistask_triangles[i], U_fluid[j][DIMS + 1],
                   SphP[SphP_index].Energy,kinetic_energy,Pressure[j]);
            terminate_program("sph energy <= 0 error.");
          }

#ifdef TREE_BASED_TIMESTEPS
          C_sound[j] = SphP[SphP_index].Csnd;
#else
          C_sound[j] = get_sound_speed(SphP_index);
#endif
          if(C_sound[j] <= 0)
          {
            terminate_program("sph Cs <= 0 error!");
          }
        }
        else
        {
          int PrimExch_index  = DP[DT[thistask_triangles[i]].p[j]].index;


          struct grad_data *grad = &GradExch[PrimExch_index];
          struct state_primitive vertex_state;
          struct state_primitive delta_time;

          vertex_state.rho = PrimExch[PrimExch_index].Density;
          vertex_state.press = PrimExch[PrimExch_index].Pressure;
          vertex_state.velx = PrimExch[PrimExch_index].VelGas[0];
          vertex_state.vely = PrimExch[PrimExch_index].VelGas[1];
          vertex_state.velz = PrimExch[PrimExch_index].VelGas[2];
          double dt_Extrapolation = All.Time - PrimExch[PrimExch_index].TimeLastPrimUpdate;

          triangle_vertex_do_time_extrapolation(&delta_time,&vertex_state,grad,dt_Extrapolation);
          triangle_vertex_add_extrapolation(&delta_time,&vertex_state);

          /*we should only use primitive variables to get U_fluid*/

#ifdef TWODIMS
          U_fluid[j][0] = vertex_state.rho;
          U_fluid[j][1] = U_fluid[j][0]*vertex_state.velx;
          U_fluid[j][2] = U_fluid[j][0]*vertex_state.vely;
          Pressure[j] = vertex_state.press;
          double kinetic_energy = 0.0;
          kinetic_energy += pow(vertex_state.velx,2)+pow(vertex_state.vely,2);
          kinetic_energy *= 0.5 * U_fluid[j][0];
          U_fluid[j][DIMS + 1] = kinetic_energy + Pressure[j] / (GAMMA_MINUS1);
#endif

          C_sound[j]           = PrimExch[PrimExch_index].Csnd;

          if(C_sound[j] <= 0)
          {
            terminate_program("primexch Cs <= 0 error!");
          }
          if(U_fluid[j][DIMS + 1] <= 0)
          {
            printf("primexch energy <= 0 error %d %d  %d %d %d    %f  %f\n", ThisTask, DT[thistask_triangles[i]].p[j],
                   DP[DT[thistask_triangles[i]].p[j]].task, DP[DT[thistask_triangles[i]].p[j]].originalindex,
                   DP[DT[thistask_triangles[i]].p[j]].ID, U_fluid[j][DIMS + 1], PrimExch[PrimExch_index].Energy);
            terminate_program("primexch energy <= 0 error.");
          }
        }
      }  // for(j = 0; j < DIMS + 1; j++) get fluid state for each vertex of this triangle








      // compute residuals: Roe Vector Z, Modified fluid state U_hat = \frac{ \partial{U(Z_avg)}}{ \partial Z} * Z
#ifdef TWODIMS  // for now we only consider 2D. 3D case needs to be included in the future
      double Z_Roe[3][4];
      double Z_avg[4];
      double Enthalpy[DIMS + 1];
      double N_X[3], N_Y[3], Mag[3];
      double U_hat[4][3];  // notice: now U_hat[4][3] and Z_Roe[3][4] We may unify their formats in the future

      for(j = 0; j < 3; j++)
        {
          Z_Roe[j][0] = sqrt(U_fluid[j][0]);
          Z_Roe[j][1] = U_fluid[j][1] / Z_Roe[j][0];
          Z_Roe[j][2] = U_fluid[j][2] / Z_Roe[j][0];
          Z_Roe[j][3] = (U_fluid[j][3] + Pressure[j]) / Z_Roe[j][0];

          N_X[j]      = tri_normals_list[i].normal[j][0];
          N_Y[j]      = tri_normals_list[i].normal[j][1];
          Mag[j]      = tri_normals_list[i].mag[j];
          Enthalpy[j] = (U_fluid[j][DIMS + 1] + Pressure[j]) / U_fluid[j][0];  // specific enthalpy = (rho E + p)/rho
        }

      for(k = 0; k < 4; k++)
        Z_avg[k] = (Z_Roe[0][k] + Z_Roe[1][k] + Z_Roe[2][k]) / 3.0;

      for(j = 0; j < 3; j++)
        {
          U_hat[0][j] = 2.0 * Z_avg[0] * Z_Roe[j][0];
          U_hat[1][j] = Z_avg[1] * Z_Roe[j][0] + Z_avg[0] * Z_Roe[j][1];
          U_hat[2][j] = Z_avg[2] * Z_Roe[j][0] + Z_avg[0] * Z_Roe[j][2];
          U_hat[3][j] = (Z_avg[3] * Z_Roe[j][0] + GAMMA_MINUS1 * Z_avg[1] * Z_Roe[j][1] + GAMMA_MINUS1 * Z_avg[2] * Z_Roe[j][2] +
                         Z_avg[0] * Z_Roe[j][3]) /
                        (GAMMA);
        }

      // compute residual: Construct average state for element
      double sum_sqrt_rho = sqrt(U_fluid[0][0]) + sqrt(U_fluid[1][0]) + sqrt(U_fluid[2][0]);
      double rho_avg      = pow(sum_sqrt_rho / 3.0, 2);
      double velx_avg = 0, vely_avg = 0, h_avg = 0;
      for(j = 0; j < 3; j++)
        {
          velx_avg += sqrt(U_fluid[j][0]) * U_fluid[j][1] / U_fluid[j][0];
          vely_avg += sqrt(U_fluid[j][0]) * U_fluid[j][2] / U_fluid[j][0];
          h_avg += sqrt(U_fluid[j][0]) * Enthalpy[j];
        }
      velx_avg /= sum_sqrt_rho;
      vely_avg /= sum_sqrt_rho;
      h_avg /= sum_sqrt_rho;

      double Cs_avg = sqrt(GAMMA_MINUS1 * (h_avg - (velx_avg * velx_avg + vely_avg * vely_avg) / 2.0));
      if(isnan(Cs_avg))
        {
          printf("cs avg nan error! %d %d    %f    %f %f %f   %f %f %f   %f %f %f\n", ThisTask, thistask_triangles[i], h_avg,
                 Enthalpy[0], Enthalpy[1], Enthalpy[2], U_fluid[0][DIMS + 1], U_fluid[1][DIMS + 1], U_fluid[2][DIMS + 1], Pressure[0],
                 Pressure[1], Pressure[2]);
          terminate_program("Cs avg nan.")
        }

      // compute residual: Reassign variables to local equivalents
      double velx_c  = velx_avg / Cs_avg;
      double vely_c  = vely_avg / Cs_avg;
      double h_c     = h_avg / Cs_avg;
      double alpha   = GAMMA_MINUS1 * (velx_avg * velx_avg + vely_avg * vely_avg) / 2.0;
      double alpha_c = alpha / Cs_avg;

      // compute residual: Calculate K+,K- and K matrices for each vertex
      double vel_dot_n;
      double Lambda[3][4], Lambda_plus[3][4], Lambda_minus[3][4];
      double Value1, Value2, Value3, Value4, Value12, Value123;
      double Kmatrix[4][4][3][3];  // Kmatrix[4][4][j=0,1,2(vertices)][p=0(K+),1(K-),2(K)]
      int kfull = 2, kplus = 0, kminus = 1;

      for(j = 0; j < 3; j++)
        {
          vel_dot_n    = velx_avg * N_X[j] + vely_avg * N_Y[j];
          Lambda[j][0] = vel_dot_n + Cs_avg;
          Lambda[j][1] = vel_dot_n - Cs_avg;
          Lambda[j][2] = vel_dot_n;
          Lambda[j][3] = vel_dot_n;

          for(k = 0; k < 4; k++)
            {
              Lambda_plus[j][k]  = dmax(0.0, Lambda[j][k]);
              Lambda_minus[j][k] = dmin(0.0, Lambda[j][k]);
            }
          // fill in Kmatrix[4][4][j][p]
          for(p = 0; p < 3; p++)
            {
              if(p == 0)
                {  // Identify and select positive eigenvalues
                  Value1 = Lambda_plus[j][0];
                  Value2 = Lambda_plus[j][1];
                  Value3 = Lambda_plus[j][2];
                  Value4 = Lambda_plus[j][3];
                }
              else if(p == 1)
                {  // Identify and select negative eigenvalues
                  Value1 = Lambda_minus[j][0];
                  Value2 = Lambda_minus[j][1];
                  Value3 = Lambda_minus[j][2];
                  Value4 = Lambda_minus[j][3];
                }
              else
                {  // Select all eigenvalues
                  Value1 = Lambda[j][0];
                  Value2 = Lambda[j][1];
                  Value3 = Lambda[j][2];
                  Value4 = Lambda[j][3];
                }

              Value12  = (Value1 - Value2) / 2.0;
              Value123 = (Value1 + Value2 - 2.0 * Value3) / 2.0;

              Kmatrix[0][0][j][p] = 0.5 * Mag[j] * (alpha_c * Value123 / Cs_avg - vel_dot_n * Value12 / Cs_avg + Value3);
              Kmatrix[0][1][j][p] = 0.5 * Mag[j] * (-1.0 * GAMMA_MINUS1 * velx_c * Value123 / Cs_avg + N_X[j] * Value12 / Cs_avg);
              Kmatrix[0][2][j][p] = 0.5 * Mag[j] * (-1.0 * GAMMA_MINUS1 * vely_c * Value123 / Cs_avg + N_Y[j] * Value12 / Cs_avg);
              Kmatrix[0][3][j][p] = 0.5 * Mag[j] * (GAMMA_MINUS1 * Value123 / (Cs_avg * Cs_avg));

              Kmatrix[1][0][j][p] =
                  0.5 * Mag[j] *
                  ((alpha_c * velx_c - vel_dot_n * N_X[j]) * Value123 + (alpha_c * N_X[j] - velx_c * vel_dot_n) * Value12);
              Kmatrix[1][1][j][p] = 0.5 * Mag[j] *
                                    ((N_X[j] * N_X[j] - GAMMA_MINUS1 * velx_c * velx_c) * Value123 -
                                     ((GAMMA - 2.0) * velx_c * N_X[j] * Value12) + Value3);
              Kmatrix[1][2][j][p] = 0.5 * Mag[j] *
                                    ((N_X[j] * N_Y[j] - GAMMA_MINUS1 * velx_c * vely_c) * Value123 +
                                     (velx_c * N_Y[j] - GAMMA_MINUS1 * vely_c * N_X[j]) * Value12);
              Kmatrix[1][3][j][p] =
                  0.5 * Mag[j] * (GAMMA_MINUS1 * velx_c * Value123 / Cs_avg + GAMMA_MINUS1 * N_X[j] * Value12 / Cs_avg);

              Kmatrix[2][0][j][p] =
                  0.5 * Mag[j] *
                  ((alpha_c * vely_c - vel_dot_n * N_Y[j]) * Value123 + (alpha_c * N_Y[j] - vely_c * vel_dot_n) * Value12);
              Kmatrix[2][1][j][p] = 0.5 * Mag[j] *
                                    ((N_X[j] * N_Y[j] - GAMMA_MINUS1 * velx_c * vely_c) * Value123 +
                                     (vely_c * N_X[j] - GAMMA_MINUS1 * velx_c * N_Y[j]) * Value12);
              Kmatrix[2][2][j][p] = 0.5 * Mag[j] *
                                    ((N_Y[j] * N_Y[j] - GAMMA_MINUS1 * vely_c * vely_c) * Value123 -
                                     ((GAMMA - 2.0) * vely_c * N_Y[j] * Value12) + Value3);
              Kmatrix[2][3][j][p] =
                  0.5 * Mag[j] * (GAMMA_MINUS1 * vely_c * Value123 / Cs_avg + GAMMA_MINUS1 * N_Y[j] * Value12 / Cs_avg);

              Kmatrix[3][0][j][p] =
                  0.5 * Mag[j] * ((alpha_c * h_c - vel_dot_n * vel_dot_n) * Value123 + vel_dot_n * (alpha_c - h_c) * Value12);
              Kmatrix[3][1][j][p] = 0.5 * Mag[j] *
                                    ((vel_dot_n * N_X[j] - velx_avg - alpha_c * velx_c) * Value123 +
                                     (h_c * N_X[j] - GAMMA_MINUS1 * velx_c * vel_dot_n) * Value12);
              Kmatrix[3][2][j][p] = 0.5 * Mag[j] *
                                    ((vel_dot_n * N_Y[j] - vely_avg - alpha_c * vely_c) * Value123 +
                                     (h_c * N_Y[j] - GAMMA_MINUS1 * vely_c * vel_dot_n) * Value12);
              Kmatrix[3][3][j][p] =
                  0.5 * Mag[j] * (GAMMA_MINUS1 * h_c * Value123 / Cs_avg + GAMMA_MINUS1 * vel_dot_n * Value12 / Cs_avg + Value3);
            }
        }
      // compute residual: get residual Phi
      double Phi[4];
      for(k = 0; k < 4; k++)
        {
          Phi[k] = 0.0;
          for(j = 0; j < 3; j++)
            {
              Phi[k] += Kmatrix[k][0][j][kfull] * U_hat[0][j] + Kmatrix[k][1][j][kfull] * U_hat[1][j] +
                        Kmatrix[k][2][j][kfull] * U_hat[2][j] + Kmatrix[k][3][j][kfull] * U_hat[3][j];
            }
        }
      // debug
      //      for(k=0;k<4;k++){
      //          if(isnan(Phi[k])){
      //              printf("phi nan error! %d %d   %f %f %f %f\n",ThisTask,i,vel_dot_n,Cs_avg,Kmatrix[0][0][0][kfull],Value123);
      //            }
      //        }

      //      for(k = 0; k < 4; k++)
      //        {
      //          Residual_List[i*4+k] = Phi[k];
      //        }

      double Kmatrix_minus_sum[4][4];

      for(k = 0; k < 4; k++)
        {
          for(p = 0; p < 4; p++)
            {
              Kmatrix_minus_sum[k][p] = 0.0;
              for(j = 0; j < 3; j++)
                {
                  Kmatrix_minus_sum[k][p] += Kmatrix[k][p][j][kminus];
                }
            }
        }

      mat_inv(&Kmatrix_minus_sum[0][0], 4);

      // residual distribution:
      double Flux_RD[4][3];
#ifdef B_SCHEME
      double Flux_LDA[4][3], Flux_N[4][3];
#endif

#if(defined(LDA_SCHEME) || defined(B_SCHEME))
      double BetaLDA[4][4][3]; /* distributed residual phi_j= BetaLDA_j * phi_total */
      for(k = 0; k < 4; k++)
        {
          for(p = 0; p < 4; p++)
            {
              for(j = 0; j < 3; j++)
                {
                  BetaLDA[k][p][j] =
                      -1.0 * (Kmatrix[k][0][j][kplus] * Kmatrix_minus_sum[0][p] + Kmatrix[k][1][j][kplus] * Kmatrix_minus_sum[1][p] +
                              Kmatrix[k][2][j][kplus] * Kmatrix_minus_sum[2][p] + Kmatrix[k][3][j][kplus] * Kmatrix_minus_sum[3][p]);
                }
            }
        }

      for(k = 0; k < 4; k++)
        {
          for(j = 0; j < 3; j++)
            {
#ifdef LDA_SCHEME
              Flux_RD[k][j] =
                  BetaLDA[k][0][j] * Phi[0] + BetaLDA[k][1][j] * Phi[1] + BetaLDA[k][2][j] * Phi[2] + BetaLDA[k][3][j] * Phi[3];
#else
              Flux_LDA[k][j] =
                  BetaLDA[k][0][j] * Phi[0] + BetaLDA[k][1][j] * Phi[1] + BetaLDA[k][2][j] * Phi[2] + BetaLDA[k][3][j] * Phi[3];
#endif
            }
        }

#endif  // LDA scheme or B scheme

#if(defined(N_SCHEME) || defined(B_SCHEME))

      double Bracket[4][3];
      double KU_Sum[4];

      for(k = 0; k < 4; k++)
        {
          KU_Sum[k] = 0.0;
          for(j = 0; j < 3; j++)
            {
              KU_Sum[k] += Kmatrix[k][0][j][kminus] * U_hat[0][j] + Kmatrix[k][1][j][kminus] * U_hat[1][j] +
                           Kmatrix[k][2][j][kminus] * U_hat[2][j] + Kmatrix[k][3][j][kminus] * U_hat[3][j];
            }
        }

      for(k = 0; k < 4; k++)
        {
          for(j = 0; j < 3; j++)
            {
              Bracket[k][j] = U_hat[k][j] - (Kmatrix_minus_sum[k][0] * KU_Sum[0] + Kmatrix_minus_sum[k][1] * KU_Sum[1] +
                                             Kmatrix_minus_sum[k][2] * KU_Sum[2] + Kmatrix_minus_sum[k][3] * KU_Sum[3]);
            }
        }

      for(k = 0; k < 4; k++)
        {
          for(j = 0; j < 3; j++)
            {
#ifdef N_SCHEME
              Flux_RD[k][j] = Kmatrix[k][0][j][kplus] * Bracket[0][j] + Kmatrix[k][1][j][kplus] * Bracket[1][j] +
                              Kmatrix[k][2][j][kplus] * Bracket[2][j] + Kmatrix[k][3][j][kplus] * Bracket[3][j];
#else
              Flux_N[k][j] = Kmatrix[k][0][j][kplus] * Bracket[0][j] + Kmatrix[k][1][j][kplus] * Bracket[1][j] +
                             Kmatrix[k][2][j][kplus] * Bracket[2][j] + Kmatrix[k][3][j][kplus] * Bracket[3][j];
#endif
            }
        }
#endif  // N scheme or B scheme

#ifdef B_SCHEME
      double Sum_Flux_N[4];
      double Theta_E[4];

      for(k = 0; k < 4; k++)
        {
          Sum_Flux_N[k] = abs(Flux_N[k][0]) + abs(Flux_N[k][1]) + abs(Flux_N[k][2]);
          if(Sum_Flux_N[i] == 0)
            {
              Theta_E[k] = 0.0;
            }
          else
            {
              Theta_E[k] = abs(Phi[k]) / Sum_Flux_N[k];
            }
          Flux_RD[k][0] = Theta_E[k] * Flux_N[k][0] + (1.0 - Theta_E[k]) * Flux_LDA[k][0];
          Flux_RD[k][1] = Theta_E[k] * Flux_N[k][1] + (1.0 - Theta_E[k]) * Flux_LDA[k][1];
          Flux_RD[k][2] = Theta_E[k] * Flux_N[k][2] + (1.0 - Theta_E[k]) * Flux_LDA[k][2];
        }

#endif  // B scheme

      /*use residuals to update fluid state of local points or export to other tasks*/

      for(j = 0; j < 3; j++)
        {
          if(DP[DT[thistask_triangles[i]].p[j]].task == ThisTask)
            {
              int SphP_index = DP[DT[thistask_triangles[i]].p[j]].index;
              //              if(SphP_index < 0)
              //                continue;   not necessary here, since the triangles we selected should not contain external points
              if(SphP_index > NumGas)
                SphP_index -= NumGas;

              int P_index = SphP_index;
              P[P_index].Mass += SphP[SphP_index].Volume * (-1.0) * triangle_dt * Flux_RD[0][j] / SphP[SphP_index].DualArea;
              SphP[SphP_index].Momentum[0] += SphP[SphP_index].Volume * (-1.0) * triangle_dt * Flux_RD[1][j] / SphP[SphP_index].DualArea;
              SphP[SphP_index].Momentum[1] += SphP[SphP_index].Volume * (-1.0) * triangle_dt * Flux_RD[2][j] / SphP[SphP_index].DualArea;
              SphP[SphP_index].Energy += SphP[SphP_index].Volume * (-1.0) * triangle_dt * Flux_RD[3][j] / SphP[SphP_index].DualArea;
            }
          else
            {
              int PrimExch_index = DP[DT[thistask_triangles[i]].p[j]].index;

              FluxRD_list[N_FluxRD_export].task  = DP[DT[thistask_triangles[i]].p[j]].task;
              FluxRD_list[N_FluxRD_export].index = DP[DT[thistask_triangles[i]].p[j]].originalindex;

              FluxRD_list[N_FluxRD_export].dMass_Dual        = PrimExch[PrimExch_index].Volume * (-1.0) * triangle_dt * Flux_RD[0][j];
              FluxRD_list[N_FluxRD_export].dMomentum_Dual[0] = PrimExch[PrimExch_index].Volume * (-1.0) * triangle_dt * Flux_RD[1][j];
              FluxRD_list[N_FluxRD_export].dMomentum_Dual[1] = PrimExch[PrimExch_index].Volume * (-1.0) * triangle_dt * Flux_RD[2][j];
              FluxRD_list[N_FluxRD_export].dMomentum_Dual[2] = 0.0;
              FluxRD_list[N_FluxRD_export].dEnergy_Dual      = PrimExch[PrimExch_index].Volume * (-1.0) * triangle_dt * Flux_RD[3][j];

              N_FluxRD_export += 1;
            }
        }

#endif  // TWO_DIMS
    }  // for loop of triangles, i= 0~ Ndt_thistask

  apply_FluxRD_list();

  // debug
  //  char res_filename[100];
  //  snprintf(res_filename, 100, "%s/residual_%03d", All.OutputDir, 0);
  //  write_residual(res_filename,Ndt_thistask, thistask_triangles,Residual_List);

  //  myfree_movable(Residual_List);
  myfree_movable(FluxRD_list);
  myfree_movable(DualArea_list);
  myfree_movable(tri_normals_list);
  myfree_movable(thistask_triangles);
  myfree_movable(boundary_triangles);
  myfree_movable(local_triangles);
  myfree_movable(DT_label);

  TIMER_STOP(CPU_RESIDUAL_DISTRIBUTION);
}

int boundary_triangle_check_responsibility_thistask(tessellation *T, int DTindex)
{
  point *DP = T->DP;
  tetra *DT = T->DT;
  int tasks[DIMS + 1]; /* tasks of the vertices in this triangle */

  int *task_count; /* counts of each task*/
  task_count = (int *)mymalloc("task_count", NTask * sizeof(int));
  for(int i = 0; i < NTask; i++)
    {
      task_count[i] = 0;
    }

  for(int i = 0; i < DIMS + 1; i++)
    {
      tasks[i] = DP[DT[DTindex].p[i]].task;
      task_count[tasks[i]] += 1;
    }
  int responsible_task = arg_imax(task_count, NTask);

  myfree(task_count);
  return responsible_task;
}

/*same triangle: return 1; different triangles: return 0 */
int boundary_triangle_compare(tessellation *T, int DTindex1, int DTindex2)
{
  if(DTindex1 < 0 || DTindex2 < 0)
    return 0;
  point *DP = T->DP;
  tetra *DT = T->DT;
  MyIDType particleIDlist1[DIMS + 1], particleIDlist2[DIMS + 1];
  for(int i = 0; i < DIMS + 1; i++)
    {
      particleIDlist1[i] = DP[DT[DTindex1].p[i]].ID;
      particleIDlist2[i] = DP[DT[DTindex2].p[i]].ID;
    }

  qsort(particleIDlist1, DIMS + 1, sizeof(int), system_compare_int);
  qsort(particleIDlist2, DIMS + 1, sizeof(int), system_compare_int);

  for(int i = 0; i < DIMS + 1; i++)
    {
      if(particleIDlist1[i] != particleIDlist2[i])
        return 0;
    }

  return 1;
}

void rd_test_func(tessellation *T)
{
  //    cout <<"test"<<endl;
  mpi_printf("numsph: %d %d %d\n", NumGas, All.MaxPart, All.MaxPartSph);

  //  output triangulation & access data
  char triangulation_name[1024];
  snprintf(triangulation_name, 100, "%s/triangulation_%03d", All.OutputDir, 0);
  write_delaunay_triangulation(T, triangulation_name, 0, 2);
  MPI_Barrier(MPI_COMM_WORLD);

  // transfer residual between tasks
  double delta_pressure = 100.0;
  SphP[10].Pressure += delta_pressure;

  point *DP = T->DP;
  //  int i = 950;
  //  int pindex = DP[i].index;
  //  if (DP[i].task != ThisTask){
  //    printf("debug: test PrimExch %d %d %d %d %d %f\n",ThisTask, DP[i].task, DP[i].index, DP[i].originalindex, DP[i].ID,
  //    PrimExch[pindex].Pressure);
  //  }

  /*
  FluxListRD = (struct fluxRD_list_data *)mymalloc_movable(&FluxListRD, "FluxListRD", 1 * sizeof(struct fluxRD_list_data));

  if(ThisTask == 0)
    {
      FluxListRD[0].task      = 1;
      FluxListRD[0].index     = 15;
      FluxListRD[0].dPressure = 500;
    }
  else if(ThisTask == 1)
    {
      FluxListRD[0].task      = 2;
      FluxListRD[0].index     = 20;
      FluxListRD[0].dPressure = 600;
    }
  else if(ThisTask == 2)
    {
      FluxListRD[0].task      = 0;
      FluxListRD[0].index     = 25;
      FluxListRD[0].dPressure = 700;
    }

  if(ThisTask == 0)
    {
      cout << DP[25].ID << " " << DP[25].index << " " << SphP[DP[25].index].Pressure << endl;
    }

  apply_fluxRD_list();

  if(ThisTask == 0)
    {
      cout << DP[25].ID << " " << DP[25].index << " " << SphP[DP[25].index].Pressure << endl;
    }

  MPI_Barrier(MPI_COMM_WORLD);
  snprintf(triangulation_name, 100, "%s/triangulation_new_%03d", All.OutputDir, 0);
  write_delaunay_triangulation(T, triangulation_name, 0, 2);
  MPI_Barrier(MPI_COMM_WORLD);
  myfree_movable(FluxListRD);

   */
}

void triangle_vertex_do_time_extrapolation(struct state_primitive *delta, struct state_primitive *st, struct grad_data *grad, double dt_Extrapolation){
  if(st->rho <= 0)
    return;

  delta->rho = -dt_Extrapolation * (st->velx * grad->drho[0] + st->rho * grad->dvel[0][0] + st->vely * grad->drho[1] +
                           st->rho * grad->dvel[1][1] + st->velz * grad->drho[2] + st->rho * grad->dvel[2][2]);

  delta->velx = -dt_Extrapolation * (1.0 / st->rho * grad->dpress[0] + st->velx * grad->dvel[0][0] + st->vely * grad->dvel[0][1] +
                            st->velz * grad->dvel[0][2]);

  delta->vely = -dt_Extrapolation * (1.0 / st->rho * grad->dpress[1] + st->velx * grad->dvel[1][0] + st->vely * grad->dvel[1][1] +
                            st->velz * grad->dvel[1][2]);

  delta->velz = -dt_Extrapolation * (1.0 / st->rho * grad->dpress[2] + st->velx * grad->dvel[2][0] + st->vely * grad->dvel[2][1] +
                            st->velz * grad->dvel[2][2]);

  delta->press = -dt_Extrapolation * (GAMMA * st->press * (grad->dvel[0][0] + grad->dvel[1][1] + grad->dvel[2][2]) +
                             st->velx * grad->dpress[0] + st->vely * grad->dpress[1] + st->velz * grad->dpress[2]);

}

void triangle_vertex_add_extrapolation(struct state_primitive *delta_time, struct state_primitive *st){
  if(st->rho <= 0)
    return;

  if(st->rho + delta_time->rho || st->press + delta_time->press< 0)
    return;

  st->rho += delta_time->rho;
  st->velx += delta_time->velx;
  st->vely += delta_time->vely;
  st->velz += delta_time->velz;
  st->press += delta_time->press;

}


void apply_FluxRD_list(void)
{
  int i, j, p, nimport, ngrp, recvTask;

  /* now exchange the flux-list and apply it when needed */

  mysort(FluxRD_list, N_FluxRD_export, sizeof(struct FluxRD_list_data), FluxRD_list_data_compare);

  for(j = 0; j < NTask; j++)
    Send_count[j] = 0;

  for(i = 0; i < N_FluxRD_export; i++)
    Send_count[FluxRD_list[i].task]++;

  if(Send_count[ThisTask] > 0)
    terminate_program("Send_count[ThisTask]");

  MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

  for(j = 0, nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
    {
      nimport += Recv_count[j];

      if(j > 0)
        {
          Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
          Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
        }
    }

  struct FluxRD_list_data *FluxListGet = (struct FluxRD_list_data *)mymalloc("FluxListGet", nimport * sizeof(struct FluxRD_list_data));

  /* exchange particle data */
  for(ngrp = 0; ngrp < (1 << PTask); ngrp++)
    {
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
            {
              /* get the particles */
              MPI_Sendrecv(&FluxRD_list[Send_offset[recvTask]], Send_count[recvTask] * sizeof(struct FluxRD_list_data), MPI_BYTE,
                           recvTask, TAG_DENS_A, &FluxListGet[Recv_offset[recvTask]],
                           Recv_count[recvTask] * sizeof(struct FluxRD_list_data), MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD,
                           MPI_STATUS_IGNORE);
            }
        }
    }

  /* apply the fluxes */

  for(i = 0; i < nimport; i++)
    {
      p = FluxListGet[i].index;

      P[p].Mass += FluxListGet[i].dMass_Dual / SphP[p].DualArea;
      SphP[p].Momentum[0] += FluxListGet[i].dMomentum_Dual[0] / SphP[p].DualArea;
      SphP[p].Momentum[1] += FluxListGet[i].dMomentum_Dual[1] / SphP[p].DualArea;
      SphP[p].Momentum[2] += FluxListGet[i].dMomentum_Dual[2] / SphP[p].DualArea;
      SphP[p].Energy += FluxListGet[i].dEnergy_Dual / SphP[p].DualArea;
    }
  myfree(FluxListGet);
}

int FluxRD_list_data_compare(const void *a, const void *b)
{
  if(((struct FluxRD_list_data *)a)->task < (((struct FluxRD_list_data *)b)->task))
    return -1;

  if(((struct FluxRD_list_data *)a)->task > (((struct FluxRD_list_data *)b)->task))
    return +1;

  return 0;
}

void apply_DualArea_list(void)
{
  int i, j, p, nimport, ngrp, recvTask;
#if defined(MAXSCALARS)
  int k;
#endif /* #if defined(MAXSCALARS) */

  /* now exchange the flux-list and apply it when needed */

  mysort(DualArea_list, N_DualArea_export, sizeof(struct DualArea_list_data), DualArea_list_data_compare);

  for(j = 0; j < NTask; j++)
    Send_count[j] = 0;

  for(i = 0; i < N_DualArea_export; i++)
    {
      Send_count[DualArea_list[i].task]++;
      if(DualArea_list[i].task == ThisTask)
        {
          printf("bug: thistask, i_inDualAreaList %d %d %d\n", ThisTask, i, N_DualArea_export);
        }
    }

  if(Send_count[ThisTask] > 0)
    terminate_program("Send_count[ThisTask]");

  MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

  for(j = 0, nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
    {
      nimport += Recv_count[j];

      if(j > 0)
        {
          Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
          Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
        }
    }

  struct DualArea_list_data *DualAreaListGet =
      (struct DualArea_list_data *)mymalloc("DualAreaListGet", nimport * sizeof(struct DualArea_list_data));

  /* exchange particle data */
  for(ngrp = 0; ngrp < (1 << PTask); ngrp++)
    {
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
            {
              /* get the particles */
              MPI_Sendrecv(&DualArea_list[Send_offset[recvTask]], Send_count[recvTask] * sizeof(struct DualArea_list_data), MPI_BYTE,
                           recvTask, TAG_DENS_A, &DualAreaListGet[Recv_offset[recvTask]],
                           Recv_count[recvTask] * sizeof(struct DualArea_list_data), MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD,
                           MPI_STATUS_IGNORE);
            }
        }
    }

  for(i = 0; i < nimport; i++)
    {
      p = DualAreaListGet[i].index;
      SphP[p].DualArea += DualAreaListGet[i].DualArea;
    }

  myfree(DualAreaListGet);
}

int DualArea_list_data_compare(const void *a, const void *b)
{
  if(((struct DualArea_list_data *)a)->task < (((struct DualArea_list_data *)b)->task))
    return -1;

  if(((struct DualArea_list_data *)a)->task > (((struct DualArea_list_data *)b)->task))
    return +1;

  return 0;
}

lapack_int mat_inv(double *A, unsigned n)
{
  int ipiv[n + 1];
  lapack_int ret;

  ret = LAPACKE_dgetrf(LAPACK_COL_MAJOR, n, n, A, n, ipiv);

  //#ifdef DEBUG
  //  std::cout << "ret =\t" << ret << "\t(0 = done, <0 = illegal arguement, >0 = singular)" << std::endl;
  //#endif

  if(ret != 0)
    {
      std::cout << "B WARNING: MATRIX CANNOT BE INVERTED\t" << std::endl;
      for(int i = 0; i < n * n; ++i)
        {
          std::cout << A[i] << "\t";
          if((i + 1) % n == 0)
            {
              std::cout << std::endl;
            }
        }
      exit(0);
    }
  ret = LAPACKE_dgetri(LAPACK_COL_MAJOR, n, A, n, ipiv);

  return ret;
}

#endif  //#ifdef RESIDUAL_DISTRIBUTION

/*
 void write_residual(char* fname,int Ndt_thistask, int* thistask_triangles,double* Residual_List){

//   DT[thistask_triangles[i]] Res[i]
//
   FILE* fdtxt;
   char msg[1000], fname_new[1000];


   snprintf(fname_new,1000,"%s_%d.txt",fname,ThisTask);

   if(!(fdtxt = fopen(fname_new, "w")))
   {
     snprintf(msg,1000,"can't open file `%s' for writing snapshot.\n", fname);
     terminate_program("%s",msg);
   }

   fprintf(fdtxt,"thistask_triangles[i](i.e. triangle index of DT),  Residual \n");
   for (int i = 0;i< Ndt_thistask;i++){
       fprintf(fdtxt, "%d   %.10g %.10g %.10g %.10g \n",
               thistask_triangles[i],Residual_List[i*4],Residual_List[i*4+1],Residual_List[i*4+2],Residual_List[i*4+3]);
     }
   fclose(fdtxt);

 }

*/
