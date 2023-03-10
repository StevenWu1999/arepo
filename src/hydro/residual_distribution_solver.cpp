
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

static struct fluxRD_list_data
{
  int task, index;
  double dPressure;
} * FluxListRD;

static struct DualArea_list_data
{
  int task, index;
  double DualArea;
} * DualArea_list;

static int Nflux, MaxNflux, N_DualArea_export;
struct triangle_normals *tri_normals_list;
extern struct primexch *PrimExch;
extern struct grad_data *GradExch;

void cpp_hello_world() { cout << "This is a test function in C++!" << endl; }

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

  // compute residual: calculate normal vectors and triangular area; assign dual area to vertices
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

              SphP[SphP_index].DualArea += tri_normals_list[thistask_triangles[i]].area / (DIMS + 1);
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
              DualArea_list[k].DualArea = tri_normals_list[thistask_triangles[i]].area / (DIMS + 1);
              k += 1;
            }
        }
    }

  apply_DualArea_list();
  MPI_Barrier(MPI_COMM_WORLD);
  char triangulation_name[1024];
  snprintf(triangulation_name, 100, "%s/triangulation_dual_%03d", All.OutputDir, 0);
  write_delaunay_triangulation(T, triangulation_name, 0, 2);



  for(i = 0; i < Ndt_thistask; i++)
    {
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

              U_fluid[j][0] = SphP[SphP_index].Density;
              for(k = 0; k < DIMS; k++)
                {
                  U_fluid[j][k + 1] =
                      SphP[SphP_index].Momentum[k] / SphP[SphP_index].Volume;  // or P[SphP_index].Vel * SphP[SphP_index].Density
                }
              U_fluid[j][DIMS + 1] = SphP[SphP_index].Energy / SphP[SphP_index].Volume;
#ifdef TREE_BASED_TIMESTEPS
              C_sound[j] = SphP[SphP_index].Csnd;
#else
              C_sound[j] = get_sound_speed(SphP_index);
#endif
              Pressure[j] = SphP[SphP_index].Pressure;
            }
          else
            {
              int PrimExch_index = DP[DT[thistask_triangles[i]].p[j]].index;
              U_fluid[j][0]      = PrimExch[PrimExch_index].Density;
              for(k = 0; k < DIMS; k++)
                {
                  U_fluid[j][k + 1] = PrimExch[PrimExch_index].VelGas[k] * U_fluid[j][0];
                }
              U_fluid[j][DIMS + 1] = PrimExch[PrimExch_index].Energy / PrimExch[PrimExch_index].Volume;
              C_sound[j]           = PrimExch[PrimExch_index].Csnd;

              Pressure[j] = PrimExch[PrimExch_index].Pressure;
            }
        }
        // compute residuals: Roe Vector Z, Modified fluid state U_hat = \frac{ \partial{U(Z_avg)}}{ \partial Z} * Z
#ifdef TWODIMS  // for now we only consider 2D. 3D case needs to be included in the future
      double Z_Roe[3][4];
      double Z_avg[4];
      double Enthalpy[DIMS + 1];
      double N_X[3], N_Y[3], Mag[3];
      double U_hat[3][4];

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
          U_hat[j][0] = 2.0 * Z_avg[0] * Z_Roe[j][0];
          U_hat[j][1] = Z_avg[1] * Z_Roe[j][0] + Z_avg[0] * Z_Roe[j][1];
          U_hat[j][2] = Z_avg[2] * Z_Roe[j][0] + Z_avg[0] * Z_Roe[j][2];
          U_hat[j][3] = (Z_avg[3] * Z_Roe[j][0] + GAMMA_MINUS1 * Z_avg[1] * Z_Roe[j][1] + GAMMA_MINUS1 * Z_avg[2] * Z_Roe[j][2] +
                         Z_avg[0] * Z_Roe[j][3]) /
                        GAMMA;
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
              Phi[k] += Kmatrix[k][0][j][kfull] * U_hat[j][0] + Kmatrix[k][1][j][kfull] * U_hat[j][1] +
                        Kmatrix[k][2][j][kfull] * U_hat[j][2] + Kmatrix[k][3][j][kfull] * U_hat[j][3];
            }
        }


#endif  // TWO_DIMS
    }  // for loop of triangles

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
}

void apply_fluxRD_list(void)
{
  int i, j, p, nimport, ngrp, recvTask;
#if defined(MAXSCALARS)
  int k;
#endif /* #if defined(MAXSCALARS) */

  /* now exchange the flux-list and apply it when needed */

  mysort(FluxListRD, 1, sizeof(struct fluxRD_list_data), fluxRD_list_data_compare);

  int Nflux = 1;

  for(j = 0; j < NTask; j++)
    Send_count[j] = 0;

  for(i = 0; i < Nflux; i++)
    Send_count[FluxListRD[i].task]++;

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

  struct fluxRD_list_data *FluxListGet = (struct fluxRD_list_data *)mymalloc("FluxListGet", nimport * sizeof(struct fluxRD_list_data));

  /* exchange particle data */
  for(ngrp = 0; ngrp < (1 << PTask); ngrp++)
    {
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
            {
              /* get the particles */
              MPI_Sendrecv(&FluxListRD[Send_offset[recvTask]], Send_count[recvTask] * sizeof(struct fluxRD_list_data), MPI_BYTE,
                           recvTask, TAG_DENS_A, &FluxListGet[Recv_offset[recvTask]],
                           Recv_count[recvTask] * sizeof(struct fluxRD_list_data), MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD,
                           MPI_STATUS_IGNORE);
            }
        }
    }

  /* apply the fluxes */

  for(i = 0; i < nimport; i++)
    {
      p = FluxListGet[i].index;

      SphP[p].Pressure += FluxListGet[i].dPressure;
    }
  myfree(FluxListGet);
}

int fluxRD_list_data_compare(const void *a, const void *b)
{
  if(((struct fluxRD_list_data *)a)->task < (((struct fluxRD_list_data *)b)->task))
    return -1;

  if(((struct fluxRD_list_data *)a)->task > (((struct fluxRD_list_data *)b)->task))
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
