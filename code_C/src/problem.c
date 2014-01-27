#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <mpi.h>
#include "parallel.h"
#include "problem.h"
#include "visual.h"

#define FILE_LINE_LENGTH 100

problem_t *problem_create(void)
{
    problem_t *pb = calloc(1, sizeof(problem_t));

    pb->alloc = 0;
    pb->temp1 = NULL;
    pb->temp2 = NULL;

    return pb;
}

int problem_config(FILE *fp_config, problem_t *pb, int alloc)
{
    double diffusivity, length_x, length_y, solution_time, initial_cond, left_bc_value, right_bc_value, bottom_bc_value, top_bc_value;
    int nb_cells_x, nb_cells_y, nb_timestep, left_bc_type, right_bc_type, bottom_bc_type, top_bc_type, write_solut;
    char line[FILE_LINE_LENGTH] = "";

    fgets(line, FILE_LINE_LENGTH, fp_config);
    fgets(line, FILE_LINE_LENGTH, fp_config);
    sscanf(line, "%lf", &diffusivity);
    problem_set_consts(pb, diffusivity);
    
    fgets(line, FILE_LINE_LENGTH, fp_config);
    fgets(line, FILE_LINE_LENGTH, fp_config);
    sscanf(line, "%lf %d %lf %d", &length_x, &nb_cells_x, &length_y, &nb_cells_y);
    problem_set_mesh(pb, length_x, nb_cells_x, length_y, nb_cells_y);

    if(alloc != 0)
    {
        problem_alloc_mesh(pb);
    }

    fgets(line, FILE_LINE_LENGTH, fp_config);
    fgets(line, FILE_LINE_LENGTH, fp_config);
    sscanf(line, "%lf %d", &solution_time, &nb_timestep);
    problem_set_tempo(pb, solution_time, nb_timestep);

    fgets(line, FILE_LINE_LENGTH, fp_config);
    fgets(line, FILE_LINE_LENGTH, fp_config);
    sscanf(line, "%lf", &initial_cond);
    problem_set_init_cond(pb, initial_cond);

    fgets(line, FILE_LINE_LENGTH, fp_config);
    fgets(line, FILE_LINE_LENGTH, fp_config);
    sscanf(line, "%d %lf", &left_bc_type, &left_bc_value);
    problem_set_bnd(pb, 0, left_bc_type, left_bc_value);

    fgets(line, FILE_LINE_LENGTH, fp_config);
    fgets(line, FILE_LINE_LENGTH, fp_config);
    sscanf(line, "%d %lf", &right_bc_type, &right_bc_value);
    problem_set_bnd(pb, 1, right_bc_type, right_bc_value);

    fgets(line, FILE_LINE_LENGTH, fp_config);
    fgets(line, FILE_LINE_LENGTH, fp_config);
    sscanf(line, "%d %lf", &bottom_bc_type, &bottom_bc_value);
    problem_set_bnd(pb, 2, bottom_bc_type, bottom_bc_value);

    fgets(line, FILE_LINE_LENGTH, fp_config);
    fgets(line, FILE_LINE_LENGTH, fp_config);
    sscanf(line, "%d %lf", &top_bc_type, &top_bc_value);
    problem_set_bnd(pb, 3, top_bc_type, top_bc_value);

    fgets(line, FILE_LINE_LENGTH, fp_config);
    fgets(line, FILE_LINE_LENGTH, fp_config);
    sscanf(line, "%d", &write_solut);
    problem_set_write(pb, write_solut);

    return 0;
}

int problem_set_consts(problem_t *pb, double alpha)
{
    pb->alpha = alpha;

    return 0;
}

int problem_set_write(problem_t *pb, int write)
{
    pb->write = write;

    return 0;
}

int problem_set_mesh(problem_t *pb, double size_x, int nb_x, double size_y, int nb_y)
{
    pb->nb_x = nb_x;
    pb->nb_y = nb_y;

    pb->dx = size_x / (double) nb_x;
    pb->dy = size_y / (double) nb_y;

    return 0;
}

int problem_alloc_mesh(problem_t *pb)
{
    int i = 0;
    double *t1, *t2;

    // Create the array containing the solution
    t1 = calloc((pb->nb_x + 2) * (pb->nb_y + 2), sizeof(*t1));
    t2 = calloc((pb->nb_x + 2) * (pb->nb_y + 2), sizeof(*t2));

    if(t1 == NULL || t2 == NULL)
    {
        return -1;
    }

    // Allocate and fill in the arrays that permit access to the data
    // as if we had allocated a 2D array
    pb->temp1 = calloc(pb->nb_x + 2, sizeof(*(pb->temp1)));
    pb->temp2 = calloc(pb->nb_x + 2, sizeof(*(pb->temp2)));

    for(i = 0 ; i < pb->nb_x + 2; i++)
    {
        pb->temp1[i] = &t1[i * (pb->nb_y + 2)];
        pb->temp2[i] = &t2[i * (pb->nb_y + 2)];
    }

    pb->alloc = 1;

    return 0;
}

int problem_set_tempo(problem_t *pb, double sim_time, int nb_t)
{
    pb->nb_t = nb_t;
    pb->dt = sim_time / (double) nb_t;

    return 0;
}

int problem_set_init_cond(problem_t *pb, double t0)
{
    int i = 0, j = 0;

    if(pb->alloc == 1)
    {
        for(i = 0 ; i < pb->nb_x + 2 ; i++)
        {
            for(j = 0 ; j < pb->nb_y + 2 ; j++)
            {
                pb->temp1[i][j] = t0;
            }
        }
    }
    else
    {
        pb->t0 = t0;
    }

    return 0;
}

int problem_set_bnd(problem_t *pb, int bnd_id, int type, double value)
{
    // bnd_id: 0 for left side, 1 for right side, 2 for bottom side and 3 for top side
    // bnd_type: 0 for Neumann (flux), 1 for Dirichlet (physical value)
    pb->bnd_type[bnd_id] = type;
    pb->bnd_value[bnd_id] = value;

    return 0;
}

int problem_check(problem_t *pb)
{
    double delta_min = pb->dx < pb->dy ? pb->dx : pb->dy;
    double cfl = pb->alpha * pb->dt / pow(delta_min, 2);

    // Stability check (CFL condition)
    fprintf(stdout, "Stability check: CFL = %f\n", cfl);

    if(cfl >= 0.25)
    {
        fprintf(stderr, "Scheme won't be stable with the settings you chose. Ensure the CFL number is under 0.25.\n");
        return -1;
    }

    return 0;
}

int problem_solve(problem_t *pb, parallel_t *par)
{
    int i = 0;
    int j = 0;
    int k = 0;
    int p = 0;
    int c = 0;
    int count[4] = {pb->nb_y, pb->nb_y, pb->nb_x, pb->nb_x};
    int pos[4] = {1, pb->nb_x, 1, pb->nb_y};
    int pos_bc[4] = {0, pb->nb_x + 1, 0, pb->nb_y + 1};
    int nrq = 0;
    double flux;
    double temp_ij = 0, c_flux = 0, dx_dy = 0, dy_dx = 0;
    double *exch_in[4] = {NULL, NULL, NULL, NULL};
    double *exch_out[4] = {NULL, NULL, NULL, NULL};
    double **ttmp;
    MPI_Request *rq = NULL;

    // Compute some values that will not change during the run
    dx_dy = pb->dx / pb->dy;
    dy_dx = pb->dy / pb->dx;
    c_flux = (pb->alpha * pb->dt) / (pb->dx * pb->dy);

    // Alloc exchange arrays
    for(p = 0 ; p < 4 ; p++)
    {
        if(par->neigh[p] >= 0)
        {
            exch_in[p] = calloc(count[p], sizeof(*exch_in[p]));
            exch_out[p] = calloc(count[p], sizeof(*exch_out[p]));
            nrq += 2;
        }
        
    }

    // Alloc MPI_Request handles
    if(nrq != 0)
    {
        rq = calloc(nrq, sizeof(*rq));
    }

    // Temporal loop
    for(k = 0 ; k < pb->nb_t ; k++)
    {
        if(par->rank == 0)
        {
            fprintf(stdout, "Iteration %d/%d (%4.2f %%)...\n", k, pb->nb_t, 100. * (double)k / (double)pb->nb_t);
        }

        // Copy temperature data in the out buffer and update ghost cells
        // containing boundary conditions values
        for(p = 0 ; p < 4 ; p++)
        {
            if(par->neigh[p] >= 0)
            {
                if(p < 2)
                {
                    for(c = 0 ; c < count[p] ; c++)
                    {
                        exch_out[p][c] = pb->temp1[pos[p]][c + 1];
                    }
                }
                else
                {
                    for(c = 0 ; c < count[p] ; c++)
                    {
                        exch_out[p][c] = pb->temp1[c + 1][pos[p]];
                    }
                }
            }
        }
        
        // Exchange temperature field between partitions
        nrq = 0;

        for(p = 0 ; p < 4 ; p++)
        {
            if(par->neigh[p] >= 0)
            {
                MPI_Irecv(exch_in[p], count[p], MPI_DOUBLE, par->neigh[p], 3000, MPI_COMM_WORLD, &rq[nrq]);
                MPI_Isend(exch_out[p], count[p], MPI_DOUBLE, par->neigh[p], 3000, MPI_COMM_WORLD, &rq[nrq + 1]);
                nrq += 2;

            }
        }

        MPI_Waitall(nrq, rq, MPI_STATUSES_IGNORE);

        // Update ghost cells with data obtained from other partitions
        // and values given by the boundary conditions
        for(p = 0 ; p < 4 ; p++)
        {
            if(par->neigh[p] >= 0)
            {
                if(p < 2)
                {
                    for(c = 0 ; c < count[p] ; c++)
                    {
                        pb->temp1[pos_bc[p]][c + 1] = exch_in[p][c];
                    }
                }
                else
                {
                    for(c = 0 ; c < count[p] ; c++)
                    {
                        pb->temp1[c + 1][pos_bc[p]] = exch_in[p][c];
                    }
                }
            }
            else
            {
                if(p < 2)
                {
                    for(c = 0 ; c < count[p] ; c++)
                    {
                        pb->temp1[pos_bc[p]][c + 1] = pb->bnd_type[p] * (2 * pb->bnd_value[p] - pb->temp1[pos[p]][c + 1])
                            + (1 - pb->bnd_type[p]) * (pb->temp1[pos[p]][c + 1] - pb->dy * pb->bnd_value[p]);
                    }
                }
                else
                {
                    for(c = 0 ; c < count[p] ; c++)
                    {
                        pb->temp1[c + 1][pos_bc[p]] = pb->bnd_type[p] * (2 * pb->bnd_value[p] - pb->temp1[c + 1][pos[p]]) 
                            + (1 - pb->bnd_type[p]) * (pb->temp1[c + 1][pos[p]] - pb->dx * pb->bnd_value[p]);
                    }
                }
            }
        }

        // Compute fluxes and cell values
        for(i = 1 ; i < pb->nb_x + 1 ; i++)
        {
            for(j = 1 ; j < pb->nb_y + 1 ; j++)
            {
                // Since temp1[i][j] is used 5 times, we copy it to optimize data access
                temp_ij = pb->temp1[i][j];

                // flux = flux[1] + flux[3] - flux[0] - flux[2]
                flux = dy_dx * (pb->temp1[i + 1][j] - temp_ij);
                flux += dx_dy * (pb->temp1[i][j + 1] - temp_ij);
                flux -= dy_dx * (temp_ij - pb->temp1[i - 1][j]);
                flux -= dx_dy * (temp_ij - pb->temp1[i][j - 1]);

                pb->temp2[i][j] = temp_ij + c_flux * flux;
            }
        }

        // Swaps pointers to the solution arrays
        ttmp = pb->temp1;
        pb->temp1 = pb->temp2; 
        pb->temp2 = ttmp;

        if(pb->write != 0)
        {
            visual_write_solution(pb, par, k);
        }
    }

    // Free MPI_Request handles
    if(rq != NULL)
    {
        free(rq);
    }

    // Free exchange arrays
    for(p = 0 ; p < 4 ; p++)
    {
        if(exch_in[p] != NULL)
        {
            free(exch_in[p]);
        }

        if(exch_out[p] != NULL)
        {
            free(exch_out[p]);
        }
    }

    return 0;
}

int problem_destroy(problem_t *pb)
{
    if(pb != NULL)
    {
        if(pb->temp1 != NULL)
        {
            free(*pb->temp1);
            free(pb->temp1);
        }

        if(pb->temp2 != NULL)
        {
            free(*pb->temp2);
            free(pb->temp2);
        }

        free(pb);
    }

    return 0;
}
