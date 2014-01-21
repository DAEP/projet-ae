#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <mpi/mpi.h>
#include "parallel.h"
#include "problem.h"

parallel_t *parallel_create(void)
{
    parallel_t *par = calloc(1, sizeof(parallel_t));

    return par;
}

int parallel_master(parallel_t *par)
{
    int np_yt = 0;
    int neigh[4], bnd_type[4];
    int i = 0;
    int col, row;
    int *nc_x, *nc_y;
    double ar = 0;
    double bnd_value[4];
    problem_t *pbt = NULL;
    FILE *fp = NULL;

    // Read the configuration file and create a temporary
    // problem data holder
    fp = fopen("config.dat", "r");

    if(fp == NULL)
    {
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    pbt = problem_create();
    problem_config(fp, pbt, 0);
    fclose(fp);

    // Mesh subdivision
    ar = (double)pbt->nb_x / (double)pbt->nb_y; // Mesh aspect ratio

    np_yt = (int)ceil(sqrt((double)par->np / ar));

    while(par->np % np_yt != 0)
    {
        np_yt--;
    }

    par->np_y = np_yt;
    par->np_x = par->np / np_yt;

    nc_x = calloc(par->np_x, sizeof(*nc_x));
    nc_y = calloc(par->np_y, sizeof(*nc_y));

    for(i = 0 ; i < par->np_x ; i++)
    {
        nc_x[i] = (int)floor((double)pbt->nb_x / (double)par->np_x);
        nc_x[i] += (pbt->nb_x % par->np_x > i) ? 1 : 0;
    }

    for(i = 0 ; i < par->np_y ; i++)
    {
        nc_y[i] = (int)floor((double)pbt->nb_y / (double)par->np_y);
        nc_y[i] += (pbt->nb_y % par->np_y > i) ? 1 : 0;
    }

    // Send data to the processes
    for(i = 0 ; i < par->np ; i++)
    {
        // Number of partitions along x and y
        MPI_Send(&par->np_x, 1, MPI_INT, i, 1001, MPI_COMM_WORLD);
        MPI_Send(&par->np_y, 1, MPI_INT, i, 1002, MPI_COMM_WORLD);

        // Send problem definition
        MPI_Send(&pbt->alpha, 1, MPI_DOUBLE, i, 2001, MPI_COMM_WORLD);
        MPI_Send(&pbt->dx, 1, MPI_DOUBLE, i, 2002, MPI_COMM_WORLD);
        MPI_Send(&pbt->dy, 1, MPI_DOUBLE, i, 2003, MPI_COMM_WORLD);
        MPI_Send(&pbt->dt, 1, MPI_DOUBLE, i, 2004, MPI_COMM_WORLD);
        MPI_Send(&pbt->nb_t, 1, MPI_INT, i, 2005, MPI_COMM_WORLD);
        MPI_Send(&pbt->t0, 1, MPI_DOUBLE, i, 2006, MPI_COMM_WORLD);

        col = i % par->np_x;
        row = (int)floor((double)i / (double)par->np_x);
        MPI_Send(&nc_x[col], 1, MPI_INT, i, 2007, MPI_COMM_WORLD);
        MPI_Send(&nc_y[row], 1, MPI_INT, i, 2008, MPI_COMM_WORLD);

        memset(bnd_value, 0., sizeof(bnd_value));
        memset(bnd_type, -1, sizeof(bnd_type));

        neigh[0] = i - 1;
        neigh[1] = i + 1;
        neigh[2] = i - par->np_x;
        neigh[3] = i + par->np_x;

        if(col == 0)
        {
            neigh[0] = -1;
            bnd_type[0] = pbt->bnd_type[0];
            bnd_value[0] = pbt->bnd_value[0];
        }

        if(col == par->np_x - 1)
        {
            neigh[1] = -1;
            bnd_type[1] = pbt->bnd_type[1];
            bnd_value[1] = pbt->bnd_value[1];
        }

        if(row == 0)
        {
            neigh[2] = -1;
            bnd_type[2] = pbt->bnd_type[2];
            bnd_value[2] = pbt->bnd_value[2];
        }

        if(row == par->np_y - 1)
        {
            neigh[3] = -1;
            bnd_type[3] = pbt->bnd_type[3];
            bnd_value[3] = pbt->bnd_value[3];
        }

        MPI_Send(&neigh, 4, MPI_INT, i, 2009, MPI_COMM_WORLD);
        MPI_Send(&bnd_type, 4, MPI_INT, i, 2010, MPI_COMM_WORLD);
        MPI_Send(&bnd_value, 4, MPI_DOUBLE, i, 2011, MPI_COMM_WORLD);

    }

    free(nc_x);
    free(nc_y);

    problem_destroy(pbt);

    return 0;
}

int parallel_init(parallel_t *par, problem_t *pb)
{
    // Receive partitionning data from master
    MPI_Recv(&par->np_x, 1, MPI_INT, 0, 1001, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&par->np_y, 1, MPI_INT, 0, 1002, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    
    // Receive problem definition from master
    MPI_Recv(&pb->alpha, 1, MPI_DOUBLE, 0, 2001, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&pb->dx, 1, MPI_DOUBLE, 0, 2002, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&pb->dy, 1, MPI_DOUBLE, 0, 2003, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&pb->dt, 1, MPI_DOUBLE, 0, 2004, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&pb->nb_t, 1, MPI_INT, 0, 2005, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&pb->t0, 1, MPI_DOUBLE, 0, 2006, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    MPI_Recv(&pb->nb_x, 1, MPI_INT, 0, 2007, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&pb->nb_y, 1, MPI_INT, 0, 2008, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    MPI_Recv(&par->neigh, 4, MPI_INT, 0, 2009, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&pb->bnd_type, 4, MPI_INT, 0, 2010, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&pb->bnd_value, 4, MPI_DOUBLE, 0, 2011, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    problem_alloc_mesh(pb);
    problem_set_init_cond(pb, pb->t0);

    return 0;
}

int parallel_destroy(parallel_t *par)
{
    free(par);

    return 0;
}
