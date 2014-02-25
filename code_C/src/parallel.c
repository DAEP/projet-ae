#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <mpi.h>
#include "parallel.h"
#include "problem.h"

parallel_t *parallel_create(void)
{
    parallel_t *par = calloc(1, sizeof(parallel_t));

    return par;
}

int parallel_master(parallel_t *par)
{
    int nrq = 16;
    int np_yt = 0;
    int neigh[4], bnd_type[4];
    int i = 0;
    int col, row;
    int *nc_x, *nc_y;
    double *cc_x, *cc_y;
    double ar = 0;
    double bnd_value[4];
    problem_t *pbt = NULL;
    FILE *fp = NULL;
    MPI_Request *rq = NULL;

    rq = calloc(nrq, sizeof(*rq));

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
    cc_x = calloc(par->np_x, sizeof(*cc_y));
    cc_y = calloc(par->np_y, sizeof(*cc_y));

    for(i = 0 ; i < par->np_x ; i++)
    {
        nc_x[i] = (int)floor((double)pbt->nb_x / (double)par->np_x);
        nc_x[i] += (pbt->nb_x % par->np_x > i) ? 1 : 0;

        if(i == 0)
        {
            cc_x[0] = 0.0;
        }
        else
        {
            cc_x[i] = cc_x[i - 1] + nc_x[i] * pbt->dx;
        }
    }

    for(i = 0 ; i < par->np_y ; i++)
    {
        nc_y[i] = (int)floor((double)pbt->nb_y / (double)par->np_y);
        nc_y[i] += (pbt->nb_y % par->np_y > i) ? 1 : 0;

        if(i == 0)
        {
            cc_y[0] = 0.0;
        }
        else
        {
            cc_y[i] = cc_y[i - 1] + nc_y[i] * pbt->dy;
        }
    }

    // Send data to the processes
    for(i = 0 ; i < par->np ; i++)
    {
        // Number of partitions along x and y
        MPI_Isend(&par->np_x, 1, MPI_INT, i, 1001, MPI_COMM_WORLD, &rq[0]);
        MPI_Isend(&par->np_y, 1, MPI_INT, i, 1002, MPI_COMM_WORLD, &rq[1]);

        // Send problem definition
        MPI_Isend(&pbt->alpha, 1, MPI_DOUBLE, i, 2001, MPI_COMM_WORLD, &rq[2]);
        MPI_Isend(&pbt->dx, 1, MPI_DOUBLE, i, 2002, MPI_COMM_WORLD, &rq[3]);
        MPI_Isend(&pbt->dy, 1, MPI_DOUBLE, i, 2003, MPI_COMM_WORLD, &rq[4]);
        MPI_Isend(&pbt->dt, 1, MPI_DOUBLE, i, 2004, MPI_COMM_WORLD, &rq[5]);
        MPI_Isend(&pbt->nb_t, 1, MPI_INT, i, 2005, MPI_COMM_WORLD, &rq[6]);
        MPI_Isend(&pbt->t0, 1, MPI_DOUBLE, i, 2006, MPI_COMM_WORLD, &rq[7]);
        MPI_Isend(&pbt->write, 1, MPI_INT, i, 2012, MPI_COMM_WORLD, &rq[8]);

        col = i % par->np_x;
        row = (int)floor((double)i / (double)par->np_x);
        MPI_Isend(&nc_x[col], 1, MPI_INT, i, 2007, MPI_COMM_WORLD, &rq[9]);
        MPI_Isend(&nc_y[row], 1, MPI_INT, i, 2008, MPI_COMM_WORLD, &rq[10]);
        MPI_Isend(&cc_x[col], 1, MPI_DOUBLE, i, 2012, MPI_COMM_WORLD, &rq[14]);
        MPI_Isend(&cc_y[row], 1, MPI_DOUBLE, i, 2013, MPI_COMM_WORLD, &rq[15]);

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

        MPI_Isend(&neigh, 4, MPI_INT, i, 2009, MPI_COMM_WORLD, &rq[11]);
        MPI_Isend(&bnd_type, 4, MPI_INT, i, 2010, MPI_COMM_WORLD, &rq[12]);
        MPI_Isend(&bnd_value, 4, MPI_DOUBLE, i, 2011, MPI_COMM_WORLD, &rq[13]);

        MPI_Waitall(nrq, rq, MPI_STATUSES_IGNORE);
    }

    free(nc_x);
    free(nc_y);
    free(cc_x);
    free(cc_y);

    free(rq);

    problem_destroy(pbt);

    return 0;
}

int parallel_init(parallel_t *par, problem_t *pb, int *nrq, MPI_Request **rq)
{
    int lnrq = 16;
    MPI_Request *lrq = NULL;

    lrq = calloc(lnrq, sizeof(*lrq));

    // Receive partitionning data from master
    MPI_Irecv(&par->np_x, 1, MPI_INT, 0, 1001, MPI_COMM_WORLD, &lrq[0]);
    MPI_Irecv(&par->np_y, 1, MPI_INT, 0, 1002, MPI_COMM_WORLD, &lrq[1]);
    
    // Receive problem definition from master
    MPI_Irecv(&pb->alpha, 1, MPI_DOUBLE, 0, 2001, MPI_COMM_WORLD, &lrq[2]);
    MPI_Irecv(&pb->dx, 1, MPI_DOUBLE, 0, 2002, MPI_COMM_WORLD, &lrq[3]);
    MPI_Irecv(&pb->dy, 1, MPI_DOUBLE, 0, 2003, MPI_COMM_WORLD, &lrq[4]);
    MPI_Irecv(&pb->dt, 1, MPI_DOUBLE, 0, 2004, MPI_COMM_WORLD, &lrq[5]);
    MPI_Irecv(&pb->nb_t, 1, MPI_INT, 0, 2005, MPI_COMM_WORLD, &lrq[6]);
    MPI_Irecv(&pb->t0, 1, MPI_DOUBLE, 0, 2006, MPI_COMM_WORLD, &lrq[7]);
    MPI_Irecv(&pb->write, 1, MPI_INT, 0, 2012, MPI_COMM_WORLD, &lrq[8]);

    MPI_Irecv(&pb->nb_x, 1, MPI_INT, 0, 2007, MPI_COMM_WORLD, &lrq[9]);
    MPI_Irecv(&pb->nb_y, 1, MPI_INT, 0, 2008, MPI_COMM_WORLD, &lrq[10]);
    MPI_Irecv(&pb->blhc_x, 1, MPI_DOUBLE, 0, 2012, MPI_COMM_WORLD, &lrq[14]);
    MPI_Irecv(&pb->blhc_y, 1, MPI_DOUBLE, 0, 2013, MPI_COMM_WORLD, &lrq[15]);

    MPI_Irecv(&par->neigh, 4, MPI_INT, 0, 2009, MPI_COMM_WORLD, &lrq[11]);
    MPI_Irecv(&pb->bnd_type, 4, MPI_INT, 0, 2010, MPI_COMM_WORLD, &lrq[12]);
    MPI_Irecv(&pb->bnd_value, 4, MPI_DOUBLE, 0, 2011, MPI_COMM_WORLD, &lrq[13]);

    *nrq = lnrq;
    *rq = lrq;

    return 0;
}

int parallel_destroy(parallel_t *par)
{
    free(par);

    return 0;
}