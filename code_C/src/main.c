#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include "parallel.h"
#include "problem.h"
#include "visual.h"

int main(int argc, char *argv[])
{
    int nrq = 0;
    int nc, ntc;
    double t;
    double tmin, tmax, tmean;
    double ti_min, ti_max, ti_mean;
    double tic_min, tic_max, tic_mean;
    problem_t *pb = NULL;
    parallel_t *par = NULL;
    MPI_Request *rq = NULL;

    MPI_Init(&argc, &argv);

    par = parallel_create();
    pb = problem_create();

    if(pb == NULL || par == NULL)
    {
        fprintf(stderr, "Could not allocate memory!\n");
        return EXIT_FAILURE;
    }

    MPI_Comm_size(MPI_COMM_WORLD, &par->np);
    MPI_Comm_rank(MPI_COMM_WORLD, &par->rank);

    parallel_init(par, pb, &nrq, &rq);

    if(par->rank == 0)
    {
        parallel_master(par);
    }

    MPI_Waitall(nrq, rq, MPI_STATUSES_IGNORE);
    free(rq);

    problem_alloc_mesh(pb);
    problem_set_init_cond(pb, pb->t0);

    if(par->rank == 0)
    {
        problem_check(pb);
    }

    if(pb->write != 0)
    {
        visual_write_case(pb, par);
        visual_write_geo(pb, par);
    }

    t = MPI_Wtime();

    problem_solve(pb, par);

    t = MPI_Wtime() - t;

    nc = pb->nb_x * pb->nb_y;

    MPI_Reduce(&t, &tmin, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&t, &tmax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&t, &tmean, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&nc, &ntc, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    tmean = tmean / (double) par->np;

    ti_min = tmin / (double)pb->nb_t;
    ti_max = tmax / (double)pb->nb_t;
    ti_mean = tmean / (double)pb->nb_t;

    tic_min = ti_min / (double)ntc;
    tic_max = ti_max / (double)ntc;
    tic_mean = ti_mean / (double)ntc;

    if(par->rank == 0)
    {
        fprintf(stdout, "\n\n");
        fprintf(stdout, "    ***********************************************************************\n");
        fprintf(stdout, "    * Elapsed time in seconds      min          max          mean         *\n");
        fprintf(stdout, "    *  - Total:                    %-12.4lf %-12.4lf %-12.4lf *\n", tmin, tmax, tmean);
        fprintf(stdout, "    *  - Per iteration:            %-12.4lf %-12.4lf %-12.4lf *\n", ti_min, ti_max, ti_mean);
        fprintf(stdout, "    *  - Per iteration and cell:   %-12.4e %-12.4e %-12.4e *\n", tic_min, tic_max, tic_mean);
        fprintf(stdout, "    ***********************************************************************\n");
        fprintf(stdout, "\n\n");
    }

    problem_destroy(pb);
    parallel_destroy(par);

    MPI_Finalize();
 
    return EXIT_SUCCESS;
}
