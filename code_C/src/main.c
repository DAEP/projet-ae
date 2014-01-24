#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include "parallel.h"
#include "problem.h"
#include "visual.h"

int main(int argc, char *argv[])
{
    problem_t *pb = NULL;
    parallel_t *par = NULL;

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

    if(par->rank == 0)
    {
        parallel_master(par);
    }
    
    parallel_init(par, pb);

    if(par->rank == 0)
    {
        problem_check(pb);
    }

    if(pb->write != 0)
    {
        visual_write_case(pb, par);
        visual_write_geo(pb, par);
    }
   
    problem_solve(pb, par);

    problem_destroy(pb);
    parallel_destroy(par);

    MPI_Finalize();
 
    return EXIT_SUCCESS;
}
