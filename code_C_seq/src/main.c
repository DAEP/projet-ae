#include <stdlib.h>
#include <stdio.h>
#include "problem.h"
#include "visual.h"

int main(int argc, char *argv[])
{
    /*
     * Constants, mesh definition and time discretization
     */
    problem_t *pb = problem_create();

    if(pb == NULL)
    {
        fprintf(stderr, "Could not allocate memory!\n");
        return EXIT_FAILURE;
    }

    problem_set_consts(pb, 20.0e-6);
//    problem_set_mesh(pb, 0.1, 100, 0.1, 1);
    problem_set_mesh(pb, 0.1, 100, 0.1, 100);
    problem_set_tempo(pb, 1., 100);
    problem_set_init_cond(pb, 100.);
    problem_set_bnd(pb, 0, 1, 1000.);
    problem_set_bnd(pb, 1, 1, 1000.);
    problem_set_bnd(pb, 2, 0, 0.);
    problem_set_bnd(pb, 3, 0, 0.);
    problem_check(pb);

    visual_write_case(pb);
    visual_write_geo(pb);

    problem_solve(pb);

    problem_destroy(pb);

    return EXIT_SUCCESS;
}
