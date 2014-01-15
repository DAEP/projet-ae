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
	FILE *fp_config = NULL;

    if(pb == NULL)
    {
        fprintf(stderr, "Could not allocate memory!\n");
        return EXIT_FAILURE;
    }

	fp_config = fopen(argv[1], "r");
	if (fp_config == NULL)
    {
		fprintf(stderr, "Could not open the configuration file!\n");
        return EXIT_FAILURE;
    }

	problem_config(fp_config, pb);

    problem_check(pb);

	fclose(fp_config);

    visual_write_case(pb);
    visual_write_geo(pb);

    problem_solve(pb);

    problem_destroy(pb);

    return EXIT_SUCCESS;
}
