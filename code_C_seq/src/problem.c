#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "problem.h"
#include "visual.h"

#define FILE_LINE_LENGTH 100

problem_t *problem_create(void)
{
    problem_t *pb = calloc(1, sizeof(problem_t));

    return pb;
}

int problem_config(FILE *fp_config, problem_t *pb)
{
	double diffusivity, length_x, length_y, solution_time, initial_cond, left_bc_value, right_bc_value, bottom_bc_value, top_bc_value;
	int nb_cells_x, nb_cells_y, nb_timestep, left_bc_type, right_bc_type, bottom_bc_type, top_bc_type;
	char line[FILE_LINE_LENGTH] = "";

	fgets(line, FILE_LINE_LENGTH, fp_config);
	fgets(line, FILE_LINE_LENGTH, fp_config);
	sscanf(line, "%lf", &diffusivity);
	problem_set_consts(pb, diffusivity);
	
	fgets(line, FILE_LINE_LENGTH, fp_config);
 	fgets(line, FILE_LINE_LENGTH, fp_config);
	sscanf(line, "%lf %d %lf %d", &length_x, &nb_cells_x, &length_y, &nb_cells_y);
    problem_set_mesh(pb, length_x, nb_cells_x, length_y, nb_cells_y);

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

	return 0;
}

int problem_set_consts(problem_t *pb, double alpha)
{
    pb->alpha = alpha;

    return 0;
}

int problem_set_mesh(problem_t *pb, double size_x, int nb_x, double size_y, int nb_y)
{
    int i = 0;

    pb->nb_x = nb_x;
    pb->nb_y = nb_y;

    pb->dx = size_x / (double) nb_x;
    pb->dy = size_y / (double) nb_y;

    // Create the array containing the solution
    pb->temp = calloc(pb->nb_x, sizeof(*(pb->temp)));
    pb->temp_old = calloc(pb->nb_x + 2, sizeof(*(pb->temp_old)));

    if(pb->temp == NULL || pb->temp_old == NULL)
    {
        return -1;
    }

    for(i = 0 ; i < pb->nb_x + 2; i++)
    {
        // pb->temp does not contain ghost cells
        if(i < pb->nb_x)
        {
            pb->temp[i] = calloc(pb->nb_y, sizeof(**(pb->temp)));
            assert(pb->temp[i] != NULL);
        }

        pb->temp_old[i] = calloc(pb->nb_y + 2, sizeof(**(pb->temp_old)));

        assert(pb->temp_old[i] != NULL);
    }

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

    for(i = 0 ; i < pb->nb_x + 2 ; i++)
    {
        for(j = 0 ; j < pb->nb_y + 2 ; j++)
        {
            pb->temp_old[i][j] = t0;
        }
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

int problem_solve(problem_t *pb)
{
    int i = 0;
    int j = 0;
    int k = 0;
    double flux;

    // Temporal loop
    for(k = 0 ; k < pb->nb_t ; k++)
    {
        fprintf(stdout, "Iteration %d...\n", k);

        // Update ghost cells on the left and on the right
        for(j = 1 ; j < pb->nb_y + 1 ; j++)
        {
            pb->temp_old[0][j] = pb->bnd_type[0] * (2 * pb->bnd_value[0] - pb->temp_old[1][j]) + (1 - pb->bnd_type[0]) * (pb->temp_old[1][j] - pb->dy * pb->bnd_value[0]);
            pb->temp_old[pb->nb_x + 1][j] = pb->bnd_type[1] * (2 * pb->bnd_value[1] - pb->temp_old[pb->nb_x][j]) + (1 - pb->bnd_type[1]) * (pb->temp_old[pb->nb_x][j] - pb->dy * pb->bnd_value[1]);
        }

        // Update ghost cells at the bottom and at the top
        for(i = 1 ; i < pb->nb_x + 1 ; i++)
        {
            pb->temp_old[i][0] = pb->bnd_type[2] * (2 * pb->bnd_value[2] - pb->temp_old[i][1]) + (1 - pb->bnd_type[2]) * (pb->temp_old[i][1] - pb->dx * pb->bnd_value[2]);
            pb->temp_old[i][pb->nb_y + 1] = pb->bnd_type[3] * (2 * pb->bnd_value[3] - pb->temp_old[i][pb->nb_y]) + (1 - pb->bnd_type[3]) * (pb->temp_old[i][pb->nb_y] - pb->dx * pb->bnd_value[3]);
        }

        // Compute fluxes and cell values
        for(i = 0 ; i < pb->nb_x ; i++)
        {
            for(j = 0 ; j < pb->nb_y ; j++)
            {
                // flux = flux[1] + flux[3] - flux[0] - flux[2]
                flux = pb->dy / pb->dx * (pb->temp_old[i + 2][j + 1] - pb->temp_old[i + 1][j + 1]);
                flux += pb->dx / pb->dy * (pb->temp_old[i + 1][j + 2] - pb->temp_old[i + 1][j + 1]);
                flux -= pb->dy / pb->dx * (pb->temp_old[i + 1][j + 1] - pb->temp_old[i][j + 1]);
                flux -= pb->dx / pb->dy * (pb->temp_old[i + 1][j + 1] - pb->temp_old[i + 1][j]);

                pb->temp[i][j] = pb->temp_old[i + 1][j + 1] + (pb->alpha * pb->dt) / (pb->dx * pb->dy) * flux;
            }
        }

        for(i = 0 ; i < pb->nb_x ; i++)
        {
            memcpy(&(pb->temp_old[i + 1][1]), &(pb->temp[i][0]), pb->nb_y * sizeof(**(pb->temp)));
        }

        visual_write_solution(pb, k);
    }

    return 0;
}

int problem_destroy(problem_t *pb)
{
    int i = 0;

    if(pb != NULL)
    {
        if(pb->temp != NULL)
        {
            for(i = 0 ; i < pb->nb_x ; i++)
            {
                if(pb->temp[i] != NULL)
                {
                    free(pb->temp[i]);
                }
            }

            free(pb->temp);
        }

        if(pb->temp_old != NULL)
        {
            for(i = 0 ; i < pb->nb_x + 2 ; i++)
            {
                if(pb->temp_old[i] != NULL)
                {
                    free(pb->temp_old[i]);
                }
            }

            free(pb->temp_old);
        }

        free(pb);
    }

    return 0;
}
