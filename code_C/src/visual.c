#include <math.h>
#include <stdlib.h>
#include "visual.h"
#include "parallel.h"
#include "problem.h"

int visual_write_case(problem_t *pb, parallel_t *par)
{
    int i = 0;
    char filename[80];
    FILE *fp = NULL;

    sprintf(filename, "solut/ensight_%04d.case", par->rank);

    fp = fopen(filename, "w");

    fprintf(fp, "FORMAT\ntype: ensight gold\n\n");
    fprintf(fp, "GEOMETRY\nmodel: ensight_%04d.geo\n\n", par->rank);
    fprintf(fp, "VARIABLE\nscalar per element: Temperature temp_*******_%04d.ensight\n\n", par->rank);
    fprintf(fp, "TIME\ntime set: 1\nnumber of steps: %d\nfilename start number: 0\nfilename increment: %d\ntime values:\n", (int)ceil((double)pb->nb_t / (double)pb->write), pb->write);

    for(i = 0 ; i < pb->nb_t ; i += pb->write)
    {
        fprintf(fp, "%12.5e\n", (i + 1) * pb->dt);
    }

    fclose(fp);

    return 0;
}

int visual_write_geo(problem_t *pb, parallel_t *par)
{
    int i = 0;
    int j = 0;
    char filename[80];
    FILE *fp = NULL;

    sprintf(filename, "solut/ensight_%04d.geo", par->rank);

    fp = fopen(filename, "w");

    fprintf(fp, "Ensight Gold Format file\n");
    fprintf(fp, "Heat equation example program\n");
    fprintf(fp, "node id off\nelement id off\n");
    fprintf(fp, "part\n1\n");
    fprintf(fp, "Mesh\nblock\n");
    fprintf(fp, "%d %d 1\n", pb->nb_x + 1, pb->nb_y + 1);


    for(j = 0 ; j < pb->nb_y + 1 ; j++)
    {
        for(i = 0 ; i < pb->nb_x + 1 ; i++)
        {
            fprintf(fp, "%10.5e\n", pb->dx * i + pb->blhc_x);
        }
    }

    for(j = 0 ; j < pb->nb_y + 1 ; j++)
    {
        for(i = 0 ; i < pb->nb_x + 1 ; i++)
        {
            fprintf(fp, "%10.5e\n", pb->dy * j + pb->blhc_y);
        }
    }

    for(j = 0 ; j < pb->nb_y + 1 ; j++)
    {
        for(i = 0 ; i < pb->nb_x + 1 ; i++)
        {
            fprintf(fp, "%10.5e\n", 0.);
        }
    }

    fclose(fp);

    return 0;
}

int visual_write_solution(problem_t *pb, parallel_t *par, int it)
{
    int i = 0;
    int j = 0;
    FILE *fp = NULL;
    char filename[80];

    sprintf(filename, "solut/temp_%07d_%04d.ensight", it, par->rank);
    fp = fopen(filename, "w");

    fprintf(fp, "Temperature at cells\n");
    fprintf(fp, "part\n1\n");
    fprintf(fp, "block\n");

    for(j = 1 ; j < pb->nb_y + 1 ; j++)
    {
        for(i = 1 ; i < pb->nb_x + 1 ; i++)
        {
            fprintf(fp, "%10.5e\n", pb->temp1[i][j]);
        }
    }

    fclose(fp);

    return 0;
}
