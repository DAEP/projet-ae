#include <stdlib.h>
#include "visual.h"
#include "problem.h"

int visual_write_case(problem_t *pb)
{
    int i = 0;
    FILE *fp = NULL;

    fp = fopen("solut/ensight.case", "w");

    fprintf(fp, "FORMAT\ntype: ensight gold\n\n");
    fprintf(fp, "GEOMETRY\nmodel: ensight.geo\n\n");
    fprintf(fp, "VARIABLE\nscalar per element: Temperature temp_*******.ensight\n\n");
    fprintf(fp, "TIME\ntime set: 1\nnumber of steps: %d\nfilename start number: 0\nfilename increment: 1\ntime values:\n", pb->nb_t);

    for(i = 0 ; i < pb->nb_t ; i++)
    {
        fprintf(fp, "%12.5e\n", (i + 1) * pb->dt);
    }

    fclose(fp);

    return 0;
}

int visual_write_geo(problem_t *pb)
{
    int i = 0;
    int j = 0;
    FILE *fp = NULL;

    fp = fopen("solut/ensight.geo", "w");

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
            fprintf(fp, "%10.5e\n", pb->dx * i);
        }
    }

    for(j = 0 ; j < pb->nb_y + 1 ; j++)
    {
        for(i = 0 ; i < pb->nb_x + 1 ; i++)
        {
            fprintf(fp, "%10.5e\n", pb->dy * j);
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

int visual_write_solution(problem_t *pb, int it)
{
    int i = 0;
    int j = 0;
    FILE *fp = NULL;
    char filename[80];

    sprintf(filename, "solut/temp_%07d.ensight", it);
    fp = fopen(filename, "w");

    fprintf(fp, "Temperature at cells\n");
    fprintf(fp, "part\n1\n");
    fprintf(fp, "block\n");

    for(j = 1 ; j < pb->nb_y + 1 ; j++)
    {
        for(i = 1 ; i < pb->nb_x + 1 ; i++)
        {
            fprintf(fp, "%10.5e\n", pb->temp[i][j]);
        }
    }

    fclose(fp);

    return 0;
}
