#ifndef PARALLEL_H_
#define PARALLEL_H_

#include <mpi.h>

struct problem_s;

typedef struct parallel_s
{
    int rank;
    int np;
    int np_x;
    int np_y;
    int neigh[4];
} parallel_t;

parallel_t *parallel_create(void);
int parallel_destroy(parallel_t *par);

int parallel_master(parallel_t *par);
int parallel_init(parallel_t *par, struct problem_s *pb, int *nrq, MPI_Request **rq);

#endif /* PARALLEL_H_ */
