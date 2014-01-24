#ifndef PROBLEM_H_
#define PROBLEM_H_

struct parallel_s;

typedef struct problem_s
{
    double alpha;

    double dx;
    double dy;
    double dt;

    int nb_x;
    int nb_y;
    int nb_t;

    int bnd_type[4];
    double bnd_value[4];

    int alloc;
    double t0;

    int write;

    double **temp;
    double **temp_old;
} problem_t;

problem_t *problem_create(void);
int problem_destroy(problem_t *pb);

int problem_config(FILE *fp_config, problem_t *pb, int alloc);

int problem_set_consts(problem_t *pb, double alpha);
int problem_set_write(problem_t *pb, int write);
int problem_set_mesh(problem_t *pb, double size_x, int nb_x, double size_y, int nb_y);
int problem_alloc_mesh(problem_t *pb);
int problem_set_tempo(problem_t *pb, double sim_time, int nb_t);
int problem_set_init_cond(problem_t *pb, double t0);
int problem_set_bnd(problem_t *pb, int bnd_id, int type, double value);
int problem_check(problem_t *pb);
int problem_solve(problem_t *pb, struct parallel_s *par);

#endif /* PROBLEM_H_ */
