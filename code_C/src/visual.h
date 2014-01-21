#ifndef VISUAL_H_
#define VISUAL_H_

#include <stdlib.h>
#include <stdio.h>

struct problem_s;
struct parallel_s;

int visual_write_case(struct problem_s *pb, struct parallel_s *par);
int visual_write_geo(struct problem_s *pb, struct parallel_s *par);
int visual_write_solution(struct problem_s *pb, struct parallel_s *par, int it);

#endif /* VISUAL_H_ */
