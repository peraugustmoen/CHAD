#define cord(r,c) ((r) + (DEPTH)*(c))
#define cord_spec(r,c, D) ((r) + (D)*(c))
#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <float.h>

//#define NEGINF -(FLT_MAX -10.0)
#define NEGINF -100000000.0


void internal_colSum(double *, int, int , double *);

void internal_threshold_matrix(double *, int, int, double, double, int, double );


