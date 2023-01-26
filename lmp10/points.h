#ifndef POINTS_H
#define POINTS_H

#include <stdio.h>

typedef struct {
		int n;
		double *x;
		double *y;
} points_t;

void zakoncz_1 ( points_t * pts );

int read_pts_failed ( FILE* inf, points_t *pts);

#endif
