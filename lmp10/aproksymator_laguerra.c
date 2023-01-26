#include "makespl.h"
#include "piv_ge_solver.h"

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>



double laguerr(int i, int alfa, double a)
{
	if (i==0)
		return 1;
	if (i==1)
		return 1 + (double)alfa - a;
	return (((2.0 * i - 1 + alfa - a)  * laguerr(i-1, alfa, a) - (i - 1 + alfa) * laguerr(i - 2, alfa ,a)) / i);
}

double pochodna(int k, int i, int alfa, double a)
{
	int h = k % 2 == 0 ? 1 : -1;
	if (k <= i)
		return h * laguerr(i - k, alfa + k, a);
	return 0;
}

void make_spl(points_t * pts, spline_t * spl, int baza)
{

	matrix_t       *eqs= NULL;
	double         *x = pts->x;
	double         *y = pts->y;
	double		a = x[0];
	double		b = x[pts->n - 1];
	int		i, j, k;
	
  
	eqs = make_matrix(baza, baza + 1);

#ifdef DEBUG
#define TESTBASE 500
	{
		FILE           *tst = fopen("debug_base_plot.txt", "w");
		double		dx = (b - a) / (TESTBASE - 1);
		for (i = 0; i < TESTBASE; i++) {
			fprintf(tst, "%g", a + i * dx);
			for (j = 0; j < baza; j++) {
				fprintf(tst, " %g", laguerr  (j, 0,  a + i * dx));
				fprintf(tst, " %g", pochodna(1, j, 0, a + i * dx));
				fprintf(tst, " %g", pochodna(2, j, 0, a + i * dx));
				fprintf(tst, " %g", pochodna(3, j, 0, a + i * dx));
			}
			fprintf(tst, "\n");
		}
		fclose(tst);
	}
#endif

	for (j = 0; j < baza; j++) {
		for (i = 0; i < baza; i++)
			for (k = 0; k < pts->n; k++)
				add_to_entry_matrix(eqs, j, i, laguerr(i, 0, x[k]) * laguerr( j, 0, x[k]));
		for (k = 0; k < pts->n; k++)
			add_to_entry_matrix(eqs, j, baza, y[k] * laguerr(j, 0, x[k]));
	}

#ifdef DEBUG
	write_matrix(eqs, stdout);
#endif

	if (piv_ge_solver(eqs)) {
		spl->n = 0;
		return;
	}
#ifdef DEBUG 
	write_matrix(eqs, stdout);
#endif

	if (alloc_spl(spl, baza) == 0) {
		double xx;
		double ck;
		for (i = 0; i < spl->n; i++) {
			xx = spl->x[i] = a + i*(b-a)/(spl->n-1);
			xx+= 10.0*DBL_EPSILON;  // zabezpieczenie przed ulokowaniem punktu w poprzednim przedziale
			spl->f[i] = 0;
			spl->f1[i] = 0;
			spl->f2[i] = 0;
			spl->f3[i] = 0;
			for (k = 0; k < baza; k++) {
				ck = get_entry_matrix(eqs, k, baza);
				spl->f[i]  += ck * laguerr (k, 0, xx);
				spl->f1[i] += ck * pochodna (1, k, 0, xx); 
				spl->f2[i] += ck * pochodna (2, k, 0, xx); 
				spl->f3[i] += ck * pochodna (3, k, 0, xx); 
			}
		}
	}

#ifdef DEBUG 
	{
		FILE           *tst = fopen("debug_spline_plot.txt", "w");
		double		dx = (b - a) / (TESTBASE - 1);
		double yi= 0;
		double dyi= 0;
		double d2yi= 0;
		double d3yi= 0;
		double xi;
		for (i = 0; i < TESTBASE; i++) {

			xi= a + i * dx;
			for( k= 0; k < baza; k++ ) {
							yi += get_entry_matrix(eqs, k, baza) *   laguerr (k, 0, xi);
							dyi += get_entry_matrix(eqs, k, baza) *  pochodna (1, k, 0, xi); 
							d2yi += get_entry_matrix(eqs, k, baza) * pochodna (2, k, 0, xi); 
							d3yi += get_entry_matrix(eqs, k, baza) * pochodna (3, k, 0, xi); 
			}
			fprintf(tst, "%g %g %g %g %g\n", xi, yi, dyi, d2yi, d3yi );
		}
		fclose(tst);
	}
#endif

	free_matrix(eqs);

}
