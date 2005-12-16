#ifndef SOBOL_H
#define SOBOL_H

#define SOBOL_MAXBIT 30
#define SOBOL_MAXDIM 40

/* Object definition */
struct _sobol {
	/* Private: */
	int          dim;			/* dimension we're set for */
	unsigned int count;
	double       recipd;
	int          lastq[SOBOL_MAXDIM];
	int          dir[SOBOL_MAXBIT][SOBOL_MAXDIM];

	/* Public: */
	/* Methods */
	int (*next)(struct _sobol *s, double *v);
	void (*reset)(struct _sobol *s);
	void (*del)(struct _sobol *s);

}; typedef struct _sobol sobol;

/* Return NULL on error */
sobol *new_sobol(int dim);

#endif /* SOBOL_H */










