#ifndef RAND_H
#define RAND_H

/*
 * Copyright 1998 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU GENERAL PUBLIC LICENCE :-
 * see the Licence.txt file for licencing details.
 */

/* Return a random number 0 and 4294967294 */
unsigned long
rand32(					/* Return 32 bit random number */
unsigned long seed);		/* Optional seed. Non-zero re-initialized with that seed */

/* Return a random integer in the range min to max inclusive */
int i_rand(int min, int max);

/* Return a random double in the range min to max */
double d_rand(double min, double max);

/* Return a random floating point number with a gausian/normal */
/* distribution, centered about 0.0, with standard deviation 1.0 */
double norm_rand(void);

#endif /* RAND_H */
