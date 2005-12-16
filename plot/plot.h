
/*
 * Simple diagnostic 2d plot function
 *
 * Copyright 1998 - 2004 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU GENERAL PUBLIC LICENCE :-
 * see the Licence.txt file for licencing details.
 */

/* Plot up to 3 X/Y Graphs. return when the user closes the window */
/* return 0 on success, -1 on error */
int do_plot(double *x, double *y1, double *y2, double *y3, int n);

/* Public routines /
/* Plot up to 6 graphs */
/* return 0 on success, -1 on error */
/* Graph order is Black = Y1, Red = Y2, Green = Y3, Blue = Y4, Yellow = Y5, Purple = Y6 */
int do_plot6(double *x, double *y1, double *y2, double *y3, double *y4, double *y5, double *y6, int n);

/* Plot a bunch of vectors + points */
/* return 0 on success, -1 on error */
int do_plot_vec(double xmin, double xmax, double ymin, double ymax,
                double *x1, double *y1, double *x2, double *y2, int n,
                int dowait, double *x3, double *y3, int m);



