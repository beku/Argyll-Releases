
#ifndef PLOT_H

/*
 * Simple diagnostic 2d plot function
 *
 * Copyright 1998 - 2004 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

/* Graph order is Black = Y1, Red = Y2, Green = Y3, Blue = Y4, Yellow = Y5, Purple = Y6 */

/* Plot up to 3 X/Y Graphs. return when the user closes the window */
/* return 0 on success, -1 on error */
int do_plot(double *x, double *y1, double *y2, double *y3, int n);

/* Plot up to 3 graphs. */
/* if dowait > 0, wait for user key */
/* if dowait < 0, wait for no seconds */
/* If xmax > xmin, use as x scale, else auto. */
/* If ymax > ymin, use as y scale, else auto. */
/* ratio is window X / Y */
/* return 0 on success, -1 on error */
/* If n is -ve, reverse the X axis */
int do_plot_x(double *x, double *y1, double *y2, double *y3, int n,
int dowait, double pxmin, double pxmax, double pymin, double pymax, double ratio);

/* Public routines */
/* Plot up to 6 graphs */
/* return 0 on success, -1 on error */
/* Graph order is Black = Y1, Red = Y2, Green = Y3, Blue = Y4, Yellow = Y5, Purple = Y6 */
int do_plot6(double *x, double *y1, double *y2, double *y3, double *y4, double *y5, double *y6, int n);

/* Plot a bunch of vectors + points */
/* return 0 on success, -1 on error */
int do_plot_vec(double xmin, double xmax, double ymin, double ymax,
                double *x1, double *y1, double *x2, double *y2, int n,
                int dowait, double *x3, double *y3, int m);

#define PLOT_H
#endif /* PLOT_H */
