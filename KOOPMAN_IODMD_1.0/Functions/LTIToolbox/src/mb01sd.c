/* MB01SD.f -- translated by f2c (version 20041007).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"
#include "slicot.h"

/* Subroutine */ int mb01sd(char *jobs, integer *m, integer *n, doublereal *
	a, integer *lda, doublereal *r__, doublereal *c__)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j;
    static doublereal cj;


/*     RELEASE 4.0, WGS COPYRIGHT 1999. */

/*     PURPOSE */

/*     To scale a general M-by-N matrix A using the row and column */
/*     scaling factors in the vectors R and C. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     JOBS    CHARACTER*1 */
/*             Specifies the scaling operation to be done, as follows: */
/*             = 'R':  row scaling, i.e., A will be premultiplied */
/*                     by diag(R); */
/*             = 'C':  column scaling, i.e., A will be postmultiplied */
/*                     by diag(C); */
/*             = 'B':  both row and column scaling, i.e., A will be */
/*                     replaced by diag(R) * A * diag(C). */

/*     Input/Output Parameters */

/*     M       (input) INTEGER */
/*             The number of rows of the matrix A.  M >= 0. */

/*     N       (input) INTEGER */
/*             The number of columns of the matrix A.  N >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the M-by-N matrix A. */
/*             On exit, the scaled matrix.  See JOBS for the form of the */
/*             scaled matrix. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= max(1,M). */

/*     R       (input) DOUBLE PRECISION array, dimension (M) */
/*             The row scale factors for A. */
/*             R is not referenced if JOBS = 'C'. */

/*     C       (input) DOUBLE PRECISION array, dimension (N) */
/*             The column scale factors for A. */
/*             C is not referenced if JOBS = 'R'. */


/*     CONTRIBUTOR */

/*     A. Varga, German Aerospace Center, */
/*     DLR Oberpfaffenhofen, April 1998. */
/*     Based on the RASP routine DMSCAL. */

/*    ****************************************************************** */

/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. Executable Statements .. */

/*     Quick return if possible. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --r__;
    --c__;

    /* Function Body */
    if (*m == 0 || *n == 0) {
	return 0;
    }

    if (lsame(jobs, "C")) {

/*        Column scaling, no row scaling. */

	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    cj = c__[j];
	    i__2 = *m;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		a[i__ + j * a_dim1] = cj * a[i__ + j * a_dim1];
/* L10: */
	    }
/* L20: */
	}
    } else if (lsame(jobs, "R")) {

/*        Row scaling, no column scaling. */

	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = *m;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		a[i__ + j * a_dim1] = r__[i__] * a[i__ + j * a_dim1];
/* L30: */
	    }
/* L40: */
	}
    } else if (lsame(jobs, "B")) {

/*        Row and column scaling. */

	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    cj = c__[j];
	    i__2 = *m;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		a[i__ + j * a_dim1] = cj * r__[i__] * a[i__ + j * a_dim1];
/* L50: */
	    }
/* L60: */
	}
    }

    return 0;
/* *** Last line of MB01SD *** */
} /* mb01sd_ */

