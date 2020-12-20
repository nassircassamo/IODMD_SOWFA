/* MB01TD.f -- translated by f2c (version 20041007).
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

/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b10 = 1.;

/* Subroutine */ int mb01td(integer *n, doublereal *a, integer *lda, 
	doublereal *b, integer *ldb, doublereal *dwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, jmin, jmnm;

/*     RELEASE 4.0, WGS COPYRIGHT 1999. */

/*     PURPOSE */

/*     To compute the matrix product A * B, where A and B are upper */
/*     quasi-triangular matrices (that is, block upper triangular with */
/*     1-by-1 or 2-by-2 diagonal blocks) with the same structure. */
/*     The result is returned in the array B. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrices A and B.  N >= 0. */

/*     A       (input) DOUBLE PRECISION array, dimension (LDA,N) */
/*             The leading N-by-N part of this array must contain the */
/*             upper quasi-triangular matrix A. The elements below the */
/*             subdiagonal are not referenced. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= max(1,N). */

/*     B       (input/output) DOUBLE PRECISION array, dimension (LDB,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the upper quasi-triangular matrix B, with the same */
/*             structure as matrix A. */
/*             On exit, the leading N-by-N part of this array contains */
/*             the computed product A * B, with the same structure as */
/*             on entry. */
/*             The elements below the subdiagonal are not referenced. */

/*     LDB     INTEGER */
/*             The leading dimension of the array B.  LDB >= max(1,N). */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (N-1) */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  if the matrices A and B have not the same structure, */
/*                   and/or A and B are not upper quasi-triangular. */

/*     METHOD */

/*     The matrix product A * B is computed column by column, using */
/*     BLAS 2 and BLAS 1 operations. */

/*     FURTHER COMMENTS */

/*     This routine can be used, for instance, for computing powers of */
/*     a real Schur form matrix. */

/*     CONTRIBUTOR */

/*     V. Sima, Katholieke Univ. Leuven, Belgium, June 1998. */

/*     REVISIONS */

/*     V. Sima, Feb. 2000. */

/*     KEYWORDS */

/*     Elementary matrix operations, matrix operations. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */

/*     .. Executable Statements .. */

/*     Test the input scalar arguments. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    --dwork;

    /* Function Body */
    *info = 0;
    if (*n < 0) {
	*info = -1;
    } else if (*lda < max(1,*n)) {
	*info = -3;
    } else if (*ldb < max(1,*n)) {
	*info = -5;
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla("MB01TD", &i__1);
	return 0;
    }

/*     Quick return, if possible. */

    if (*n == 0) {
	return 0;
    } else if (*n == 1) {
	b[b_dim1 + 1] = a[a_dim1 + 1] * b[b_dim1 + 1];
	return 0;
    }

/*     Test the upper quasi-triangular structure of A and B for identity. */

    i__1 = *n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (a[i__ + 1 + i__ * a_dim1] == 0.) {
	    if (b[i__ + 1 + i__ * b_dim1] != 0.) {
		*info = 1;
		return 0;
	    }
	} else if (i__ < *n - 1) {
	    if (a[i__ + 2 + (i__ + 1) * a_dim1] != 0.) {
		*info = 1;
		return 0;
	    }
	}
/* L10: */
    }

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
	i__2 = j + 1;
	jmin = min(i__2,*n);
/* Computing MIN */
	i__2 = jmin, i__3 = *n - 1;
	jmnm = min(i__2,i__3);

/*        Compute the contribution of the subdiagonal of A to the */
/*        j-th column of the product. */

	i__2 = jmnm;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    dwork[i__] = a[i__ + 1 + i__ * a_dim1] * b[i__ + j * b_dim1];
/* L20: */
	}

/*        Multiply the upper triangle of A by the j-th column of B, */
/*        and add to the above result. */

	dtrmv("Upper", "No transpose", "Non-unit", &jmin, &a[a_offset], lda, 
		&b[j * b_dim1 + 1], &c__1);
	daxpy(&jmnm, &c_b10, &dwork[1], &c__1, &b[j * b_dim1 + 2], &c__1);
/* L30: */
    }

    return 0;
/* *** Last line of MB01TD *** */
} /* mb01td_ */

