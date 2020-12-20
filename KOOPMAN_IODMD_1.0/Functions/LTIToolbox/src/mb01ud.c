/* MB01UD.f -- translated by f2c (version 20041007).
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

static doublereal c_b9 = 0.;
static integer c__1 = 1;

/* Subroutine */ int mb01ud(char *side, char *trans, integer *m, integer *n, 
	doublereal *alpha, doublereal *h__, integer *ldh, doublereal *a, 
	integer *lda, doublereal *b, integer *ldb, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, h_dim1, h_offset, i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static integer i__, j;
    static logical lside;
    static logical ltrans;


/*     RELEASE 4.0, WGS COPYRIGHT 1999. */

/*     PURPOSE */

/*     To compute one of the matrix products */

/*        B = alpha*op( H ) * A, or B = alpha*A * op( H ), */

/*     where alpha is a scalar, A and B are m-by-n matrices, H is an */
/*     upper Hessenberg matrix, and op( H ) is one of */

/*        op( H ) = H   or   op( H ) = H',  the transpose of H. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     SIDE    CHARACTER*1 */
/*             Specifies whether the Hessenberg matrix H appears on the */
/*             left or right in the matrix product as follows: */
/*             = 'L':  B = alpha*op( H ) * A; */
/*             = 'R':  B = alpha*A * op( H ). */

/*     TRANS   CHARACTER*1 */
/*             Specifies the form of op( H ) to be used in the matrix */
/*             multiplication as follows: */
/*             = 'N':  op( H ) = H; */
/*             = 'T':  op( H ) = H'; */
/*             = 'C':  op( H ) = H'. */

/*     Input/Output Parameters */

/*     M       (input) INTEGER */
/*             The number of rows of the matrices A and B.  M >= 0. */

/*     N       (input) INTEGER */
/*             The number of columns of the matrices A and B.  N >= 0. */

/*     ALPHA   (input) DOUBLE PRECISION */
/*             The scalar alpha. When alpha is zero then H is not */
/*             referenced and A need not be set before entry. */

/*     H       (input) DOUBLE PRECISION array, dimension (LDH,k) */
/*             where k is M when SIDE = 'L' and is N when SIDE = 'R'. */
/*             On entry with SIDE = 'L', the leading M-by-M upper */
/*             Hessenberg part of this array must contain the upper */
/*             Hessenberg matrix H. */
/*             On entry with SIDE = 'R', the leading N-by-N upper */
/*             Hessenberg part of this array must contain the upper */
/*             Hessenberg matrix H. */
/*             The elements below the subdiagonal are not referenced, */
/*             except possibly for those in the first column, which */
/*             could be overwritten, but are restored on exit. */

/*     LDH     INTEGER */
/*             The leading dimension of the array H.  LDH >= max(1,k), */
/*             where k is M when SIDE = 'L' and is N when SIDE = 'R'. */

/*     A       (input) DOUBLE PRECISION array, dimension (LDA,N) */
/*             The leading M-by-N part of this array must contain the */
/*             matrix A. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= max(1,M). */

/*     B       (output) DOUBLE PRECISION array, dimension (LDB,N) */
/*             The leading M-by-N part of this array contains the */
/*             computed product. */

/*     LDB     INTEGER */
/*             The leading dimension of the array B.  LDB >= max(1,M). */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     The required matrix product is computed in two steps. In the first */
/*     step, the upper triangle of H is used; in the second step, the */
/*     contribution of the subdiagonal is added. A fast BLAS 3 DTRMM */
/*     operation is used in the first step. */

/*     CONTRIBUTOR */

/*     V. Sima, Katholieke Univ. Leuven, Belgium, January 1999. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Elementary matrix operations, matrix operations. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */

/*     .. Executable Statements .. */

/*     Test the input scalar arguments. */

    /* Parameter adjustments */
    h_dim1 = *ldh;
    h_offset = 1 + h_dim1;
    h__ -= h_offset;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;

    /* Function Body */
    *info = 0;
    lside = lsame(side, "L");
    ltrans = lsame(trans, "T") || lsame(trans, "C");

    if (! lside && ! lsame(side, "R")) {
	*info = -1;
    } else if (! ltrans && ! lsame(trans, "N")) {
	*info = -2;
    } else if (*m < 0) {
	*info = -3;
    } else if (*n < 0) {
	*info = -4;
    } else if (*ldh < 1 || lside && *ldh < *m || ! lside && *ldh < *n) {
	*info = -7;
    } else if (*lda < max(1,*m)) {
	*info = -9;
    } else if (*ldb < max(1,*m)) {
	*info = -11;
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla("MB01UD", &i__1);
	return 0;
    }

/*     Quick return, if possible. */

    if (min(*m,*n) == 0) {
	return 0;
    }

    if (*alpha == 0.) {

/*        Set B to zero and return. */

	dlaset("Full", m, n, &c_b9, &c_b9, &b[b_offset], ldb);
	return 0;
    }

/*     Copy A in B and compute one of the matrix products */
/*       B = alpha*op( triu( H ) ) * A, or */
/*       B = alpha*A * op( triu( H ) ), */
/*     involving the upper triangle of H. */

    dlacpy("Full", m, n, &a[a_offset], lda, &b[b_offset], ldb);
    dtrmm(side, "Upper", trans, "Non-unit", m, n, alpha, &h__[h_offset], ldh,
	     &b[b_offset], ldb);

/*     Add the contribution of the subdiagonal of H. */
/*     If SIDE = 'L', the subdiagonal of H is swapped with the */
/*     corresponding elements in the first column of H, and the */
/*     calculations are organized for column operations. */

    if (lside) {
	if (*m > 2) {
	    i__1 = *m - 2;
	    i__2 = *ldh + 1;
	    dswap(&i__1, &h__[(h_dim1 << 1) + 3], &i__2, &h__[h_dim1 + 3], &
		    c__1);
	}
	if (ltrans) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *m - 1;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    b[i__ + j * b_dim1] += *alpha * h__[i__ + 1 + h_dim1] * a[
			    i__ + 1 + j * a_dim1];
/* L10: */
		}
/* L20: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *m;
		for (i__ = 2; i__ <= i__2; ++i__) {
		    b[i__ + j * b_dim1] += *alpha * h__[i__ + h_dim1] * a[i__ 
			    - 1 + j * a_dim1];
/* L30: */
		}
/* L40: */
	    }
	}
	if (*m > 2) {
	    i__1 = *m - 2;
	    i__2 = *ldh + 1;
	    dswap(&i__1, &h__[(h_dim1 << 1) + 3], &i__2, &h__[h_dim1 + 3], &
		    c__1);
	}

    } else {

	if (ltrans) {
	    i__1 = *n - 1;
	    for (j = 1; j <= i__1; ++j) {
		if (h__[j + 1 + j * h_dim1] != 0.) {
		    d__1 = *alpha * h__[j + 1 + j * h_dim1];
		    daxpy(m, &d__1, &a[j * a_dim1 + 1], &c__1, &b[(j + 1) * 
			    b_dim1 + 1], &c__1);
		}
/* L50: */
	    }
	} else {
	    i__1 = *n - 1;
	    for (j = 1; j <= i__1; ++j) {
		if (h__[j + 1 + j * h_dim1] != 0.) {
		    d__1 = *alpha * h__[j + 1 + j * h_dim1];
		    daxpy(m, &d__1, &a[(j + 1) * a_dim1 + 1], &c__1, &b[j * 
			    b_dim1 + 1], &c__1);
		}
/* L60: */
	    }
	}
    }

    return 0;
/* *** Last line of MB01UD *** */
} /* mb01ud_ */

