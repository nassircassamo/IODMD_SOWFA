/* MB02SZ.f -- translated by f2c (version 20041007).
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

/* Subroutine */ int mb02sz(integer *n, doublecomplex *h__, integer *ldh, 
	integer *ipiv, integer *info)
{
    /* System generated locals */
    integer h_dim1, h_offset, i__1, i__2, i__3;
    doublecomplex z__1;

    /* Builtin functions */
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer j, jp;

/*     RELEASE 4.0, WGS COPYRIGHT 1999. */

/*     PURPOSE */

/*     To compute an LU factorization of a complex n-by-n upper */
/*     Hessenberg matrix H using partial pivoting with row interchanges. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrix H.  N >= 0. */

/*     H       (input/output) COMPLEX*16 array, dimension (LDH,N) */
/*             On entry, the n-by-n upper Hessenberg matrix to be */
/*             factored. */
/*             On exit, the factors L and U from the factorization */
/*             H = P*L*U; the unit diagonal elements of L are not stored, */
/*             and L is lower bidiagonal. */

/*     LDH     INTEGER */
/*             The leading dimension of the array H.  LDH >= max(1,N). */

/*     IPIV    (output) INTEGER array, dimension (N) */
/*             The pivot indices; for 1 <= i <= N, row i of the matrix */
/*             was interchanged with row IPIV(i). */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             > 0:  if INFO = i, U(i,i) is exactly zero. The */
/*                   factorization has been completed, but the factor U */
/*                   is exactly singular, and division by zero will occur */
/*                   if it is used to solve a system of equations. */

/*     METHOD */

/*     The factorization has the form */
/*        H = P * L * U */
/*     where P is a permutation matrix, L is lower triangular with unit */
/*     diagonal elements (and one nonzero subdiagonal), and U is upper */
/*     triangular. */

/*     This is the right-looking Level 2 BLAS version of the algorithm */
/*     (adapted after ZGETF2). */

/*     REFERENCES */

/*     - */

/*     NUMERICAL ASPECTS */
/*                                2 */
/*     The algorithm requires 0( N ) complex operations. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Dec. 1996. */
/*     Supersedes Release 2.0 routine TB01FX by A.J. Laub, University of */
/*     Southern California, United States of America, May 1980. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Oct. 2000. */

/*     KEYWORDS */

/*     Frequency response, Hessenberg form, matrix algebra. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Check the scalar input parameters. */

    /* Parameter adjustments */
    h_dim1 = *ldh;
    h_offset = 1 + h_dim1;
    h__ -= h_offset;
    --ipiv;

    /* Function Body */
    *info = 0;
    if (*n < 0) {
	*info = -1;
    } else if (*ldh < max(1,*n)) {
	*info = -3;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla("MB02SZ", &i__1);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0) {
	return 0;
    }

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {

/*        Find pivot and test for singularity. */

	jp = j;
	if (j < *n) {
	    if (dcabs1(&h__[j + 1 + j * h_dim1]) > dcabs1(&h__[j + j * 
		    h_dim1])) {
		jp = j + 1;
	    }
	}
	ipiv[j] = jp;
	i__2 = jp + j * h_dim1;
	if (h__[i__2].r != 0. || h__[i__2].i != 0.) {

/*           Apply the interchange to columns J:N. */

	    if (jp != j) {
		i__2 = *n - j + 1;
		zswap(&i__2, &h__[j + j * h_dim1], ldh, &h__[jp + j * h_dim1]
			, ldh);
	    }

/*           Compute element J+1 of J-th column. */

	    if (j < *n) {
		i__2 = j + 1 + j * h_dim1;
		z_div(&z__1, &h__[j + 1 + j * h_dim1], &h__[j + j * h_dim1]);
		h__[i__2].r = z__1.r, h__[i__2].i = z__1.i;
	    }

	} else if (*info == 0) {

	    *info = j;
	}

	if (j < *n) {

/*           Update trailing submatrix. */

	    i__2 = *n - j;
	    i__3 = j + 1 + j * h_dim1;
	    z__1.r = -h__[i__3].r, z__1.i = -h__[i__3].i;
	    zaxpy(&i__2, &z__1, &h__[j + (j + 1) * h_dim1], ldh, &h__[j + 1 
		    + (j + 1) * h_dim1], ldh);
	}
/* L10: */
    }
    return 0;
/* *** Last line of MB02SZ *** */
} /* mb02sz_ */

