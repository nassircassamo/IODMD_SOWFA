/* MB04ID.f -- translated by f2c (version 20041007).
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

/* Subroutine */ int mb04id(integer *n, integer *m, integer *p, integer *l, 
	doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *
	tau, doublereal *dwork, integer *ldwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__;
    static doublereal first;
    static doublereal wrkopt;


/*     RELEASE 4.0, WGS COPYRIGHT 1999. */

/*     PURPOSE */

/*     To compute a QR factorization of an n-by-m matrix A (A = Q * R), */
/*     having a p-by-min(p,m) zero triangle in the lower left-hand side */
/*     corner, as shown below, for n = 8, m = 7, and p = 2: */

/*            [ x x x x x x x ] */
/*            [ x x x x x x x ] */
/*            [ x x x x x x x ] */
/*            [ x x x x x x x ] */
/*        A = [ x x x x x x x ], */
/*            [ x x x x x x x ] */
/*            [ 0 x x x x x x ] */
/*            [ 0 0 x x x x x ] */

/*     and optionally apply the transformations to an n-by-l matrix B */
/*     (from the left). The problem structure is exploited. This */
/*     computation is useful, for instance, in combined measurement and */
/*     time update of one iteration of the time-invariant Kalman filter */
/*     (square root information filter). */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The number of rows of the matrix A.  N >= 0. */

/*     M       (input) INTEGER */
/*             The number of columns of the matrix A.  M >= 0. */

/*     P       (input) INTEGER */
/*             The order of the zero triagle.  P >= 0. */

/*     L       (input) INTEGER */
/*             The number of columns of the matrix B.  L >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,M) */
/*             On entry, the leading N-by-M part of this array must */
/*             contain the matrix A. The elements corresponding to the */
/*             zero P-by-MIN(P,M) lower trapezoidal/triangular part */
/*             (if P > 0) are not referenced. */
/*             On exit, the elements on and above the diagonal of this */
/*             array contain the MIN(N,M)-by-M upper trapezoidal matrix */
/*             R (R is upper triangular, if N >= M) of the QR */
/*             factorization, and the relevant elements below the */
/*             diagonal contain the trailing components (the vectors v, */
/*             see Method) of the elementary reflectors used in the */
/*             factorization. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,N). */

/*     B       (input/output) DOUBLE PRECISION array, dimension (LDB,L) */
/*             On entry, the leading N-by-L part of this array must */
/*             contain the matrix B. */
/*             On exit, the leading N-by-L part of this array contains */
/*             the updated matrix B. */
/*             If L = 0, this array is not referenced. */

/*     LDB     INTEGER */
/*             The leading dimension of array B. */
/*             LDB >= MAX(1,N) if L > 0; */
/*             LDB >= 1        if L = 0. */

/*     TAU     (output) DOUBLE PRECISION array, dimension MIN(N,M) */
/*             The scalar factors of the elementary reflectors used. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK. */

/*     LDWORK  The length of the array DWORK. */
/*             LDWORK >= MAX(1,M-1,M-P,L). */
/*             For optimum performance LDWORK should be larger. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     The routine uses min(N,M) Householder transformations exploiting */
/*     the zero pattern of the matrix.  A Householder matrix has the form */

/*                                     ( 1 ), */
/*        H  = I - tau *u *u',    u  = ( v ) */
/*         i          i  i  i      i   (  i) */

/*     where v  is an (N-P+I-2)-vector.  The components of v  are stored */
/*            i                                             i */
/*     in the i-th column of A, beginning from the location i+1, and */
/*     tau  is stored in TAU(i). */
/*        i */

/*     NUMERICAL ASPECTS */

/*     The algorithm is backward stable. */

/*     CONTRIBUTORS */

/*     V. Sima, Katholieke Univ. Leuven, Belgium, Feb. 1997. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Elementary reflector, QR factorization, orthogonal transformation. */

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
    --tau;
    --dwork;

    /* Function Body */
    *info = 0;
    if (*n < 0) {
	*info = -1;
    } else if (*m < 0) {
	*info = -2;
    } else if (*p < 0) {
	*info = -3;
    } else if (*l < 0) {
	*info = -4;
    } else if (*lda < max(1,*n)) {
	*info = -6;
    } else if (*l == 0 && *ldb < 1 || *l > 0 && *ldb < max(1,*n)) {
	*info = -8;
    } else /* if(complicated condition) */ {
/* Computing MAX */
	i__1 = 1, i__2 = *m - 1, i__1 = max(i__1,i__2), i__2 = *m - *p, i__1 =
		 max(i__1,i__2);
	if (*ldwork < max(i__1,*l)) {
	    *info = -11;
	}
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla("MB04ID", &i__1);
	return 0;
    }

/*     Quick return if possible. */

    if (min(*m,*n) == 0) {
	dwork[1] = 1.;
	return 0;
    } else if (*n <= *p + 1) {
	i__1 = min(*n,*m);
	for (i__ = 1; i__ <= i__1; ++i__) {
	    tau[i__] = 0.;
/* L5: */
	}
	dwork[1] = 1.;
	return 0;
    }

/*     Annihilate the subdiagonal elements of A and apply the */
/*     transformations to B, if L > 0. */
/*     Workspace: need MAX(M-1,L). */

/*     (Note: Comments in the code beginning "Workspace:" describe the */
/*     minimal amount of real workspace needed at that point in the */
/*     code, as well as the preferred amount for good performance. */
/*     NB refers to the optimal block size for the immediately */
/*     following subroutine, as returned by ILAENV.) */

    i__1 = min(*p,*m);
    for (i__ = 1; i__ <= i__1; ++i__) {

/*        Exploit the structure of the I-th column of A. */

	i__2 = *n - *p;
	dlarfg(&i__2, &a[i__ + i__ * a_dim1], &a[i__ + 1 + i__ * a_dim1], &
		c__1, &tau[i__]);
	if (tau[i__] != 0.) {

	    first = a[i__ + i__ * a_dim1];
	    a[i__ + i__ * a_dim1] = 1.;

	    if (i__ < *m) {
		i__2 = *n - *p;
		i__3 = *m - i__;
		dlarf("Left", &i__2, &i__3, &a[i__ + i__ * a_dim1], &c__1, &
			tau[i__], &a[i__ + (i__ + 1) * a_dim1], lda, &dwork[1]
			);
	    }
	    if (*l > 0) {
		i__2 = *n - *p;
		dlarf("Left", &i__2, l, &a[i__ + i__ * a_dim1], &c__1, &tau[
			i__], &b[i__ + b_dim1], ldb, &dwork[1]);
	    }

	    a[i__ + i__ * a_dim1] = first;
	}
/* L10: */
    }

/* Computing MAX */
    d__1 = 1., d__2 = (doublereal) (*m - 1), d__1 = max(d__1,d__2), d__2 = (
	    doublereal) (*l);
    wrkopt = max(d__1,d__2);

/*     Fast QR factorization of the remaining right submatrix, if any. */
/*     Workspace: need M-P;  prefer (M-P)*NB. */

    if (*m > *p) {
	i__1 = *n - *p;
	i__2 = *m - *p;
	dgeqrf(&i__1, &i__2, &a[*p + 1 + (*p + 1) * a_dim1], lda, &tau[*p + 
		1], &dwork[1], ldwork, info);
	wrkopt = max(wrkopt,dwork[1]);

	if (*l > 0) {

/*           Apply the transformations to B. */
/*           Workspace: need L;  prefer L*NB. */

	    i__1 = *n - *p;
	    i__2 = min(*n,*m) - *p;
	    dormqr("Left", "Transpose", &i__1, l, &i__2, &a[*p + 1 + (*p + 1)
		     * a_dim1], lda, &tau[*p + 1], &b[*p + 1 + b_dim1], ldb, &
		    dwork[1], ldwork, info);
	    wrkopt = max(wrkopt,dwork[1]);
	}
    }

    dwork[1] = wrkopt;
    return 0;
/* *** Last line of MB04ID *** */
} /* mb04id_ */

