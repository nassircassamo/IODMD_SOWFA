/* MB01RU.f -- translated by f2c (version 20041007).
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

static doublereal c_b8 = 0.;
static integer c__0 = 0;
static doublereal c_b12 = 1.;
static doublereal c_b13 = .5;
static doublereal c_b24 = 2.;

/* Subroutine */ int mb01ru(char *uplo, char *trans, integer *m, integer *n, 
	doublereal *alpha, doublereal *beta, doublereal *r__, integer *ldr, 
	doublereal *a, integer *lda, doublereal *x, integer *ldx, doublereal *
	dwork, integer *ldwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, r_dim1, r_offset, x_dim1, x_offset, i__1;

    /* Local variables */
    static logical luplo;
    static logical ltrans;


/*     RELEASE 4.0, WGS COPYRIGHT 1999. */

/*     PURPOSE */

/*     To compute the matrix formula */
/*        _ */
/*        R = alpha*R + beta*op( A )*X*op( A )', */
/*                                                 _ */
/*     where alpha and beta are scalars, R, X, and R are symmetric */
/*     matrices, A is a general matrix, and op( A ) is one of */

/*        op( A ) = A   or   op( A ) = A'. */

/*     The result is overwritten on R. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     UPLO    CHARACTER*1 */
/*             Specifies which triangles of the symmetric matrices R */
/*             and X are given as follows: */
/*             = 'U':  the upper triangular part is given; */
/*             = 'L':  the lower triangular part is given. */

/*     TRANS   CHARACTER*1 */
/*             Specifies the form of op( A ) to be used in the matrix */
/*             multiplication as follows: */
/*             = 'N':  op( A ) = A; */
/*             = 'T':  op( A ) = A'; */
/*             = 'C':  op( A ) = A'. */

/*     Input/Output Parameters */

/*     M       (input) INTEGER           _ */
/*             The order of the matrices R and R and the number of rows */
/*             of the matrix op( A ).  M >= 0. */

/*     N       (input) INTEGER */
/*             The order of the matrix X and the number of columns of the */
/*             the matrix op( A ).  N >= 0. */

/*     ALPHA   (input) DOUBLE PRECISION */
/*             The scalar alpha. When alpha is zero then R need not be */
/*             set before entry, except when R is identified with X in */
/*             the call. */

/*     BETA    (input) DOUBLE PRECISION */
/*             The scalar beta. When beta is zero then A and X are not */
/*             referenced. */

/*     R       (input/output) DOUBLE PRECISION array, dimension (LDR,M) */
/*             On entry with UPLO = 'U', the leading M-by-M upper */
/*             triangular part of this array must contain the upper */
/*             triangular part of the symmetric matrix R. */
/*             On entry with UPLO = 'L', the leading M-by-M lower */
/*             triangular part of this array must contain the lower */
/*             triangular part of the symmetric matrix R. */
/*             On exit, the leading M-by-M upper triangular part (if */
/*             UPLO = 'U'), or lower triangular part (if UPLO = 'L'), of */
/*             this array contains the corresponding triangular part of */
/*                                 _ */
/*             the computed matrix R. */

/*     LDR     INTEGER */
/*             The leading dimension of array R.  LDR >= MAX(1,M). */

/*     A       (input) DOUBLE PRECISION array, dimension (LDA,k) */
/*             where k is N when TRANS = 'N' and is M when TRANS = 'T' or */
/*             TRANS = 'C'. */
/*             On entry with TRANS = 'N', the leading M-by-N part of this */
/*             array must contain the matrix A. */
/*             On entry with TRANS = 'T' or TRANS = 'C', the leading */
/*             N-by-M part of this array must contain the matrix A. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,k), */
/*             where k is M when TRANS = 'N' and is N when TRANS = 'T' or */
/*             TRANS = 'C'. */

/*     X       (input) DOUBLE PRECISION array, dimension (LDX,N) */
/*             On entry, if UPLO = 'U', the leading N-by-N upper */
/*             triangular part of this array must contain the upper */
/*             triangular part of the symmetric matrix X and the strictly */
/*             lower triangular part of the array is not referenced. */
/*             On entry, if UPLO = 'L', the leading N-by-N lower */
/*             triangular part of this array must contain the lower */
/*             triangular part of the symmetric matrix X and the strictly */
/*             upper triangular part of the array is not referenced. */
/*             The diagonal elements of this array are modified */
/*             internally, but are restored on exit. */

/*     LDX     INTEGER */
/*             The leading dimension of array X.  LDX >= MAX(1,N). */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             This array is not referenced when beta = 0, or M*N = 0. */

/*     LDWORK  The length of the array DWORK. */
/*             LDWORK >= M*N, if  beta <> 0; */
/*             LDWORK >= 0,   if  beta =  0. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -k, the k-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     The matrix expression is efficiently evaluated taking the symmetry */
/*     into account. Specifically, let X = T + T', with T an upper or */
/*     lower triangular matrix, defined by */

/*        T = triu( X ) - (1/2)*diag( X ),  if UPLO = 'U', */
/*        T = tril( X ) - (1/2)*diag( X ),  if UPLO = 'L', */

/*     where triu, tril, and diag denote the upper triangular part, lower */
/*     triangular part, and diagonal part of X, respectively. Then, */

/*        A*X*A' = ( A*T )*A' + A*( A*T )',  for TRANS = 'N', */
/*        A'*X*A = A'*( T*A ) + ( T*A )'*A,  for TRANS = 'T', or 'C', */

/*     which involve BLAS 3 operations (DTRMM and DSYR2K). */

/*     NUMERICAL ASPECTS */

/*     The algorithm requires approximately */

/*                   2         2 */
/*        3/2 x M x N + 1/2 x M */

/*     operations. */

/*     FURTHER COMMENTS */

/*     This is a simpler version for MB01RD. */

/*     CONTRIBUTORS */

/*     V. Sima, Katholieke Univ. Leuven, Belgium, Jan. 1999. */

/*     REVISIONS */

/*     A. Varga, German Aerospace Center, Oberpfaffenhofen, March 2004. */
/*     V. Sima, Research Institute for Informatics, Bucharest, Mar. 2004. */

/*     KEYWORDS */

/*     Elementary matrix operations, matrix algebra, matrix operations. */

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
    r_dim1 = *ldr;
    r_offset = 1 + r_dim1;
    r__ -= r_offset;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    x_dim1 = *ldx;
    x_offset = 1 + x_dim1;
    x -= x_offset;
    --dwork;

    /* Function Body */
    *info = 0;
    luplo = lsame(uplo, "U");
    ltrans = lsame(trans, "T") || lsame(trans, "C");

    if (! luplo && ! lsame(uplo, "L")) {
	*info = -1;
    } else if (! ltrans && ! lsame(trans, "N")) {
	*info = -2;
    } else if (*m < 0) {
	*info = -3;
    } else if (*n < 0) {
	*info = -4;
    } else if (*ldr < max(1,*m)) {
	*info = -8;
    } else if (*lda < 1 || ltrans && *lda < *n || ! ltrans && *lda < *m) {
	*info = -10;
    } else if (*ldx < max(1,*n)) {
	*info = -12;
    } else if (*beta != 0. && *ldwork < *m * *n || *beta == 0. && *ldwork < 0)
	     {
	*info = -14;
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla("MB01RU", &i__1);
	return 0;
    }

/*     Quick return if possible. */

    if (*m == 0) {
	return 0;
    }

    if (*beta == 0. || *n == 0) {
	if (*alpha == 0.) {

/*           Special case alpha = 0. */

	    dlaset(uplo, m, m, &c_b8, &c_b8, &r__[r_offset], ldr);
	} else {

/*           Special case beta = 0 or N = 0. */

	    if (*alpha != 1.) {
		dlascl(uplo, &c__0, &c__0, &c_b12, alpha, m, m, &r__[
			r_offset], ldr, info);
	    }
	}
	return 0;
    }

/*     General case: beta <> 0. */
/*     Compute W = op( A )*T or W = T*op( A ) in DWORK, and apply the */
/*     updating formula (see METHOD section). */
/*     Workspace: need M*N. */

    i__1 = *ldx + 1;
    dscal(n, &c_b13, &x[x_offset], &i__1);

    if (ltrans) {

	dlacpy("Full", n, m, &a[a_offset], lda, &dwork[1], n);
	dtrmm("Left", uplo, "NoTranspose", "Non-unit", n, m, &c_b12, &x[
		x_offset], ldx, &dwork[1], n);
	dsyr2k(uplo, trans, m, n, beta, &dwork[1], n, &a[a_offset], lda, 
		alpha, &r__[r_offset], ldr);

    } else {

	dlacpy("Full", m, n, &a[a_offset], lda, &dwork[1], m);
	dtrmm("Right", uplo, "NoTranspose", "Non-unit", m, n, &c_b12, &x[
		x_offset], ldx, &dwork[1], m);
	dsyr2k(uplo, trans, m, n, beta, &dwork[1], m, &a[a_offset], lda, 
		alpha, &r__[r_offset], ldr);

    }

    i__1 = *ldx + 1;
    dscal(n, &c_b24, &x[x_offset], &i__1);

    return 0;
/* *** Last line of MB01RU *** */
} /* mb01ru_ */

