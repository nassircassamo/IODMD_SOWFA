/* MB04KD.f -- translated by f2c (version 20041007).
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
static doublereal c_b7 = 1.;
static doublereal c_b14 = 0.;

/* Subroutine */ int mb04kd(char *uplo, integer *n, integer *m, integer *p, 
	doublereal *r__, integer *ldr, doublereal *a, integer *lda, 
	doublereal *b, integer *ldb, doublereal *c__, integer *ldc, 
	doublereal *tau, doublereal *dwork)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, r_dim1, 
	    r_offset, i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static integer i__, im;
    static logical luplo;


/*     RELEASE 4.0, WGS COPYRIGHT 1999. */

/*     PURPOSE */

/*     To calculate a QR factorization of the first block column and */
/*     apply the orthogonal transformations (from the left) also to the */
/*     second block column of a structured matrix, as follows */
/*                          _ */
/*            [ R   0 ]   [ R   C ] */
/*       Q' * [       ] = [       ] */
/*            [ A   B ]   [ 0   D ] */
/*                 _ */
/*     where R and R are upper triangular. The matrix A can be full or */
/*     upper trapezoidal/triangular. The problem structure is exploited. */
/*     This computation is useful, for instance, in combined measurement */
/*     and time update of one iteration of the Kalman filter (square */
/*     root information filter). */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     UPLO    CHARACTER*1 */
/*             Indicates if the matrix A is or not triangular as follows: */
/*             = 'U':  Matrix A is upper trapezoidal/triangular; */
/*             = 'F':  Matrix A is full. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER                 _ */
/*             The order of the matrices R and R.  N >= 0. */

/*     M       (input) INTEGER */
/*             The number of columns of the matrices B, C and D.  M >= 0. */

/*     P       (input) INTEGER */
/*             The number of rows of the matrices A, B and D.  P >= 0. */

/*     R       (input/output) DOUBLE PRECISION array, dimension (LDR,N) */
/*             On entry, the leading N-by-N upper triangular part of this */
/*             array must contain the upper triangular matrix R. */
/*             On exit, the leading N-by-N upper triangular part of this */
/*                                                        _ */
/*             array contains the upper triangular matrix R. */
/*             The strict lower triangular part of this array is not */
/*             referenced. */

/*     LDR     INTEGER */
/*             The leading dimension of array R.  LDR >= MAX(1,N). */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, if UPLO = 'F', the leading P-by-N part of this */
/*             array must contain the matrix A. If UPLO = 'U', the */
/*             leading MIN(P,N)-by-N part of this array must contain the */
/*             upper trapezoidal (upper triangular if P >= N) matrix A, */
/*             and the elements below the diagonal are not referenced. */
/*             On exit, the leading P-by-N part (upper trapezoidal or */
/*             triangular, if UPLO = 'U') of this array contains the */
/*             trailing components (the vectors v, see Method) of the */
/*             elementary reflectors used in the factorization. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,P). */

/*     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M) */
/*             On entry, the leading P-by-M part of this array must */
/*             contain the matrix B. */
/*             On exit, the leading P-by-M part of this array contains */
/*             the computed matrix D. */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= MAX(1,P). */

/*     C       (output) DOUBLE PRECISION array, dimension (LDC,M) */
/*             The leading N-by-M part of this array contains the */
/*             computed matrix C. */

/*     LDC     INTEGER */
/*             The leading dimension of array C.  LDC >= MAX(1,N). */

/*     TAU     (output) DOUBLE PRECISION array, dimension (N) */
/*             The scalar factors of the elementary reflectors used. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (N) */

/*     METHOD */

/*     The routine uses N Householder transformations exploiting the zero */
/*     pattern of the block matrix.  A Householder matrix has the form */

/*                                     ( 1 ), */
/*        H  = I - tau *u *u',    u  = ( v ) */
/*         i          i  i  i      i   (  i) */

/*     where v  is a P-vector, if UPLO = 'F', or an min(i,P)-vector, if */
/*            i */
/*     UPLO = 'U'.  The components of v  are stored in the i-th column */
/*                                     i */
/*     of A, and tau  is stored in TAU(i). */
/*                  i */

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
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    r_dim1 = *ldr;
    r_offset = 1 + r_dim1;
    r__ -= r_offset;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --tau;
    --dwork;

    /* Function Body */
    if (min(*n,*p) == 0) {
	return 0;
    }

    luplo = lsame(uplo, "U");
    im = *p;

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {

/*        Annihilate the I-th column of A and apply the transformations */
/*        to the entire block matrix, exploiting its structure. */

	if (luplo) {
	    im = min(i__,*p);
	}
	i__2 = im + 1;
	dlarfg(&i__2, &r__[i__ + i__ * r_dim1], &a[i__ * a_dim1 + 1], &c__1, 
		&tau[i__]);
	if (tau[i__] != 0.) {

/*                                      [ R(I,I+1:N)        0     ] */
/*           [ w C(I,:) ] := [ 1 v' ] * [                         ] */
/*                                      [ A(1:IM,I+1:N) B(1:IM,:) ] */

	    if (i__ < *n) {
		i__2 = *n - i__;
		dcopy(&i__2, &r__[i__ + (i__ + 1) * r_dim1], ldr, &dwork[1], 
			&c__1);
		i__2 = *n - i__;
		dgemv("Transpose", &im, &i__2, &c_b7, &a[(i__ + 1) * a_dim1 
			+ 1], lda, &a[i__ * a_dim1 + 1], &c__1, &c_b7, &dwork[
			1], &c__1);
	    }
	    dgemv("Transpose", &im, m, &c_b7, &b[b_offset], ldb, &a[i__ * 
		    a_dim1 + 1], &c__1, &c_b14, &c__[i__ + c_dim1], ldc);

/*           [ R(I,I+1:N)      C(I,:)  ]    [ R(I,I+1:N)        0     ] */
/*           [                         ] := [                         ] */
/*           [ A(1:IM,I+1:N) D(1:IM,:) ]    [ A(1:IM,I+1:N) B(1:IM,:) ] */

/*                                                 [ 1 ] */
/*                                         - tau * [   ] * [ w C(I,:) ] */
/*                                                 [ v ] */

	    if (i__ < *n) {
		i__2 = *n - i__;
		d__1 = -tau[i__];
		daxpy(&i__2, &d__1, &dwork[1], &c__1, &r__[i__ + (i__ + 1) * 
			r_dim1], ldr);
		i__2 = *n - i__;
		d__1 = -tau[i__];
		dger(&im, &i__2, &d__1, &a[i__ * a_dim1 + 1], &c__1, &dwork[
			1], &c__1, &a[(i__ + 1) * a_dim1 + 1], lda);
	    }
	    d__1 = -tau[i__];
	    dscal(m, &d__1, &c__[i__ + c_dim1], ldc);
	    dger(&im, m, &c_b7, &a[i__ * a_dim1 + 1], &c__1, &c__[i__ + 
		    c_dim1], ldc, &b[b_offset], ldb);
	}
/* L10: */
    }

    return 0;
/* *** Last line of MB04KD *** */
} /* mb04kd_ */

