/* SB02MT.f -- translated by f2c (version 20041007).
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
static doublereal c_b28 = 1.;
static doublereal c_b31 = 0.;
static doublereal c_b37 = -1.;

/* Subroutine */ int sb02mt(char *jobg, char *jobl, char *fact, char *uplo, 
	integer *n, integer *m, doublereal *a, integer *lda, doublereal *b, 
	integer *ldb, doublereal *q, integer *ldq, doublereal *r__, integer *
	ldr, doublereal *l, integer *ldl, integer *ipiv, integer *oufact, 
	doublereal *g, integer *ldg, integer *iwork, doublereal *dwork, 
	integer *ldwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, g_dim1, g_offset, l_dim1, 
	    l_offset, q_dim1, q_offset, r_dim1, r_offset, i__1, i__2, i__3, 
	    i__4;

    /* Local variables */
    static integer i__, j;
    static doublereal eps;
    static logical ljobg;
    static logical ljobl;
    static doublereal rcond;
    static char trans[1];
    static doublereal rnorm;
    static logical lfacta, lfactc;
    static logical lfactu;
    static logical luplou;
    static integer wrkopt;


/*     RELEASE 4.0, WGS COPYRIGHT 1999. */

/*     PURPOSE */

/*     To compute the following matrices */

/*                -1 */
/*         G = B*R  *B', */

/*         -          -1 */
/*         A = A - B*R  *L', */

/*         -          -1 */
/*         Q = Q - L*R  *L', */

/*     where A, B, Q, R, L, and G are N-by-N, N-by-M, N-by-N, M-by-M, */
/*     N-by-M, and N-by-N matrices, respectively, with Q, R and G */
/*     symmetric matrices. */

/*     When R is well-conditioned with respect to inversion, standard */
/*     algorithms for solving linear-quadratic optimization problems will */
/*     then also solve optimization problems with coupling weighting */
/*     matrix L. Moreover, a gain in efficiency is possible using matrix */
/*     G in the deflating subspace algorithms (see SLICOT Library routine */
/*     SB02OD). */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     JOBG    CHARACTER*1 */
/*             Specifies whether or not the matrix G is to be computed, */
/*             as follows: */
/*             = 'G':  Compute G; */
/*             = 'N':  Do not compute G. */

/*     JOBL    CHARACTER*1 */
/*             Specifies whether or not the matrix L is zero, as follows: */
/*             = 'Z':  L is zero; */
/*             = 'N':  L is nonzero. */

/*     FACT    CHARACTER*1 */
/*             Specifies how the matrix R is given (factored or not), as */
/*             follows: */
/*             = 'N':  Array R contains the matrix R; */
/*             = 'C':  Array R contains the Cholesky factor of R; */
/*             = 'U':  Array R contains the symmetric indefinite UdU' or */
/*                     LdL' factorization of R. */

/*     UPLO    CHARACTER*1 */
/*             Specifies which triangle of the matrices R and Q (if */
/*             JOBL = 'N') is stored, as follows: */
/*             = 'U':  Upper triangle is stored; */
/*             = 'L':  Lower triangle is stored. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrices A, Q, and G, and the number of */
/*             rows of the matrices B and L.  N >= 0. */

/*     M       (input) INTEGER */
/*             The order of the matrix R, and the number of columns of */
/*             the matrices B and L.  M >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, if JOBL = 'N', the leading N-by-N part of this */
/*             array must contain the matrix A. */
/*             On exit, if JOBL = 'N', and INFO = 0, the leading N-by-N */
/*                                                    -          -1 */
/*             part of this array contains the matrix A = A - B*R  L'. */
/*             If JOBL = 'Z', this array is not referenced. */

/*     LDA     INTEGER */
/*             The leading dimension of array A. */
/*             LDA >= MAX(1,N) if JOBL = 'N'; */
/*             LDA >= 1        if JOBL = 'Z'. */

/*     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M) */
/*             On entry, the leading N-by-M part of this array must */
/*             contain the matrix B. */
/*             On exit, if OUFACT = 1, and INFO = 0, the leading N-by-M */
/*                                                             -1 */
/*             part of this array contains the matrix B*chol(R)  . */
/*             On exit, B is unchanged if OUFACT = 2 (hence also when */
/*             FACT = 'U'). */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= MAX(1,N). */

/*     Q       (input/output) DOUBLE PRECISION array, dimension (LDQ,N) */
/*             On entry, if JOBL = 'N', the leading N-by-N upper */
/*             triangular part (if UPLO = 'U') or lower triangular part */
/*             (if UPLO = 'L') of this array must contain the upper */
/*             triangular part or lower triangular part, respectively, of */
/*             the symmetric matrix Q. The stricly lower triangular part */
/*             (if UPLO = 'U') or stricly upper triangular part (if */
/*             UPLO = 'L') is not referenced. */
/*             On exit, if JOBL = 'N' and INFO = 0, the leading N-by-N */
/*             upper triangular part (if UPLO = 'U') or lower triangular */
/*             part (if UPLO = 'L') of this array contains the upper */
/*             triangular part or lower triangular part, respectively, of */
/*                                  -          -1 */
/*             the symmetric matrix Q = Q - L*R  *L'. */
/*             If JOBL = 'Z', this array is not referenced. */

/*     LDQ     INTEGER */
/*             The leading dimension of array Q. */
/*             LDQ >= MAX(1,N) if JOBL = 'N'; */
/*             LDQ >= 1        if JOBL = 'Z'. */

/*     R       (input/output) DOUBLE PRECISION array, dimension (LDR,M) */
/*             On entry, if FACT = 'N', the leading M-by-M upper */
/*             triangular part (if UPLO = 'U') or lower triangular part */
/*             (if UPLO = 'L') of this array must contain the upper */
/*             triangular part or lower triangular part, respectively, */
/*             of the symmetric input weighting matrix R. */
/*             On entry, if FACT = 'C', the leading M-by-M upper */
/*             triangular part (if UPLO = 'U') or lower triangular part */
/*             (if UPLO = 'L') of this array must contain the Cholesky */
/*             factor of the positive definite input weighting matrix R */
/*             (as produced by LAPACK routine DPOTRF). */
/*             On entry, if FACT = 'U', the leading M-by-M upper */
/*             triangular part (if UPLO = 'U') or lower triangular part */
/*             (if UPLO = 'L') of this array must contain the factors of */
/*             the UdU' or LdL' factorization, respectively, of the */
/*             symmetric indefinite input weighting matrix R (as produced */
/*             by LAPACK routine DSYTRF). */
/*             If FACT = 'N', the stricly lower triangular part (if UPLO */
/*             = 'U') or stricly upper triangular part (if UPLO = 'L') of */
/*             this array is used as workspace. */
/*             On exit, if OUFACT = 1, and INFO = 0 (or INFO = M+1), */
/*             the leading M-by-M upper triangular part (if UPLO = 'U') */
/*             or lower triangular part (if UPLO = 'L') of this array */
/*             contains the Cholesky factor of the given input weighting */
/*             matrix. */
/*             On exit, if OUFACT = 2, and INFO = 0 (or INFO = M+1), */
/*             the leading M-by-M upper triangular part (if UPLO = 'U') */
/*             or lower triangular part (if UPLO = 'L') of this array */
/*             contains the factors of the UdU' or LdL' factorization, */
/*             respectively, of the given input weighting matrix. */
/*             On exit R is unchanged if FACT = 'C' or 'U'. */

/*     LDR     INTEGER */
/*             The leading dimension of array R.  LDR >= MAX(1,M). */

/*     L       (input/output) DOUBLE PRECISION array, dimension (LDL,M) */
/*             On entry, if JOBL = 'N', the leading N-by-M part of this */
/*             array must contain the matrix L. */
/*             On exit, if JOBL = 'N', OUFACT = 1, and INFO = 0, the */
/*             leading N-by-M part of this array contains the matrix */
/*                      -1 */
/*             L*chol(R)  . */
/*             On exit, L is unchanged if OUFACT = 2 (hence also when */
/*             FACT = 'U'). */
/*             L is not referenced if JOBL = 'Z'. */

/*     LDL     INTEGER */
/*             The leading dimension of array L. */
/*             LDL >= MAX(1,N) if JOBL = 'N'; */
/*             LDL >= 1        if JOBL = 'Z'. */

/*     IPIV    (input/output) INTEGER array, dimension (M) */
/*             On entry, if FACT = 'U', this array must contain details */
/*             of the interchanges performed and the block structure of */
/*             the d factor in the UdU' or LdL' factorization of matrix R */
/*             (as produced by LAPACK routine DSYTRF). */
/*             On exit, if OUFACT = 2, this array contains details of */
/*             the interchanges performed and the block structure of the */
/*             d factor in the UdU' or LdL' factorization of matrix R, */
/*             as produced by LAPACK routine DSYTRF. */
/*             This array is not referenced if FACT = 'C'. */

/*     OUFACT  (output) INTEGER */
/*             Information about the factorization finally used. */
/*             OUFACT = 1:  Cholesky factorization of R has been used; */
/*             OUFACT = 2:  UdU' (if UPLO = 'U') or LdL' (if UPLO = 'L') */
/*                          factorization of R has been used. */

/*     G       (output) DOUBLE PRECISION array, dimension (LDG,N) */
/*             If JOBG = 'G', and INFO = 0, the leading N-by-N upper */
/*             triangular part (if UPLO = 'U') or lower triangular part */
/*             (if UPLO = 'L') of this array contains the upper */
/*             triangular part (if UPLO = 'U') or lower triangular part */
/*                                                                 -1 */
/*             (if UPLO = 'L'), respectively, of the matrix G = B*R  B'. */
/*             If JOBG = 'N', this array is not referenced. */

/*     LDG     INTEGER */
/*             The leading dimension of array G. */
/*             LDG >= MAX(1,N) if JOBG = 'G', */
/*             LDG >= 1        if JOBG = 'N'. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (M) */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK; if FACT = 'N', DWORK(2) contains the reciprocal */
/*             condition number of the given matrix R. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= 1              if FACT = 'C'; */
/*             LDWORK >= MAX(2,3*M,N*M) if FACT = 'N'; */
/*             LDWORK >= MAX(1,N*M)     if FACT = 'U'. */
/*             For optimum performance LDWORK should be larger than 3*M, */
/*             if FACT = 'N'. */
/*             The N*M workspace is not needed for FACT = 'N', if matrix */
/*             R is positive definite. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = i:  if the i-th element (1 <= i <= M) of the d factor is */
/*                   exactly zero; the UdU' (or LdL') factorization has */
/*                   been completed, but the block diagonal matrix d is */
/*                   exactly singular; */
/*             = M+1:  if the matrix R is numerically singular. */

/*     METHOD */
/*                            -     - */
/*     The matrices G, and/or A and Q are evaluated using the given or */
/*     computed symmetric factorization of R. */

/*     NUMERICAL ASPECTS */

/*     The routine should not be used when R is ill-conditioned. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Sep. 1997. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Algebraic Riccati equation, closed loop system, continuous-time */
/*     system, discrete-time system, optimal regulator, Schur form. */

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
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    q_dim1 = *ldq;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    r_dim1 = *ldr;
    r_offset = 1 + r_dim1;
    r__ -= r_offset;
    l_dim1 = *ldl;
    l_offset = 1 + l_dim1;
    l -= l_offset;
    --ipiv;
    g_dim1 = *ldg;
    g_offset = 1 + g_dim1;
    g -= g_offset;
    --iwork;
    --dwork;

    /* Function Body */
    *info = 0;
    ljobg = lsame(jobg, "G");
    ljobl = lsame(jobl, "N");
    lfactc = lsame(fact, "C");
    lfactu = lsame(fact, "U");
    luplou = lsame(uplo, "U");
    lfacta = lfactc || lfactu;

/*     Test the input scalar arguments. */

    if (! ljobg && ! lsame(jobg, "N")) {
	*info = -1;
    } else if (! ljobl && ! lsame(jobl, "Z")) {
	*info = -2;
    } else if (! lfacta && ! lsame(fact, "N")) {
	*info = -3;
    } else if (! luplou && ! lsame(uplo, "L")) {
	*info = -4;
    } else if (*n < 0) {
	*info = -5;
    } else if (*m < 0) {
	*info = -6;
    } else if (*lda < 1 || ljobl && *lda < *n) {
	*info = -8;
    } else if (*ldb < max(1,*n)) {
	*info = -10;
    } else if (*ldq < 1 || ljobl && *ldq < *n) {
	*info = -12;
    } else if (*ldr < max(1,*m)) {
	*info = -14;
    } else if (*ldl < 1 || ljobl && *ldl < *n) {
	*info = -16;
    } else if (*ldg < 1 || ljobg && *ldg < *n) {
	*info = -20;
    } else /* if(complicated condition) */ {
/* Computing MAX */
	i__1 = 1, i__2 = *n * *m;
/* Computing MAX */
	i__3 = 2, i__4 = *n * *m, i__3 = max(i__3,i__4), i__4 = *m * 3;
	if (lfactc && *ldwork < 1 || lfactu && *ldwork < max(i__1,i__2) || ! 
		lfacta && *ldwork < max(i__3,i__4)) {
	    *info = -23;
	}
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla("SB02MT", &i__1);
	return 0;
    }

    if (lfactc) {
	*oufact = 1;
    } else if (lfactu) {
	*oufact = 2;
    }

/*     Quick return if possible. */

    if (*n == 0 || *m == 0 || ! (ljobl || ljobg)) {
	dwork[1] = 1.;
	if (! lfacta) {
	    dwork[2] = 1.;
	}
	return 0;
    }

/*     (Note: Comments in the code beginning "Workspace:" describe the */
/*     minimal amount of workspace needed at that point in the code, */
/*     as well as the preferred amount for good performance. */
/*     NB refers to the optimal block size for the immediately */
/*     following subroutine, as returned by ILAENV.) */

    wrkopt = 1;

/*     Set relative machine precision. */

    eps = dlamch("Epsilon");

    if (! lfacta) {

/*        Compute the norm of the matrix R, which is not factored. */
/*        Then save the given triangle of R in the other strict triangle */
/*        and the diagonal in the workspace, and try Cholesky */
/*        factorization. */
/*        Workspace: need M. */

	rnorm = dlansy("1-norm", uplo, m, &r__[r_offset], ldr, &dwork[1]);
	i__1 = *ldr + 1;
	dcopy(m, &r__[r_offset], &i__1, &dwork[1], &c__1);
	if (luplou) {

	    i__1 = *m;
	    for (j = 2; j <= i__1; ++j) {
		i__2 = j - 1;
		dcopy(&i__2, &r__[j * r_dim1 + 1], &c__1, &r__[j + r_dim1], 
			ldr);
/* L20: */
	    }

	} else {

	    i__1 = *m;
	    for (j = 2; j <= i__1; ++j) {
		i__2 = j - 1;
		dcopy(&i__2, &r__[j + r_dim1], ldr, &r__[j * r_dim1 + 1], &
			c__1);
/* L40: */
	    }

	}
	dpotrf(uplo, m, &r__[r_offset], ldr, info);
	if (*info == 0) {

/*           Compute the reciprocal of the condition number of R. */
/*           Workspace: need 3*M. */

	    dpocon(uplo, m, &r__[r_offset], ldr, &rnorm, &rcond, &dwork[1], &
		    iwork[1], info);

/*           Return if the matrix is singular to working precision. */

	    *oufact = 1;
	    dwork[2] = rcond;
	    if (rcond < eps) {
		*info = *m + 1;
		return 0;
	    }
/* Computing MAX */
	    i__1 = wrkopt, i__2 = *m * 3;
	    wrkopt = max(i__1,i__2);
	} else {

/*           Use UdU' or LdL' factorization, first restoring the saved */
/*           triangle. */

	    i__1 = *ldr + 1;
	    dcopy(m, &dwork[1], &c__1, &r__[r_offset], &i__1);
	    if (luplou) {

		i__1 = *m;
		for (j = 2; j <= i__1; ++j) {
		    i__2 = j - 1;
		    dcopy(&i__2, &r__[j + r_dim1], ldr, &r__[j * r_dim1 + 1],
			     &c__1);
/* L60: */
		}

	    } else {

		i__1 = *m;
		for (j = 2; j <= i__1; ++j) {
		    i__2 = j - 1;
		    dcopy(&i__2, &r__[j * r_dim1 + 1], &c__1, &r__[j + 
			    r_dim1], ldr);
/* L80: */
		}

	    }

/*           Compute the UdU' or LdL' factorization. */
/*           Workspace: need   1, */
/*                      prefer M*NB. */

	    dsytrf(uplo, m, &r__[r_offset], ldr, &ipiv[1], &dwork[1], ldwork,
		     info);
	    *oufact = 2;
	    if (*info > 0) {
		dwork[2] = 1.;
		return 0;
	    }
/* Computing MAX */
	    i__1 = wrkopt, i__2 = (integer) dwork[1];
	    wrkopt = max(i__1,i__2);

/*           Compute the reciprocal of the condition number of R. */
/*           Workspace: need 2*M. */

	    dsycon(uplo, m, &r__[r_offset], ldr, &ipiv[1], &rnorm, &rcond, &
		    dwork[1], &iwork[1], info);

/*           Return if the matrix is singular to working precision. */

	    dwork[2] = rcond;
	    if (rcond < eps) {
		*info = *m + 1;
		return 0;
	    }
	}
    }

    if (*oufact == 1) {

/*        Solve positive definite linear system(s). */

	if (luplou) {
	    *(unsigned char *)trans = 'N';
	} else {
	    *(unsigned char *)trans = 'T';
	}

/*        Solve the system X*U = B, overwriting B with X. */

	dtrsm("Right", uplo, trans, "Non-unit", n, m, &c_b28, &r__[r_offset],
		 ldr, &b[b_offset], ldb);

	if (ljobg) {
/*                                      -1 */
/*           Compute the matrix  G = B*R  *B', multiplying X*X' in G. */

	    dsyrk(uplo, "No transpose", n, m, &c_b28, &b[b_offset], ldb, &
		    c_b31, &g[g_offset], ldg);
	}

	if (ljobl) {

/*           Update matrices A and Q. */

/*           Solve the system Y*U = L, overwriting L with Y. */

	    dtrsm("Right", uplo, trans, "Non-unit", n, m, &c_b28, &r__[
		    r_offset], ldr, &l[l_offset], ldl);

/*           Compute A <- A - X*Y'. */

	    dgemm("No transpose", "Transpose", n, n, m, &c_b37, &b[b_offset],
		     ldb, &l[l_offset], ldl, &c_b28, &a[a_offset], lda);

/*           Compute Q <- Q - Y*Y'. */

	    dsyrk(uplo, "No transpose", n, m, &c_b37, &l[l_offset], ldl, &
		    c_b28, &q[q_offset], ldq);
	}
    } else {

/*        Solve indefinite linear system(s). */

/*        Solve the system UdU'*X = B' (or LdL'*X = B'). */
/*        Workspace: need N*M. */

	i__1 = *m;
	for (j = 1; j <= i__1; ++j) {
	    dcopy(n, &b[j * b_dim1 + 1], &c__1, &dwork[j], m);
/* L100: */
	}

	dsytrs(uplo, m, n, &r__[r_offset], ldr, &ipiv[1], &dwork[1], m, info);

	if (ljobg) {
/*                                                    -1 */
/*           Compute a triangle of the matrix  G = B*R  *B' = B*X. */

	    if (luplou) {
		i__ = 1;

		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    dgemv("No transpose", &j, m, &c_b28, &b[b_offset], ldb, &
			    dwork[i__], &c__1, &c_b31, &g[j * g_dim1 + 1], &
			    c__1);
		    i__ += *m;
/* L120: */
		}

	    } else {

		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    dgemv("Transpose", m, &j, &c_b28, &dwork[1], m, &b[j + 
			    b_dim1], ldb, &c_b31, &g[j + g_dim1], ldg);
/* L140: */
		}

	    }
	}

	if (ljobl) {

/*           Update matrices A and Q. */

/*           Solve the system UdU'*Y = L' (or LdL'*Y = L'). */

	    i__1 = *m;
	    for (j = 1; j <= i__1; ++j) {
		dcopy(n, &l[j * l_dim1 + 1], &c__1, &dwork[j], m);
/* L160: */
	    }

	    dsytrs(uplo, m, n, &r__[r_offset], ldr, &ipiv[1], &dwork[1], m, 
		    info);

/*           A <- A - B*Y. */

	    dgemm("No transpose", "No transpose", n, n, m, &c_b37, &b[
		    b_offset], ldb, &dwork[1], m, &c_b28, &a[a_offset], lda);
/*                                            -          -1 */
/*           Compute a triangle of the matrix Q = Q - L*R  *L' = Q - L*Y. */

	    if (luplou) {
		i__ = 1;

		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    dgemv("No transpose", &j, m, &c_b37, &l[l_offset], ldl, &
			    dwork[i__], &c__1, &c_b28, &q[j * q_dim1 + 1], &
			    c__1);
		    i__ += *m;
/* L180: */
		}

	    } else {

		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    dgemv("Transpose", m, &j, &c_b37, &dwork[1], m, &l[j + 
			    l_dim1], ldl, &c_b28, &q[j + q_dim1], ldq);
/* L200: */
		}

	    }
	}
    }

    dwork[1] = (doublereal) wrkopt;
    if (! lfacta) {
	dwork[2] = rcond;
    }

/* *** Last line of SB02MT *** */
    return 0;
} /* sb02mt_ */

