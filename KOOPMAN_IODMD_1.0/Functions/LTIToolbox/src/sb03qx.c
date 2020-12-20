/* SB03QX.f -- translated by f2c (version 20041007).
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

static doublereal c_b24 = 0.;
static doublereal c_b25 = 1.;
static doublereal c_b26 = .5;

/* Subroutine */ int sb03qx(char *trana, char *uplo, char *lyapun, integer *
	n, doublereal *xanorm, doublereal *t, integer *ldt, doublereal *u, 
	integer *ldu, doublereal *r__, integer *ldr, doublereal *ferr, 
	integer *iwork, doublereal *dwork, integer *ldwork, integer *info)
{
    /* System generated locals */
    integer r_dim1, r_offset, t_dim1, t_offset, u_dim1, u_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j, ij, nn;
    static doublereal est;
    static integer kase;
    static doublereal temp;
    static integer itmp, info2;
    static doublereal scale;
    static logical lower;
    static char uplow[1];
    static logical update;
    static char tranat[1];
    static logical notrna;


/*     RELEASE 4.0, WGS COPYRIGHT 1999. */

/*     PURPOSE */

/*     To estimate a forward error bound for the solution X of a real */
/*     continuous-time Lyapunov matrix equation, */

/*            op(A)'*X + X*op(A) = C, */

/*     where op(A) = A or A' (A**T) and C is symmetric (C = C**T). The */
/*     matrix A, the right hand side C, and the solution X are N-by-N. */
/*     An absolute residual matrix, which takes into account the rounding */
/*     errors in forming it, is given in the array R. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     TRANA   CHARACTER*1 */
/*             Specifies the form of op(A) to be used, as follows: */
/*             = 'N':  op(A) = A    (No transpose); */
/*             = 'T':  op(A) = A**T (Transpose); */
/*             = 'C':  op(A) = A**T (Conjugate transpose = Transpose). */

/*     UPLO    CHARACTER*1 */
/*             Specifies which part of the symmetric matrix R is to be */
/*             used, as follows: */
/*             = 'U':  Upper triangular part; */
/*             = 'L':  Lower triangular part. */

/*     LYAPUN  CHARACTER*1 */
/*             Specifies whether or not the original Lyapunov equations */
/*             should be solved, as follows: */
/*             = 'O':  Solve the original Lyapunov equations, updating */
/*                     the right-hand sides and solutions with the */
/*                     matrix U, e.g., X <-- U'*X*U; */
/*             = 'R':  Solve reduced Lyapunov equations only, without */
/*                     updating the right-hand sides and solutions. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrices A and R.  N >= 0. */

/*     XANORM  (input) DOUBLE PRECISION */
/*             The absolute (maximal) norm of the symmetric solution */
/*             matrix X of the Lyapunov equation.  XANORM >= 0. */

/*     T       (input) DOUBLE PRECISION array, dimension (LDT,N) */
/*             The leading N-by-N upper Hessenberg part of this array */
/*             must contain the upper quasi-triangular matrix T in Schur */
/*             canonical form from a Schur factorization of A. */

/*     LDT     INTEGER */
/*             The leading dimension of array T.  LDT >= MAX(1,N). */

/*     U       (input) DOUBLE PRECISION array, dimension (LDU,N) */
/*             The leading N-by-N part of this array must contain the */
/*             orthogonal matrix U from a real Schur factorization of A. */
/*             If LYAPUN = 'R', the array U is not referenced. */

/*     LDU     INTEGER */
/*             The leading dimension of array U. */
/*             LDU >= 1,        if LYAPUN = 'R'; */
/*             LDU >= MAX(1,N), if LYAPUN = 'O'. */

/*     R       (input/output) DOUBLE PRECISION array, dimension (LDR,N) */
/*             On entry, if UPLO = 'U', the leading N-by-N upper */
/*             triangular part of this array must contain the upper */
/*             triangular part of the absolute residual matrix R, with */
/*             bounds on rounding errors added. */
/*             On entry, if UPLO = 'L', the leading N-by-N lower */
/*             triangular part of this array must contain the lower */
/*             triangular part of the absolute residual matrix R, with */
/*             bounds on rounding errors added. */
/*             On exit, the leading N-by-N part of this array contains */
/*             the symmetric absolute residual matrix R (with bounds on */
/*             rounding errors added), fully stored. */

/*     LDR     INTEGER */
/*             The leading dimension of array R.  LDR >= MAX(1,N). */

/*     FERR    (output) DOUBLE PRECISION */
/*             An estimated forward error bound for the solution X. */
/*             If XTRUE is the true solution, FERR bounds the magnitude */
/*             of the largest entry in (X - XTRUE) divided by the */
/*             magnitude of the largest entry in X. */
/*             If N = 0 or XANORM = 0, FERR is set to 0, without any */
/*             calculations. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (N*N) */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK.  LDWORK >= 2*N*N. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = N+1:  if the matrices T and -T' have common or very */
/*                   close eigenvalues; perturbed values were used to */
/*                   solve Lyapunov equations (but the matrix T is */
/*                   unchanged). */

/*     METHOD */

/*     The forward error bound is estimated using a practical error bound */
/*     similar to the one proposed in [1], based on the 1-norm estimator */
/*     in [2]. */

/*     REFERENCES */

/*     [1] Higham, N.J. */
/*         Perturbation theory and backward error for AX-XB=C. */
/*         BIT, vol. 33, pp. 124-136, 1993. */

/*     [2] Higham, N.J. */
/*         FORTRAN codes for estimating the one-norm of a real or */
/*         complex matrix, with applications to condition estimation. */
/*         ACM Trans. Math. Softw., 14, pp. 381-396, 1988. */

/*     NUMERICAL ASPECTS */
/*                               3 */
/*     The algorithm requires 0(N ) operations. */

/*     FURTHER COMMENTS */

/*     The option LYAPUN = 'R' may occasionally produce slightly worse */
/*     or better estimates, and it is much faster than the option 'O'. */
/*     The routine can be also used as a final step in estimating a */
/*     forward error bound for the solution of a continuous-time */
/*     algebraic matrix Riccati equation. */

/*     CONTRIBUTOR */

/*     V. Sima, Research Institute for Informatics, Bucharest, Romania, */
/*     Oct. 1998. Partly based on DGLSVX (and then SB03QD) by P. Petkov, */
/*     Tech. University of Sofia, March 1998 (and December 1998). */

/*     REVISIONS */

/*     February 6, 1999, V. Sima, Katholieke Univ. Leuven, Belgium. */
/*     V. Sima, Research Institute for Informatics, Bucharest, Oct. 2004. */

/*     KEYWORDS */

/*     Lyapunov equation, orthogonal transformation, real Schur form. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Decode and Test input parameters. */

    /* Parameter adjustments */
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;
    u_dim1 = *ldu;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    r_dim1 = *ldr;
    r_offset = 1 + r_dim1;
    r__ -= r_offset;
    --iwork;
    --dwork;

    /* Function Body */
    notrna = lsame(trana, "N");
    update = lsame(lyapun, "O");

    nn = *n * *n;
    *info = 0;
    if (! (notrna || lsame(trana, "T") || lsame(trana, "C"))) {
	*info = -1;
    } else if (! (lsame(uplo, "L") || lsame(uplo, "U"))) {
	*info = -2;
    } else if (! (update || lsame(lyapun, "R"))) {
	*info = -3;
    } else if (*n < 0) {
	*info = -4;
    } else if (*xanorm < 0.) {
	*info = -5;
    } else if (*ldt < max(1,*n)) {
	*info = -7;
    } else if (*ldu < 1 || update && *ldu < *n) {
	*info = -9;
    } else if (*ldr < max(1,*n)) {
	*info = -11;
    } else if (*ldwork < nn << 1) {
	*info = -15;
    }

    if (*info != 0) {
	i__1 = -(*info);
	xerbla("SB03QX", &i__1);
	return 0;
    }

/*     Quick return if possible. */

    *ferr = 0.;
    if (*n == 0 || *xanorm == 0.) {
	return 0;
    }

    itmp = nn + 1;

    if (notrna) {
	*(unsigned char *)tranat = 'T';
    } else {
	*(unsigned char *)tranat = 'N';
    }

/*     Fill in the remaining triangle of the symmetric residual matrix. */

    ma02ed(uplo, n, &r__[r_offset], ldr);

    kase = 0;

/*     REPEAT */
L10:
    dlacon(&nn, &dwork[itmp], &dwork[1], &iwork[1], &est, &kase);
    if (kase != 0) {

/*        Select the triangular part of symmetric matrix to be used. */

	if (dlansy("1-norm", "Upper", n, &dwork[1], n, &dwork[itmp]) >= dlansy("1-norm", "Lower", n, &dwork[1], n, &
		dwork[itmp])) {
	    *(unsigned char *)uplow = 'U';
	    lower = FALSE_;
	} else {
	    *(unsigned char *)uplow = 'L';
	    lower = TRUE_;
	}

	if (kase == 2) {
	    ij = 0;
	    if (lower) {

/*              Scale the lower triangular part of symmetric matrix */
/*              by the residual matrix. */

		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *n;
		    for (i__ = j; i__ <= i__2; ++i__) {
			++ij;
			dwork[ij] *= r__[i__ + j * r_dim1];
/* L20: */
		    }
		    ij += j;
/* L30: */
		}
	    } else {

/*              Scale the upper triangular part of symmetric matrix */
/*              by the residual matrix. */

		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = j;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			++ij;
			dwork[ij] *= r__[i__ + j * r_dim1];
/* L40: */
		    }
		    ij = ij + *n - j;
/* L50: */
		}
	    }
	}

	if (update) {

/*           Transform the right-hand side: RHS := U'*RHS*U. */

	    mb01ru(uplow, "Transpose", n, n, &c_b24, &c_b25, &dwork[1], n, &
		    u[u_offset], ldu, &dwork[1], n, &dwork[itmp], &nn, &info2);
	    i__1 = *n + 1;
	    dscal(n, &c_b26, &dwork[1], &i__1);
	}
	ma02ed(uplow, n, &dwork[1], n);

	if (kase == 2) {

/*           Solve op(T)'*Y + Y*op(T) = scale*RHS. */

	    sb03my(trana, n, &t[t_offset], ldt, &dwork[1], n, &scale, &info2);
	} else {

/*           Solve op(T)*W + W*op(T)' = scale*RHS. */

	    sb03my(tranat, n, &t[t_offset], ldt, &dwork[1], n, &scale, &
		    info2);
	}

	if (info2 > 0) {
	    *info = *n + 1;
	}

	if (update) {

/*           Transform back to obtain the solution: Z := U*Z*U', with */
/*           Z = Y or Z = W. */

	    mb01ru(uplow, "No transpose", n, n, &c_b24, &c_b25, &dwork[1], n,
		     &u[u_offset], ldu, &dwork[1], n, &dwork[itmp], &nn, &
		    info2);
	    i__1 = *n + 1;
	    dscal(n, &c_b26, &dwork[1], &i__1);
	}

	if (kase == 1) {
	    ij = 0;
	    if (lower) {

/*              Scale the lower triangular part of symmetric matrix */
/*              by the residual matrix. */

		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *n;
		    for (i__ = j; i__ <= i__2; ++i__) {
			++ij;
			dwork[ij] *= r__[i__ + j * r_dim1];
/* L60: */
		    }
		    ij += j;
/* L70: */
		}
	    } else {

/*              Scale the upper triangular part of symmetric matrix */
/*              by the residual matrix. */

		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = j;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			++ij;
			dwork[ij] *= r__[i__ + j * r_dim1];
/* L80: */
		    }
		    ij = ij + *n - j;
/* L90: */
		}
	    }
	}

/*        Fill in the remaining triangle of the symmetric matrix. */

	ma02ed(uplow, n, &dwork[1], n);
	goto L10;
    }

/*     UNTIL KASE = 0 */

/*     Compute the estimate of the relative error. */

    temp = *xanorm * scale;
    if (temp > est) {
	*ferr = est / temp;
    } else {
	*ferr = 1.;
    }

    return 0;

/* *** Last line of SB03QX *** */
} /* sb03qx_ */

