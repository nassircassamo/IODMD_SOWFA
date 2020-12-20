/* MB03UD.f -- translated by f2c (version 20041007).
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
static integer c_n1 = -1;
static integer c__0 = 0;
static doublereal c_b32 = 0.;

/* Subroutine */ int mb03ud(char *jobq, char *jobp, integer *n, doublereal *
	a, integer *lda, doublereal *q, integer *ldq, doublereal *sv, 
	doublereal *dwork, integer *ldwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, q_dim1, q_offset, i__1, i__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, ie;
    static doublereal dum[1], eps;
    static integer iscl;
    static doublereal anrm;
    static integer ncolp, ncolq, itaup, itauq;
    static logical wantp, wantq;
    static integer jwork;
    static doublereal bignum;
    static integer minwrk, maxwrk;
    static doublereal smlnum;


/*     RELEASE 4.0, WGS COPYRIGHT 1999. */

/*     PURPOSE */

/*     To compute all, or part, of the singular value decomposition of a */
/*     real upper triangular matrix. */

/*     The N-by-N upper triangular matrix A is factored as  A = Q*S*P', */
/*     where Q and P are N-by-N orthogonal matrices and S is an */
/*     N-by-N diagonal matrix with non-negative diagonal elements, */
/*     SV(1), SV(2), ..., SV(N), ordered such that */

/*        SV(1) >= SV(2) >= ... >= SV(N) >= 0. */

/*     The columns of Q are the left singular vectors of A, the diagonal */
/*     elements of S are the singular values of A and the columns of P */
/*     are the right singular vectors of A. */

/*     Either or both of Q and P' may be requested. */
/*     When P' is computed, it is returned in A. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     JOBQ    CHARACTER*1 */
/*             Specifies whether the user wishes to compute the matrix Q */
/*             of left singular vectors as follows: */
/*             = 'V':  Left singular vectors are computed; */
/*             = 'N':  No left singular vectors are computed. */

/*     JOBP    CHARACTER*1 */
/*             Specifies whether the user wishes to compute the matrix P' */
/*             of right singular vectors as follows: */
/*             = 'V':  Right singular vectors are computed; */
/*             = 'N':  No right singular vectors are computed. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrix A.  N >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N upper triangular part of this */
/*             array must contain the upper triangular matrix A. */
/*             On exit, if JOBP = 'V', the leading N-by-N part of this */
/*             array contains the N-by-N orthogonal matrix  P'; otherwise */
/*             the N-by-N upper triangular part of A is used as internal */
/*             workspace. The strictly lower triangular part of A is set */
/*             internally to zero before the reduction to bidiagonal form */
/*             is performed. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,N). */

/*     Q       (output) DOUBLE PRECISION array, dimension (LDQ,N) */
/*             If JOBQ = 'V', the leading N-by-N part of this array */
/*             contains the orthogonal matrix Q. */
/*             If JOBQ = 'N', Q is not referenced. */

/*     LDQ     INTEGER */
/*             The leading dimension of array Q. */
/*             LDQ >= 1,  and when JOBQ = 'V',  LDQ >= MAX(1,N). */

/*     SV      (output) DOUBLE PRECISION array, dimension (N) */
/*             The N singular values of the matrix A, sorted in */
/*             descending order. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal LDWORK; */
/*             if INFO > 0, DWORK(2:N) contains the unconverged */
/*             superdiagonal elements of an upper bidiagonal matrix B */
/*             whose diagonal is in SV (not necessarily sorted). */
/*             B satisfies A = Q*B*P', so it has the same singular */
/*             values as A, and singular vectors related by Q and P'. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= MAX(1,5*N). */
/*             For optimum performance LDWORK should be larger. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             > 0:  the QR algorithm has failed to converge. In this */
/*                   case INFO specifies how many superdiagonals did not */
/*                   converge (see the description of DWORK). */
/*                   This failure is not likely to occur. */

/*     METHOD */

/*     The routine reduces A to bidiagonal form by means of elementary */
/*     reflectors and then uses the QR algorithm on the bidiagonal form. */

/*     CONTRIBUTOR */

/*     V. Sima, Research Institute of Informatics, Bucharest, and */
/*     A. Varga, German Aerospace Center, DLR Oberpfaffenhofen, */
/*     March 1998. Based on the RASP routine DTRSVD. */

/*     REVISIONS */

/*     V. Sima, Feb. 2000. */

/*     KEYWORDS */

/*     Bidiagonalization, orthogonal transformation, singular value */
/*     decomposition, singular values, triangular form. */

/*    ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. Local Arrays .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

/*     Check the input scalar arguments. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    q_dim1 = *ldq;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    --sv;
    --dwork;

    /* Function Body */
    *info = 0;
    wantq = lsame(jobq, "V");
    wantp = lsame(jobp, "V");
    minwrk = 1;
    if (! wantq && ! lsame(jobq, "N")) {
	*info = -1;
    } else if (! wantp && ! lsame(jobp, "N")) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*lda < max(1,*n)) {
	*info = -5;
    } else if (wantq && *ldq < max(1,*n) || ! wantq && *ldq < 1) {
	*info = -7;
    }

/*     Compute workspace */
/*     (Note: Comments in the code beginning "Workspace:" describe the */
/*     minimal amount of workspace needed at that point in the code, */
/*     as well as the preferred amount for good performance. */
/*     NB refers to the optimal block size for the immediately following */
/*     subroutine, as returned by ILAENV.) */

    if (*info == 0 && *ldwork >= 1 && *n > 0) {
	maxwrk = *n * 3 + (*n << 1) * ilaenv(&c__1, "DGEBRD", " ", n, n, &
		c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
	if (wantq) {
/* Computing MAX */
	    i__1 = maxwrk, i__2 = *n * 3 + *n * ilaenv(&c__1, "DORGBR", 
		    "Q", n, n, n, &c_n1, (ftnlen)6, (ftnlen)1);
	    maxwrk = max(i__1,i__2);
	}
	if (wantp) {
/* Computing MAX */
	    i__1 = maxwrk, i__2 = *n * 3 + *n * ilaenv(&c__1, "DORGBR", 
		    "P", n, n, n, &c_n1, (ftnlen)6, (ftnlen)1);
	    maxwrk = max(i__1,i__2);
	}
	minwrk = *n * 5;
	maxwrk = max(maxwrk,minwrk);
	dwork[1] = (doublereal) maxwrk;
    }

    if (*ldwork < minwrk) {
	*info = -10;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla("MB03UD", &i__1);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0) {
	dwork[1] = 1.;
	return 0;
    }

/*     Get machine constants. */

    eps = dlamch("P");
    smlnum = sqrt(dlamch("S")) / eps;
    bignum = 1. / smlnum;

/*     Scale A if max entry outside range [SMLNUM,BIGNUM]. */

    anrm = dlantr("Max", "Upper", "Non-unit", n, n, &a[a_offset], lda, dum);
    iscl = 0;
    if (anrm > 0. && anrm < smlnum) {
	iscl = 1;
	dlascl("Upper", &c__0, &c__0, &anrm, &smlnum, n, n, &a[a_offset], 
		lda, info);
    } else if (anrm > bignum) {
	iscl = 1;
	dlascl("Upper", &c__0, &c__0, &anrm, &bignum, n, n, &a[a_offset], 
		lda, info);
    }

/*     Zero out below. */

    if (*n > 1) {
	i__1 = *n - 1;
	i__2 = *n - 1;
	dlaset("Lower", &i__1, &i__2, &c_b32, &c_b32, &a[a_dim1 + 2], lda);
    }

/*     Find the singular values and optionally the singular vectors */
/*     of the upper triangular matrix A. */

    ie = 1;
    itauq = ie + *n;
    itaup = itauq + *n;
    jwork = itaup + *n;

/*     First reduce the matrix to bidiagonal form. The diagonal */
/*     elements will be in SV and the superdiagonals in DWORK(IE). */
/*     (Workspace: need 4*N, prefer 3*N+2*N*NB) */

    i__1 = *ldwork - jwork + 1;
    dgebrd(n, n, &a[a_offset], lda, &sv[1], &dwork[ie], &dwork[itauq], &
	    dwork[itaup], &dwork[jwork], &i__1, info);
    if (wantq) {

/*        Generate the transformation matrix Q corresponding to the */
/*        left singular vectors. */
/*        (Workspace: need 4*N, prefer 3*N+N*NB) */

	ncolq = *n;
	dlacpy("Lower", n, n, &a[a_offset], lda, &q[q_offset], ldq);
	i__1 = *ldwork - jwork + 1;
	dorgbr("Q", n, n, n, &q[q_offset], ldq, &dwork[itauq], &dwork[jwork],
		 &i__1, info);
    } else {
	ncolq = 0;
    }
    if (wantp) {

/*        Generate the transformation matrix P' corresponding to the */
/*        right singular vectors. */
/*        (Workspace: need 4*N, prefer 3*N+N*NB) */

	ncolp = *n;
	i__1 = *ldwork - jwork + 1;
	dorgbr("P", n, n, n, &a[a_offset], lda, &dwork[itaup], &dwork[jwork],
		 &i__1, info);
    } else {
	ncolp = 0;
    }
    jwork = ie + *n;

/*     Perform bidiagonal QR iteration, to obtain all or part of the */
/*     singular value decomposition of A. */
/*     (Workspace: need 5*N) */

    dbdsqr("U", n, &ncolp, &ncolq, &c__0, &sv[1], &dwork[ie], &a[a_offset], 
	    lda, &q[q_offset], ldq, dum, &c__1, &dwork[jwork], info);

/*     If DBDSQR failed to converge, copy unconverged superdiagonals */
/*     to DWORK(2:N). */

    if (*info != 0) {
	for (i__ = *n - 1; i__ >= 1; --i__) {
	    dwork[i__ + 1] = dwork[i__ + ie - 1];
/* L10: */
	}
    }

/*     Undo scaling if necessary. */

    if (iscl == 1) {
	if (anrm > bignum) {
	    dlascl("G", &c__0, &c__0, &bignum, &anrm, n, &c__1, &sv[1], n, info);
	}
	if (*info != 0 && anrm > bignum) {
	    i__1 = *n - 1;
	    dlascl("G", &c__0, &c__0, &bignum, &anrm, &i__1, &c__1, &dwork[2]
		    , n, info);
	}
	if (anrm < smlnum) {
	    dlascl("G", &c__0, &c__0, &smlnum, &anrm, n, &c__1, &sv[1], n, 
		    info);
	}
	if (*info != 0 && anrm < smlnum) {
	    i__1 = *n - 1;
	    dlascl("G", &c__0, &c__0, &smlnum, &anrm, &i__1, &c__1, &dwork[2]
		    , n, info);
	}
    }

/*     Return optimal workspace in DWORK(1). */

    dwork[1] = (doublereal) maxwrk;

    return 0;
/* *** Last line of MB03UD *** */
} /* mb03ud_ */

