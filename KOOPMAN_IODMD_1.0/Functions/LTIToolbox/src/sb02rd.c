/* SB02RD.f -- translated by f2c (version 20041007).
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

static integer c__0 = 0;
static integer c__1 = 1;
static doublereal c_b61 = 0.;
static doublereal c_b65 = 1.;
static doublereal c_b67 = .5;
static doublereal c_b85 = -1.;

/* Subroutine */ int sb02rd(char *job, char *dico, char *hinv, char *trana, 
	char *uplo, char *scal, char *sort, char *fact, char *lyapun, integer 
	*n, doublereal *a, integer *lda, doublereal *t, integer *ldt, 
	doublereal *v, integer *ldv, doublereal *g, integer *ldg, doublereal *
	q, integer *ldq, doublereal *x, integer *ldx, doublereal *sep, 
	doublereal *rcond, doublereal *ferr, doublereal *wr, doublereal *wi, 
	doublereal *s, integer *lds, integer *iwork, doublereal *dwork, 
	integer *ldwork, logical *bwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, g_dim1, g_offset, q_dim1, q_offset, s_dim1, 
	    s_offset, t_dim1, t_offset, v_dim1, v_offset, x_dim1, x_offset, 
	    i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, n2, nn, iu, iw, np1, iwb, iwc, iwf, iwi, ldw, lwe, 
	    lwn, iwr, lws;
    static logical joba, jobc, jobe, jbxa, lscl;
    static char jobs[1];
    static integer ierr;
    static logical jobx;
    static char loup[1];
    static integer nrot;
    static logical lscal;
    static char equed[1];
    static logical discr;
    static logical lhinv;
    static doublereal gnorm;
    static logical luplo;
    static doublereal qnorm;
    static logical lsort;
    static doublereal rconda;
    static char lofact[1];
    static logical nofact;
    static logical update, colequ;
    static char tranat[1];
    static doublereal rcondu;
    static logical notrna;
    static doublereal pivota;
    static logical rowequ;
    static doublereal pivotu, wrkopt;


/*     RELEASE 4.0, WGS COPYRIGHT 1999. */

/*     PURPOSE */

/*     To solve for X either the continuous-time algebraic Riccati */
/*     equation */
/*                                          -1 */
/*        Q + op(A)'*X + X*op(A) - X*op(B)*R  op(B)'*X = 0,           (1) */

/*     or the discrete-time algebraic Riccati equation */
/*                                                                -1 */
/*        X = op(A)'*X*op(A) - op(A)'*X*op(B)*(R + op(B)'*X*op(B))  * */
/*                             op(B)'*X*op(A) + Q,                    (2) */

/*     where op(M) = M or M' (M**T), A, op(B), Q, and R are N-by-N, */
/*     N-by-M, N-by-N, and M-by-M matrices respectively, with Q symmetric */
/*     and R symmetric nonsingular; X is an N-by-N symmetric matrix. */
/*                           -1 */
/*     The matrix G = op(B)*R  *op(B)' must be provided on input, instead */
/*     of B and R, that is, the continuous-time equation */

/*        Q + op(A)'*X + X*op(A) - X*G*X = 0,                         (3) */

/*     or the discrete-time equation */
/*                                -1 */
/*        Q + op(A)'*X*(I_n + G*X)  *op(A) - X = 0,                   (4) */

/*     are solved, where G is an N-by-N symmetric matrix. SLICOT Library */
/*     routine SB02MT should be used to compute G, given B and R. SB02MT */
/*     also enables to solve Riccati equations corresponding to optimal */
/*     problems with coupling terms. */

/*     The routine also returns the computed values of the closed-loop */
/*     spectrum of the optimal system, i.e., the stable eigenvalues */
/*     lambda(1),...,lambda(N) of the corresponding Hamiltonian or */
/*     symplectic matrix associated to the optimal problem. It is assumed */
/*     that the matrices A, G, and Q are such that the associated */
/*     Hamiltonian or symplectic matrix has N stable eigenvalues, i.e., */
/*     with negative real parts, in the continuous-time case, and with */
/*     moduli less than one, in the discrete-time case. */

/*     Optionally, estimates of the conditioning and error bound on the */
/*     solution of the Riccati equation (3) or (4) are returned. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     JOB     CHARACTER*1 */
/*             Specifies the computation to be performed, as follows: */
/*             = 'X':  Compute the solution only; */
/*             = 'C':  Compute the reciprocal condition number only; */
/*             = 'E':  Compute the error bound only; */
/*             = 'A':  Compute all: the solution, reciprocal condition */
/*                     number, and the error bound. */

/*     DICO    CHARACTER*1 */
/*             Specifies the type of Riccati equation to be solved or */
/*             analyzed, as follows: */
/*             = 'C':  Equation (3), continuous-time case; */
/*             = 'D':  Equation (4), discrete-time case. */

/*     HINV    CHARACTER*1 */
/*             If DICO = 'D' and JOB = 'X' or JOB = 'A', specifies which */
/*             symplectic matrix is to be constructed, as follows: */
/*             = 'D':  The matrix H in (6) (see METHOD) is constructed; */
/*             = 'I':  The inverse of the matrix H in (6) is constructed. */
/*             HINV is not used if DICO = 'C', or JOB = 'C' or 'E'. */

/*     TRANA   CHARACTER*1 */
/*             Specifies the form of op(A) to be used, as follows: */
/*             = 'N':  op(A) = A    (No transpose); */
/*             = 'T':  op(A) = A**T (Transpose); */
/*             = 'C':  op(A) = A**T (Conjugate transpose = Transpose). */

/*     UPLO    CHARACTER*1 */
/*             Specifies which triangle of the matrices G and Q is */
/*             stored, as follows: */
/*             = 'U':  Upper triangle is stored; */
/*             = 'L':  Lower triangle is stored. */

/*     SCAL    CHARACTER*1 */
/*             If JOB = 'X' or JOB = 'A', specifies whether or not a */
/*             scaling strategy should be used, as follows: */
/*             = 'G':  General scaling should be used; */
/*             = 'N':  No scaling should be used. */
/*             SCAL is not used if JOB = 'C' or 'E'. */

/*     SORT    CHARACTER*1 */
/*             If JOB = 'X' or JOB = 'A', specifies which eigenvalues */
/*             should be obtained in the top of the Schur form, as */
/*             follows: */
/*             = 'S':  Stable   eigenvalues come first; */
/*             = 'U':  Unstable eigenvalues come first. */
/*             SORT is not used if JOB = 'C' or 'E'. */

/*     FACT    CHARACTER*1 */
/*             If JOB <> 'X', specifies whether or not a real Schur */
/*             factorization of the closed-loop system matrix Ac is */
/*             supplied on entry, as follows: */
/*             = 'F':  On entry, T and V contain the factors from a real */
/*                     Schur factorization of the matrix Ac; */
/*             = 'N':  A Schur factorization of Ac will be computed */
/*                     and the factors will be stored in T and V. */
/*             For a continuous-time system, the matrix Ac is given by */
/*                Ac = A - G*X, if TRANA = 'N', or */
/*                Ac = A - X*G, if TRANA = 'T' or 'C', */
/*             and for a discrete-time system, the matrix Ac is given by */
/*                Ac = inv(I_n + G*X)*A, if TRANA = 'N', or */
/*                Ac = A*inv(I_n + X*G), if TRANA = 'T' or 'C'. */
/*             FACT is not used if JOB = 'X'. */

/*     LYAPUN  CHARACTER*1 */
/*             If JOB <> 'X', specifies whether or not the original or */
/*             "reduced" Lyapunov equations should be solved for */
/*             estimating reciprocal condition number and/or the error */
/*             bound, as follows: */
/*             = 'O':  Solve the original Lyapunov equations, updating */
/*                     the right-hand sides and solutions with the */
/*                     matrix V, e.g., X <-- V'*X*V; */
/*             = 'R':  Solve reduced Lyapunov equations only, without */
/*                     updating the right-hand sides and solutions. */
/*                     This means that a real Schur form T of Ac appears */
/*                     in the equations, instead of Ac. */
/*             LYAPUN is not used if JOB = 'X'. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrices A, Q, G, and X.  N >= 0. */

/*     A       (input) DOUBLE PRECISION array, dimension (LDA,N) */
/*             If JOB = 'X' or JOB = 'A' or FACT = 'N' or LYAPUN = 'O', */
/*             the leading N-by-N part of this array must contain the */
/*             coefficient matrix A of the equation. */
/*             If JOB = 'C' or 'E' and FACT = 'F' and LYAPUN = 'R', A is */
/*             not referenced. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A. */
/*             LDA >= MAX(1,N), if JOB  = 'X' or JOB = 'A' or */
/*                                 FACT = 'N' or LYAPUN = 'O'. */
/*             LDA >= 1,        otherwise. */

/*     T       (input or output) DOUBLE PRECISION array, dimension */
/*             (LDT,N) */
/*             If JOB <> 'X' and FACT = 'F', then T is an input argument */
/*             and on entry, the leading N-by-N upper Hessenberg part of */
/*             this array must contain the upper quasi-triangular matrix */
/*             T in Schur canonical form from a Schur factorization of Ac */
/*             (see argument FACT). */
/*             If JOB <> 'X' and FACT = 'N', then T is an output argument */
/*             and on exit, if INFO = 0 or INFO = 7, the leading N-by-N */
/*             upper Hessenberg part of this array contains the upper */
/*             quasi-triangular matrix T in Schur canonical form from a */
/*             Schur factorization of Ac (see argument FACT). */
/*             If JOB = 'X', the array T is not referenced. */

/*     LDT     INTEGER */
/*             The leading dimension of the array T. */
/*             LDT >= 1,        if JOB =  'X'; */
/*             LDT >= MAX(1,N), if JOB <> 'X'. */

/*     V       (input or output) DOUBLE PRECISION array, dimension */
/*             (LDV,N) */
/*             If JOB <> 'X' and FACT = 'F', then V is an input argument */
/*             and on entry, the leading N-by-N part of this array must */
/*             contain the orthogonal matrix V from a real Schur */
/*             factorization of Ac (see argument FACT). */
/*             If JOB <> 'X' and FACT = 'N', then V is an output argument */
/*             and on exit, if INFO = 0 or INFO = 7, the leading N-by-N */
/*             part of this array contains the orthogonal N-by-N matrix */
/*             from a real Schur factorization of Ac (see argument FACT). */
/*             If JOB = 'X', the array V is not referenced. */

/*     LDV     INTEGER */
/*             The leading dimension of the array V. */
/*             LDV >= 1,        if JOB =  'X'; */
/*             LDV >= MAX(1,N), if JOB <> 'X'. */

/*     G       (input/output) DOUBLE PRECISION array, dimension (LDG,N) */
/*             On entry, the leading N-by-N upper triangular part (if */
/*             UPLO = 'U') or lower triangular part (if UPLO = 'L') of */
/*             this array must contain the upper triangular part or lower */
/*             triangular part, respectively, of the symmetric matrix G. */
/*             On exit, if JOB = 'X' and DICO = 'D', or JOB <> 'X' and */
/*             LYAPUN = 'R', the leading N-by-N part of this array */
/*             contains the symmetric matrix G fully stored. */
/*             If JOB <> 'X' and LYAPUN = 'R', this array is modified */
/*             internally, but restored on exit. */

/*     LDG     INTEGER */
/*             The leading dimension of the array G.  LDG >= MAX(1,N). */

/*     Q       (input/output) DOUBLE PRECISION array, dimension (LDQ,N) */
/*             On entry, the leading N-by-N upper triangular part (if */
/*             UPLO = 'U') or lower triangular part (if UPLO = 'L') of */
/*             this array must contain the upper triangular part or lower */
/*             triangular part, respectively, of the symmetric matrix Q. */
/*             On exit, if JOB = 'X' and DICO = 'D', or JOB <> 'X' and */
/*             LYAPUN = 'R', the leading N-by-N part of this array */
/*             contains the symmetric matrix Q fully stored. */
/*             If JOB <> 'X' and LYAPUN = 'R', this array is modified */
/*             internally, but restored on exit. */

/*     LDQ     INTEGER */
/*             The leading dimension of the array Q.  LDQ >= MAX(1,N). */

/*     X       (input or output) DOUBLE PRECISION array, dimension */
/*             (LDX,N) */
/*             If JOB = 'C' or JOB = 'E', then X is an input argument */
/*             and on entry, the leading N-by-N part of this array must */
/*             contain the symmetric solution matrix of the algebraic */
/*             Riccati equation. If LYAPUN = 'R', this array is modified */
/*             internally, but restored on exit; however, it could differ */
/*             from the input matrix at the round-off error level. */
/*             If JOB = 'X' or JOB = 'A', then X is an output argument */
/*             and on exit, if INFO = 0 or INFO >= 6, the leading N-by-N */
/*             part of this array contains the symmetric solution matrix */
/*             X of the algebraic Riccati equation. */

/*     LDX     INTEGER */
/*             The leading dimension of the array X.  LDX >= MAX(1,N). */

/*     SEP     (output) DOUBLE PRECISION */
/*             If JOB = 'C' or JOB = 'A', and INFO = 0 or INFO = 7, the */
/*             estimated quantity */
/*                sep(op(Ac),-op(Ac)'), if DICO = 'C', or */
/*                sepd(op(Ac),op(Ac)'), if DICO = 'D'. (See METHOD.) */
/*             If JOB = 'C' or JOB = 'A' and X = 0, or JOB = 'E', SEP is */
/*             not referenced. */
/*             If JOB = 'X', and INFO = 0, INFO = 5 or INFO = 7, */
/*             SEP contains the scaling factor used, which should */
/*             multiply the (2,1) submatrix of U to recover X from the */
/*             first N columns of U (see METHOD). If SCAL = 'N', SEP is */
/*             set to 1. */

/*     RCOND   (output) DOUBLE PRECISION */
/*             If JOB = 'C' or JOB = 'A', and INFO = 0 or INFO = 7, an */
/*             estimate of the reciprocal condition number of the */
/*             algebraic Riccati equation. */
/*             If N = 0 or X = 0, RCOND is set to 1 or 0, respectively. */
/*             If JOB = 'X', or JOB = 'E', RCOND is not referenced. */

/*     FERR    (output) DOUBLE PRECISION */
/*             If JOB = 'E' or JOB = 'A', and INFO = 0 or INFO = 7, an */
/*             estimated forward error bound for the solution X. If XTRUE */
/*             is the true solution, FERR bounds the magnitude of the */
/*             largest entry in (X - XTRUE) divided by the magnitude of */
/*             the largest entry in X. */
/*             If N = 0 or X = 0, FERR is set to 0. */
/*             If JOB = 'X', or JOB = 'C', FERR is not referenced. */

/*     WR      (output) DOUBLE PRECISION array, dimension (2*N) */
/*     WI      (output) DOUBLE PRECISION array, dimension (2*N) */
/*             If JOB = 'X' or JOB = 'A', and INFO = 0 or INFO >= 5, */
/*             these arrays contain the real and imaginary parts, */
/*             respectively, of the eigenvalues of the 2N-by-2N matrix S, */
/*             ordered as specified by SORT (except for the case */
/*             HINV = 'D', when the order is opposite to that specified */
/*             by SORT). The leading N elements of these arrays contain */
/*             the closed-loop spectrum of the system matrix Ac (see */
/*             argument FACT). Specifically, */
/*                lambda(k) = WR(k) + j*WI(k), for k = 1,2,...,N. */
/*             If JOB = 'C' or JOB = 'E', these arrays are not */
/*             referenced. */

/*     S       (output) DOUBLE PRECISION array, dimension (LDS,2*N) */
/*             If JOB = 'X' or JOB = 'A', and INFO = 0 or INFO >= 5, the */
/*             leading 2N-by-2N part of this array contains the ordered */
/*             real Schur form S of the (scaled, if SCAL = 'G') */
/*             Hamiltonian or symplectic matrix H. That is, */

/*                    ( S    S   ) */
/*                    (  11   12 ) */
/*                S = (          ), */
/*                    ( 0    S   ) */
/*                    (       22 ) */

/*             where S  , S   and S   are N-by-N matrices. */
/*                    11   12      22 */
/*             If JOB = 'C' or JOB = 'E', this array is not referenced. */

/*     LDS     INTEGER */
/*             The leading dimension of the array S. */
/*             LDS >= MAX(1,2*N), if JOB = 'X' or JOB = 'A'; */
/*             LDS >= 1,          if JOB = 'C' or JOB = 'E'. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (LIWORK) */
/*             LIWORK >= 2*N,          if JOB = 'X'; */
/*             LIWORK >= N*N,          if JOB = 'C' or JOB = 'E'; */
/*             LIWORK >= MAX(2*N,N*N), if JOB = 'A'. */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, or INFO = 7, DWORK(1) returns the */
/*             optimal value of LDWORK. If INFO = 0, or INFO >= 5, and */
/*             JOB = 'X', or JOB = 'A', then DWORK(2) returns an estimate */
/*             RCONDU of the reciprocal of the condition number (in the */
/*             1-norm) of the N-th order system of algebraic equations */
/*             from which the solution matrix X is obtained, and DWORK(3) */
/*             returns the reciprocal pivot growth factor for the LU */
/*             factorization of the coefficient matrix of that system */
/*             (see SLICOT Library routine MB02PD); if DWORK(3) is much */
/*             less than 1, then the computed X and RCONDU could be */
/*             unreliable. */
/*             If DICO = 'D', and JOB = 'X', or JOB = 'A', then DWORK(4) */
/*             returns the reciprocal condition number RCONDA of the */
/*             given matrix A, and DWORK(5) returns the reciprocal pivot */
/*             growth factor for A or for its leading columns, if A is */
/*             singular (see SLICOT Library routine MB02PD); if DWORK(5) */
/*             is much less than 1, then the computed S and RCONDA could */
/*             be unreliable. */
/*             On exit, if INFO = 0, or INFO >= 4, and JOB = 'X', the */
/*             elements DWORK(6:5+4*N*N) contain the 2*N-by-2*N */
/*             transformation matrix  U  which reduced the Hamiltonian or */
/*             symplectic matrix  H  to the ordered real Schur form  S. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= 5+MAX(1,4*N*N+8*N), if JOB = 'X' or JOB = 'A'; */
/*             This may also be used for JOB = 'C' or JOB = 'E', but */
/*             exact bounds are as follows: */
/*             LDWORK >= 5 + MAX(1,LWS,LWE) + LWN, where */
/*             LWS = 0,       if FACT = 'F' or  LYAPUN = 'R'; */
/*                 = 5*N,     if FACT = 'N' and LYAPUN = 'O' and */
/*                                              DICO = 'C' and JOB = 'C'; */
/*                 = 5*N+N*N, if FACT = 'N' and LYAPUN = 'O' and */
/*                                              DICO = 'C' and JOB = 'E'; */
/*                 = 5*N+N*N, if FACT = 'N' and LYAPUN = 'O' and */
/*                                              DICO = 'D'; */
/*             LWE = 2*N*N,                if DICO = 'C' and JOB = 'C'; */
/*                 = 4*N*N,                if DICO = 'C' and JOB = 'E'; */
/*                 = MAX(3,2*N*N) + N*N,   if DICO = 'D' and JOB = 'C'; */
/*                 = MAX(3,2*N*N) + 2*N*N, if DICO = 'D' and JOB = 'E'; */
/*             LWN = 0,   if LYAPUN = 'O' or   JOB = 'C'; */
/*                 = 2*N, if LYAPUN = 'R' and DICO = 'C' and JOB = 'E'; */
/*                 = 3*N, if LYAPUN = 'R' and DICO = 'D' and JOB = 'E'. */
/*             For optimum performance LDWORK should sometimes be larger. */

/*     BWORK   LOGICAL array, dimension (LBWORK) */
/*             LBWORK >= 2*N,          if JOB = 'X' or JOB = 'A'; */
/*             LBWORK >= 1,            if JOB = 'C' or JOB = 'E', and */
/*                                     FACT = 'N' and LYAPUN = 'R'; */
/*             LBWORK >= 0,            otherwise. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  if matrix A is (numerically) singular in discrete- */
/*                   time case; */
/*             = 2:  if the Hamiltonian or symplectic matrix H cannot be */
/*                   reduced to real Schur form; */
/*             = 3:  if the real Schur form of the Hamiltonian or */
/*                   symplectic matrix H cannot be appropriately ordered; */
/*             = 4:  if the Hamiltonian or symplectic matrix H has less */
/*                   than N stable eigenvalues; */
/*             = 5:  if the N-th order system of linear algebraic */
/*                   equations, from which the solution matrix X would */
/*                   be obtained, is singular to working precision; */
/*             = 6:  if the QR algorithm failed to complete the reduction */
/*                   of the matrix Ac to Schur canonical form, T; */
/*             = 7:  if T and -T' have some almost equal eigenvalues, if */
/*                   DICO = 'C', or T has almost reciprocal eigenvalues, */
/*                   if DICO = 'D'; perturbed values were used to solve */
/*                   Lyapunov equations, but the matrix T, if given (for */
/*                   FACT = 'F'), is unchanged. (This is a warning */
/*                   indicator.) */

/*     METHOD */

/*     The method used is the Schur vector approach proposed by Laub [1], */
/*     but with an optional scaling, which enhances the numerical */
/*     stability [6]. It is assumed that [A,B] is a stabilizable pair */
/*     (where for (3) or (4), B is any matrix such that B*B' = G with */
/*     rank(B) = rank(G)), and [E,A] is a detectable pair, where E is any */
/*     matrix such that E*E' = Q with rank(E) = rank(Q). Under these */
/*     assumptions, any of the algebraic Riccati equations (1)-(4) is */
/*     known to have a unique non-negative definite solution. See [2]. */
/*     Now consider the 2N-by-2N Hamiltonian or symplectic matrix */

/*                 ( op(A)   -G    ) */
/*            H =  (               ),                                 (5) */
/*                 (  -Q   -op(A)' ), */

/*     for continuous-time equation, and */
/*                         -1              -1 */
/*                 (  op(A)           op(A)  *G       ) */
/*            H =  (        -1                   -1   ),              (6) */
/*                 ( Q*op(A)     op(A)' + Q*op(A)  *G ) */

/*     for discrete-time equation, respectively, where */
/*                       -1 */
/*            G = op(B)*R  *op(B)'. */
/*     The assumptions guarantee that H in (5) has no pure imaginary */
/*     eigenvalues, and H in (6) has no eigenvalues on the unit circle. */
/*     If Y is an N-by-N matrix then there exists an orthogonal matrix U */
/*     such that U'*Y*U is an upper quasi-triangular matrix. Moreover, U */
/*     can be chosen so that the 2-by-2 and 1-by-1 diagonal blocks */
/*     (corresponding to the complex conjugate eigenvalues and real */
/*     eigenvalues respectively) appear in any desired order. This is the */
/*     ordered real Schur form. Thus, we can find an orthogonal */
/*     similarity transformation U which puts (5) or (6) in ordered real */
/*     Schur form */

/*            U'*H*U = S = (S(1,1)  S(1,2)) */
/*                         (  0     S(2,2)) */

/*     where S(i,j) is an N-by-N matrix and the eigenvalues of S(1,1) */
/*     have negative real parts in case of (5), or moduli greater than */
/*     one in case of (6). If U is conformably partitioned into four */
/*     N-by-N blocks */

/*               U = (U(1,1)  U(1,2)) */
/*                   (U(2,1)  U(2,2)) */

/*     with respect to the assumptions we then have */
/*     (a) U(1,1) is invertible and X = U(2,1)*inv(U(1,1)) solves (1), */
/*         (2), (3), or (4) with X = X' and non-negative definite; */
/*     (b) the eigenvalues of S(1,1) (if DICO = 'C') or S(2,2) (if */
/*         DICO = 'D') are equal to the eigenvalues of optimal system */
/*         (the 'closed-loop' spectrum). */

/*     [A,B] is stabilizable if there exists a matrix F such that (A-BF) */
/*     is stable. [E,A] is detectable if [A',E'] is stabilizable. */

/*     The condition number of a Riccati equation is estimated as */

/*     cond = ( norm(Theta)*norm(A) + norm(inv(Omega))*norm(Q) + */
/*                 norm(Pi)*norm(G) ) / norm(X), */

/*     where Omega, Theta and Pi are linear operators defined by */

/*     Omega(W) = op(Ac)'*W + W*op(Ac), */
/*     Theta(W) = inv(Omega(op(W)'*X + X*op(W))), */
/*        Pi(W) = inv(Omega(X*W*X)), */

/*     in the continuous-time case, and */

/*     Omega(W) = op(Ac)'*W*op(Ac) - W, */
/*     Theta(W) = inv(Omega(op(W)'*X*op(Ac) + op(Ac)'X*op(W))), */
/*        Pi(W) = inv(Omega(op(Ac)'*X*W*X*op(Ac))), */

/*     in the discrete-time case, and Ac has been defined (see argument */
/*     FACT). Details are given in the comments of SLICOT Library */
/*     routines SB02QD and SB02SD. */

/*     The routine estimates the quantities */

/*     sep(op(Ac),-op(Ac)') = 1 / norm(inv(Omega)), */
/*     sepd(op(Ac),op(Ac)') = 1 / norm(inv(Omega)), */

/*     norm(Theta) and norm(Pi) using 1-norm condition estimator. */

/*     The forward error bound is estimated using a practical error bound */
/*     similar to the one proposed in [5]. */

/*     REFERENCES */

/*     [1] Laub, A.J. */
/*         A Schur Method for Solving Algebraic Riccati equations. */
/*         IEEE Trans. Auto. Contr., AC-24, pp. 913-921, 1979. */

/*     [2] Wonham, W.M. */
/*         On a matrix Riccati equation of stochastic control. */
/*         SIAM J. Contr., 6, pp. 681-697, 1968. */

/*     [3] Sima, V. */
/*         Algorithms for Linear-Quadratic Optimization. */
/*         Pure and Applied Mathematics: A Series of Monographs and */
/*         Textbooks, vol. 200, Marcel Dekker, Inc., New York, 1996. */

/*     [4] Ghavimi, A.R. and Laub, A.J. */
/*         Backward error, sensitivity, and refinement of computed */
/*         solutions of algebraic Riccati equations. */
/*         Numerical Linear Algebra with Applications, vol. 2, pp. 29-49, */
/*         1995. */

/*     [5] Higham, N.J. */
/*         Perturbation theory and backward error for AX-XB=C. */
/*         BIT, vol. 33, pp. 124-136, 1993. */

/*     [6] Petkov, P.Hr., Konstantinov, M.M., and Mehrmann, V. */
/*         DGRSVX and DMSRIC: Fortran 77 subroutines for solving */
/*         continuous-time matrix algebraic Riccati equations with */
/*         condition and accuracy estimates. */
/*         Preprint SFB393/98-16, Fak. f. Mathematik, Tech. Univ. */
/*         Chemnitz, May 1998. */

/*     NUMERICAL ASPECTS */
/*                               3 */
/*     The algorithm requires 0(N ) operations. The solution accuracy */
/*     can be controlled by the output parameter FERR. */

/*     FURTHER COMMENTS */

/*     To obtain a stabilizing solution of the algebraic Riccati */
/*     equation for DICO = 'D', set SORT = 'U', if HINV = 'D', or set */
/*     SORT = 'S', if HINV = 'I'. */

/*     The routine can also compute the anti-stabilizing solutions of */
/*     the algebraic Riccati equations, by specifying */
/*         SORT = 'U' if DICO = 'D' and HINV = 'I', or DICO = 'C', or */
/*         SORT = 'S' if DICO = 'D' and HINV = 'D'. */

/*     Usually, the combinations HINV = 'D' and SORT = 'U', or HINV = 'I' */
/*     and SORT = 'U', for stabilizing and anti-stabilizing solutions, */
/*     respectively, will be faster then the other combinations [3]. */

/*     The option LYAPUN = 'R' may produce slightly worse or better */
/*     estimates, and it is faster than the option 'O'. */

/*     This routine is a functionally extended and more accurate */
/*     version of the SLICOT Library routine SB02MD. Transposed problems */
/*     can be dealt with as well. Iterative refinement is used whenever */
/*     useful to solve linear algebraic systems. Condition numbers and */
/*     error bounds on the solutions are optionally provided. */

/*     CONTRIBUTOR */

/*     V. Sima, Research Institute for Informatics, Bucharest, Apr. 1999. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Oct. 2001, */
/*     Dec. 2002, Oct. 2004. */

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

/*     Decode the input parameters. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;
    v_dim1 = *ldv;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    g_dim1 = *ldg;
    g_offset = 1 + g_dim1;
    g -= g_offset;
    q_dim1 = *ldq;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    x_dim1 = *ldx;
    x_offset = 1 + x_dim1;
    x -= x_offset;
    --wr;
    --wi;
    s_dim1 = *lds;
    s_offset = 1 + s_dim1;
    s -= s_offset;
    --iwork;
    --dwork;
    --bwork;

    /* Function Body */
    n2 = *n + *n;
    nn = *n * *n;
    np1 = *n + 1;
    *info = 0;
    joba = lsame(job, "A");
    jobc = lsame(job, "C");
    jobe = lsame(job, "E");
    jobx = lsame(job, "X");
    nofact = lsame(fact, "N");
    notrna = lsame(trana, "N");
    discr = lsame(dico, "D");
    luplo = lsame(uplo, "U");
    lscal = lsame(scal, "G");
    lsort = lsame(sort, "S");
    update = lsame(lyapun, "O");
    jbxa = jobx || joba;
    lhinv = FALSE_;
    if (discr && jbxa) {
	lhinv = lsame(hinv, "D");
    }

/*     Test the input scalar arguments. */

    if (! (jbxa || jobc || jobe)) {
	*info = -1;
    } else if (! (discr || lsame(dico, "C"))) {
	*info = -2;
    } else if (discr && jbxa) {
	if (! (lhinv || lsame(hinv, "I"))) {
	    *info = -3;
	}
    }
    if (*info == 0) {
	if (! (notrna || lsame(trana, "T") || lsame(
		trana, "C"))) {
	    *info = -4;
	} else if (! (luplo || lsame(uplo, "L"))) {
	    *info = -5;
	} else if (jbxa) {
	    if (! (lscal || lsame(scal, "N"))) {
		*info = -6;
	    } else if (! (lsort || lsame(sort, "U"))) {
		*info = -7;
	    }
	}
	if (*info == 0 && ! jobx) {
	    if (! (nofact || lsame(fact, "F"))) {
		*info = -8;
	    } else if (! (update || lsame(lyapun, "R"))
		    ) {
		*info = -9;
	    }
	}
	if (*info == 0) {
	    if (*n < 0) {
		*info = -10;
	    } else if (*lda < 1 || (jbxa || nofact || update) && *lda < *n) {
		*info = -12;
	    } else if (*ldt < 1 || ! jobx && *ldt < *n) {
		*info = -14;
	    } else if (*ldv < 1 || ! jobx && *ldv < *n) {
		*info = -16;
	    } else if (*ldg < max(1,*n)) {
		*info = -18;
	    } else if (*ldq < max(1,*n)) {
		*info = -20;
	    } else if (*ldx < max(1,*n)) {
		*info = -22;
	    } else if (*lds < 1 || jbxa && *lds < n2) {
		*info = -29;
	    } else {
		if (jbxa) {
/* Computing MAX */
		    i__1 = 1, i__2 = (nn << 2) + (*n << 3);
		    if (*ldwork < max(i__1,i__2) + 5) {
			*info = -32;
		    }
		} else {
		    if (nofact && update) {
			if (! discr && jobc) {
			    lws = *n * 5;
			} else {
			    lws = *n * 5 + nn;
			}
		    } else {
			lws = 0;
		    }
		    if (discr) {
			if (jobc) {
/* Computing MAX */
			    i__1 = 3, i__2 = nn << 1;
			    lwe = max(i__1,i__2) + nn;
			} else {
/* Computing MAX */
			    i__1 = 3, i__2 = nn << 1;
			    lwe = max(i__1,i__2) + (nn << 1);
			}
		    } else {
			if (jobc) {
			    lwe = nn << 1;
			} else {
			    lwe = nn << 2;
			}
		    }
		    if (update || jobc) {
			lwn = 0;
		    } else {
			if (discr) {
			    lwn = *n * 3;
			} else {
			    lwn = *n << 1;
			}
		    }
/* Computing MAX */
		    i__1 = max(1,lws);
		    if (*ldwork < max(i__1,lwe) + 5 + lwn) {
			*info = -32;
		    }
		}
	    }
	}
    }

    if (*info != 0) {

/*        Error return. */

	i__1 = -(*info);
	xerbla("SB02RD", &i__1);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0) {
	if (jobx) {
	    *sep = 1.;
	}
	if (jobc || joba) {
	    *rcond = 1.;
	}
	if (jobe || joba) {
	    *ferr = 0.;
	}
	dwork[1] = 1.;
	dwork[2] = 1.;
	dwork[3] = 1.;
	if (discr) {
	    dwork[4] = 1.;
	    dwork[5] = 1.;
	}
	return 0;
    }

    if (jbxa) {

/*        Compute the solution matrix X. */

/*        Initialise the Hamiltonian or symplectic matrix associated with */
/*        the problem. */
/*        Workspace:  need   0    if DICO = 'C'; */
/*                           6*N, if DICO = 'D'. */

	sb02ru(dico, hinv, trana, uplo, n, &a[a_offset], lda, &g[g_offset], 
		ldg, &q[q_offset], ldq, &s[s_offset], lds, &iwork[1], &dwork[
		1], ldwork, &ierr);

	if (ierr != 0) {
	    *info = 1;
	    if (discr) {
		dwork[4] = dwork[1];
		dwork[5] = dwork[2];
	    }
	    return 0;
	}

	if (discr) {
	    wrkopt = (doublereal) (*n * 6);
	    rconda = dwork[1];
	    pivota = dwork[2];
	} else {
	    wrkopt = 0.;
	}

	if (lscal) {

/*           Scale the Hamiltonian or symplectic matrix S, using the */
/*           square roots of the norms of the matrices Q and G. */

	    qnorm = sqrt(dlansy("1-norm", uplo, n, &q[q_offset], ldq, &dwork[
		    1]));
	    gnorm = sqrt(dlansy("1-norm", uplo, n, &g[g_offset], ldg, &dwork[
		    1]));

	    lscl = qnorm > gnorm && gnorm > 0.;
	    if (lscl) {
		dlascl("G", &c__0, &c__0, &qnorm, &gnorm, n, n, &s[np1 + 
			s_dim1], lds, &ierr);
		dlascl("G", &c__0, &c__0, &gnorm, &qnorm, n, n, &s[np1 * 
			s_dim1 + 1], lds, &ierr);
	    }
	} else {
	    lscl = FALSE_;
	}

/*        Find the ordered Schur factorization of S,  S = U*H*U'. */
/*        Workspace:  need   5 + 4*N*N + 6*N; */
/*                    prefer larger. */

	iu = 6;
	iw = iu + (nn << 2);
	ldw = *ldwork - iw + 1;
	if (! discr) {
	    if (lsort) {
		dgees("Vectors", "Sorted", (L_fp)sb02mv, &n2, &s[s_offset], 
			lds, &nrot, &wr[1], &wi[1], &dwork[iu], &n2, &dwork[
			iw], &ldw, &bwork[1], &ierr);
	    } else {
		dgees("Vectors", "Sorted", (L_fp)sb02mr, &n2, &s[s_offset], 
			lds, &nrot, &wr[1], &wi[1], &dwork[iu], &n2, &dwork[
			iw], &ldw, &bwork[1], &ierr);
	    }
	} else {
	    if (lsort) {
		dgees("Vectors", "Sorted", (L_fp)sb02mw, &n2, &s[s_offset], 
			lds, &nrot, &wr[1], &wi[1], &dwork[iu], &n2, &dwork[
			iw], &ldw, &bwork[1], &ierr);
	    } else {
		dgees("Vectors", "Sorted", (L_fp)sb02ms, &n2, &s[s_offset], 
			lds, &nrot, &wr[1], &wi[1], &dwork[iu], &n2, &dwork[
			iw], &ldw, &bwork[1], &ierr);
	    }
	    if (lhinv) {
		dswap(n, &wr[1], &c__1, &wr[np1], &c__1);
		dswap(n, &wi[1], &c__1, &wi[np1], &c__1);
	    }
	}
	if (ierr > n2) {
	    *info = 3;
	} else if (ierr > 0) {
	    *info = 2;
	} else if (nrot != *n) {
	    *info = 4;
	}
	if (*info != 0) {
	    if (discr) {
		dwork[4] = rconda;
		dwork[5] = pivota;
	    }
	    return 0;
	}

/* Computing MAX */
	d__1 = wrkopt, d__2 = dwork[iw] + (doublereal) (iw - 1);
	wrkopt = max(d__1,d__2);

/*        Compute the solution of X*U(1,1) = U(2,1) using */
/*        LU factorization and iterative refinement. The (2,1) block of S */
/*        is used as a workspace for factoring U(1,1). */
/*        Workspace:  need   5 + 4*N*N + 8*N. */

/*        First transpose U(2,1) in-situ. */

	i__1 = *n - 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = *n - i__;
	    dswap(&i__2, &dwork[iu + *n + i__ * (n2 + 1) - 1], &n2, &dwork[
		    iu + *n + (i__ - 1) * (n2 + 1) + 1], &c__1);
/* L20: */
	}

	iwr = iw;
	iwc = iwr + *n;
	iwf = iwc + *n;
	iwb = iwf + *n;
	iw = iwb + *n;

	mb02pd("Equilibrate", "Transpose", n, n, &dwork[iu], &n2, &s[np1 + 
		s_dim1], lds, &iwork[1], equed, &dwork[iwr], &dwork[iwc], &
		dwork[iu + *n], &n2, &x[x_offset], ldx, &rcondu, &dwork[iwf], 
		&dwork[iwb], &iwork[np1], &dwork[iw], &ierr);
	if (jobx) {

/*           Restore U(2,1) back in-situ. */

	    i__1 = *n - 1;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		i__2 = *n - i__;
		dswap(&i__2, &dwork[iu + *n + i__ * (n2 + 1) - 1], &n2, &
			dwork[iu + *n + (i__ - 1) * (n2 + 1) + 1], &c__1);
/* L40: */
	    }

	    if (! lsame(equed, "N")) {

/*              Undo the equilibration of U(1,1) and U(2,1). */

		rowequ = lsame(equed, "R") || lsame(
			equed, "B");
		colequ = lsame(equed, "C") || lsame(
			equed, "B");

		if (rowequ) {

		    i__1 = *n;
		    for (i__ = 1; i__ <= i__1; ++i__) {
			dwork[iwr + i__ - 1] = 1. / dwork[iwr + i__ - 1];
/* L60: */
		    }

		    mb01sd("Row scaling", n, n, &dwork[iu], &n2, &dwork[iwr],
			     &dwork[iwc]);
		}

		if (colequ) {

		    i__1 = *n;
		    for (i__ = 1; i__ <= i__1; ++i__) {
			dwork[iwc + i__ - 1] = 1. / dwork[iwc + i__ - 1];
/* L80: */
		    }

		    mb01sd("Column scaling", n, n, &dwork[iu], &n2, &dwork[
			    iwr], &dwork[iwc]);
		    mb01sd("Column scaling", n, n, &dwork[iu + *n], &n2, &
			    dwork[iwr], &dwork[iwc]);
		}
	    }

/*           Set S(2,1) to zero. */

	    dlaset("Full", n, n, &c_b61, &c_b61, &s[np1 + s_dim1], lds);
	}

	pivotu = dwork[iw];

	if (ierr > 0) {

/*           Singular matrix. Set INFO and DWORK for error return. */

	    *info = 5;
	    goto L160;
	}

/*        Make sure the solution matrix X is symmetric. */

	i__1 = *n - 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = *n - i__;
	    daxpy(&i__2, &c_b65, &x[i__ + (i__ + 1) * x_dim1], ldx, &x[i__ + 
		    1 + i__ * x_dim1], &c__1);
	    i__2 = *n - i__;
	    dscal(&i__2, &c_b67, &x[i__ + 1 + i__ * x_dim1], &c__1);
	    i__2 = *n - i__;
	    dcopy(&i__2, &x[i__ + 1 + i__ * x_dim1], &c__1, &x[i__ + (i__ + 
		    1) * x_dim1], ldx);
/* L100: */
	}

	if (lscal) {

/*           Undo scaling for the solution matrix. */

	    if (lscl) {
		dlascl("G", &c__0, &c__0, &gnorm, &qnorm, n, n, &x[x_offset],
			 ldx, &ierr);
	    }
	}
    }

    if (! jobx) {
	if (! joba) {
	    wrkopt = 0.;
	}

/*        Estimate the conditioning and compute an error bound on the */
/*        solution of the algebraic Riccati equation. */

	iw = 6;
	*(unsigned char *)lofact = *(unsigned char *)fact;
	if (nofact && ! update) {

/*           Compute Ac and its Schur factorization. */

	    if (discr) {
		dlaset("Full", n, n, &c_b61, &c_b65, &dwork[iw], n);
		dsymm("Left", uplo, n, n, &c_b65, &g[g_offset], ldg, &x[
			x_offset], ldx, &c_b65, &dwork[iw], n);
		if (notrna) {

/*                 Compute Ac = inv(I_n + G*X)*A. */

		    dlacpy("Full", n, n, &a[a_offset], lda, &t[t_offset], 
			    ldt);
		    dgesv(n, n, &dwork[iw], n, &iwork[1], &t[t_offset], ldt, 
			    &ierr);
		} else {

/*                 Compute Ac = A*inv(I_n + X*G). */

		    ma02ad("Full", n, n, &a[a_offset], lda, &t[t_offset], 
			    ldt);
		    dgesv(n, n, &dwork[iw], n, &iwork[1], &t[t_offset], ldt, 
			    &ierr);
		    i__1 = *n;
		    for (i__ = 2; i__ <= i__1; ++i__) {
			i__2 = i__ - 1;
			dswap(&i__2, &t[i__ * t_dim1 + 1], &c__1, &t[i__ + 
				t_dim1], ldt);
/* L120: */
		    }
		}

	    } else {

		dlacpy("Full", n, n, &a[a_offset], lda, &t[t_offset], ldt);
		if (notrna) {

/*                 Compute Ac = A - G*X. */

		    dsymm("Left", uplo, n, n, &c_b85, &g[g_offset], ldg, &x[
			    x_offset], ldx, &c_b65, &t[t_offset], ldt);
		} else {

/*                 Compute Ac = A - X*G. */

		    dsymm("Right", uplo, n, n, &c_b85, &g[g_offset], ldg, &x[
			    x_offset], ldx, &c_b65, &t[t_offset], ldt);
		}
	    }

/*           Compute the Schur factorization of Ac, Ac = V*T*V'. */
/*           Workspace:  need   5 + 5*N. */
/*                       prefer larger. */

	    iwr = iw;
	    iwi = iwr + *n;
	    iw = iwi + *n;
	    ldw = *ldwork - iw + 1;

	    dgees("Vectors", "Not ordered", (L_fp)sb02ms, n, &t[t_offset], 
		    ldt, &nrot, &dwork[iwr], &dwork[iwi], &v[v_offset], ldv, &
		    dwork[iw], &ldw, &bwork[1], &ierr);

	    if (ierr != 0) {
		*info = 6;
		goto L160;
	    }

/* Computing MAX */
	    d__1 = wrkopt, d__2 = dwork[iw] + (doublereal) (iw - 1);
	    wrkopt = max(d__1,d__2);
	    *(unsigned char *)lofact = 'F';
	    iw = 6;
	}

	if (! update) {

/*           Update G, Q, and X using the orthogonal matrix V. */

	    *(unsigned char *)tranat = 'T';

/*           Save the diagonal elements of G and Q. */

	    i__1 = *ldg + 1;
	    dcopy(n, &g[g_offset], &i__1, &dwork[iw], &c__1);
	    i__1 = *ldq + 1;
	    dcopy(n, &q[q_offset], &i__1, &dwork[iw + *n], &c__1);
	    iw += n2;

	    if (joba) {
		dlacpy("Full", n, n, &x[x_offset], ldx, &s[np1 + s_dim1], 
			lds);
	    }
	    mb01ru(uplo, tranat, n, n, &c_b61, &c_b65, &x[x_offset], ldx, &v[
		    v_offset], ldv, &x[x_offset], ldx, &dwork[iw], &nn, &ierr);
	    i__1 = *ldx + 1;
	    dscal(n, &c_b67, &x[x_offset], &i__1);
	    ma02ed(uplo, n, &x[x_offset], ldx);
	    if (! discr) {
		ma02ed(uplo, n, &g[g_offset], ldg);
		ma02ed(uplo, n, &q[q_offset], ldq);
	    }
	    mb01ru(uplo, tranat, n, n, &c_b61, &c_b65, &g[g_offset], ldg, &v[
		    v_offset], ldv, &g[g_offset], ldg, &dwork[iw], &nn, &ierr);
	    i__1 = *ldg + 1;
	    dscal(n, &c_b67, &g[g_offset], &i__1);
	    mb01ru(uplo, tranat, n, n, &c_b61, &c_b65, &q[q_offset], ldq, &v[
		    v_offset], ldv, &q[q_offset], ldq, &dwork[iw], &nn, &ierr);
	    i__1 = *ldq + 1;
	    dscal(n, &c_b67, &q[q_offset], &i__1);
	}

/*        Estimate the conditioning and/or the error bound. */
/*        Workspace: 5 + MAX(1,LWS,LWE) + LWN, where */

/*           LWS = 0,       if FACT = 'F' or  LYAPUN = 'R'; */
/*               = 5*N,     if FACT = 'N' and LYAPUN = 'O' and DICO = 'C' */
/*                                                         and JOB = 'C'; */
/*               = 5*N+N*N, if FACT = 'N' and LYAPUN = 'O' and DICO = 'C' */
/*                                          and (JOB = 'E' or JOB = 'A'); */
/*               = 5*N+N*N, if FACT = 'N' and LYAPUN = 'O' and */
/*                                                         DICO = 'D'; */
/*           LWE = 2*N*N,                if DICO = 'C' and  JOB = 'C'; */
/*               = 4*N*N,                if DICO = 'C' and (JOB = 'E' or */
/*                                                          JOB = 'A'); */
/*               = MAX(3,2*N*N) + N*N,   if DICO = 'D' and  JOB = 'C'; */
/*               = MAX(3,2*N*N) + 2*N*N, if DICO = 'D' and (JOB = 'E' or */
/*                                                          JOB = 'A'); */
/*           LWN = 0,   if LYAPUN = 'O' or   JOB = 'C'; */
/*               = 2*N, if LYAPUN = 'R' and DICO = 'C' and (JOB = 'E' or */
/*                                                          JOB = 'A'); */
/*               = 3*N, if LYAPUN = 'R' and DICO = 'D' and (JOB = 'E' or */
/*                                                          JOB = 'A'). */

	ldw = *ldwork - iw + 1;
	if (joba) {
	    *(unsigned char *)jobs = 'B';
	} else {
	    *(unsigned char *)jobs = *(unsigned char *)job;
	}

	if (discr) {
	    sb02sd(jobs, lofact, trana, uplo, lyapun, n, &a[a_offset], lda, &
		    t[t_offset], ldt, &v[v_offset], ldv, &g[g_offset], ldg, &
		    q[q_offset], ldq, &x[x_offset], ldx, sep, rcond, ferr, &
		    iwork[1], &dwork[iw], &ldw, &ierr);
	} else {
	    sb02qd(jobs, lofact, trana, uplo, lyapun, n, &a[a_offset], lda, &
		    t[t_offset], ldt, &v[v_offset], ldv, &g[g_offset], ldg, &
		    q[q_offset], ldq, &x[x_offset], ldx, sep, rcond, ferr, &
		    iwork[1], &dwork[iw], &ldw, &ierr);
	}

/* Computing MAX */
	d__1 = wrkopt, d__2 = dwork[iw] + (doublereal) (iw - 1);
	wrkopt = max(d__1,d__2);
	if (ierr == np1) {
	    *info = 7;
	} else if (ierr > 0) {
	    *info = 6;
	    goto L160;
	}

	if (! update) {

/*           Restore X, G, and Q and set S(2,1) to zero, if needed. */

	    if (joba) {
		dlacpy("Full", n, n, &s[np1 + s_dim1], lds, &x[x_offset], ldx);
		dlaset("Full", n, n, &c_b61, &c_b61, &s[np1 + s_dim1], lds);
	    } else {
		mb01ru(uplo, trana, n, n, &c_b61, &c_b65, &x[x_offset], ldx, 
			&v[v_offset], ldv, &x[x_offset], ldx, &dwork[iw], &nn,
			 &ierr);
		i__1 = *ldx + 1;
		dscal(n, &c_b67, &x[x_offset], &i__1);
		ma02ed(uplo, n, &x[x_offset], ldx);
	    }
	    if (luplo) {
		*(unsigned char *)loup = 'L';
	    } else {
		*(unsigned char *)loup = 'U';
	    }

	    iw = 6;
	    i__1 = *ldg + 1;
	    dcopy(n, &dwork[iw], &c__1, &g[g_offset], &i__1);
	    ma02ed(loup, n, &g[g_offset], ldg);
	    i__1 = *ldq + 1;
	    dcopy(n, &dwork[iw + *n], &c__1, &q[q_offset], &i__1);
	    ma02ed(loup, n, &q[q_offset], ldq);
	}

    }

/*     Set the optimal workspace and other details. */

    dwork[1] = wrkopt;
L160:
    if (jbxa) {
	dwork[2] = rcondu;
	dwork[3] = pivotu;
	if (discr) {
	    dwork[4] = rconda;
	    dwork[5] = pivota;
	}
	if (jobx) {
	    if (lscl) {
		*sep = qnorm / gnorm;
	    } else {
		*sep = 1.;
	    }
	}
    }

    return 0;
/* *** Last line of SB02RD *** */
} /* sb02rd_ */

