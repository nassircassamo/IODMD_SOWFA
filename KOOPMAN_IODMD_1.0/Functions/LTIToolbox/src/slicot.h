/* Starting from version 7.5, MATLAB BLAS is seperated */
#if MATLAB_VERSION >= 0x0705
#include <blas.h>
#endif
#include <lapack.h>

#if MATLAB_VERSION <= 0x0705
#define dcabs1 FORTRAN_WRAPPER(dcabs1)
extern doublereal dcabs1(
        doublereal *z
        );
		
#define ilaenv FORTRAN_WRAPPER(ilaenv)
extern integer ilaenv(
		integer *ispec,
		char *name,
		char *opts,
		integer *n1,
		integer *n2,
		integer *n3,
		integer *n4,
		integer name_len,
		integer opts_len
		);	
#endif

int ib01bd(char *meth, char *job, char *jobck, integer *
        nobr, integer *n, integer *m, integer *l, integer *nsmpl, doublereal *
        r__, integer *ldr, doublereal *a, integer *lda, doublereal *c__,
        integer *ldc, doublereal *b, integer *ldb, doublereal *d__, integer *
        ldd, doublereal *q, integer *ldq, doublereal *ry, integer *ldry,
        doublereal *s, integer *lds, doublereal *k, integer *ldk, doublereal *
        tol, integer *iwork, doublereal *dwork, integer *ldwork, logical *
        bwork, integer *iwarn, integer *info);

int ib01cd(char *jobx0, char *comuse, char *job, integer *n,
        integer *m, integer *l, integer *nsmp, doublereal *a, integer *lda,
        doublereal *b, integer *ldb, doublereal *c__, integer *ldc,
        doublereal *d__, integer *ldd, doublereal *u, integer *ldu,
        doublereal *y, integer *ldy, doublereal *x0, doublereal *v, integer *
        ldv, doublereal *tol, integer *iwork, doublereal *dwork, integer *
        ldwork, integer *iwarn, integer *info);

int ib01md(char *meth, char *alg, char *batch, char *conct,
        integer *nobr, integer *m, integer *l, integer *nsmp, doublereal *u,
        integer *ldu, doublereal *y, integer *ldy, doublereal *r__, integer *
        ldr, integer *iwork, doublereal *dwork, integer *ldwork, integer *
        iwarn, integer *info);

int ib01my(char *meth, char *batch, char *conct, integer *
        nobr, integer *m, integer *l, integer *nsmp, doublereal *u, integer *
        ldu, doublereal *y, integer *ldy, doublereal *r__, integer *ldr,
        integer *iwork, doublereal *dwork, integer *ldwork, integer *iwarn,
        integer *info);

int ib01nd(char *meth, char *jobd, integer *nobr, integer *
        m, integer *l, doublereal *r__, integer *ldr, doublereal *sv,
        doublereal *tol, integer *iwork, doublereal *dwork, integer *ldwork,
        integer *iwarn, integer *info);

int ib01pd(char *meth, char *job, char *jobcv, integer *
        nobr, integer *n, integer *m, integer *l, integer *nsmpl, doublereal *
        r__, integer *ldr, doublereal *a, integer *lda, doublereal *c__,
        integer *ldc, doublereal *b, integer *ldb, doublereal *d__, integer *
        ldd, doublereal *q, integer *ldq, doublereal *ry, integer *ldry,
        doublereal *s, integer *lds, doublereal *o, integer *ldo, doublereal *
        tol, integer *iwork, doublereal *dwork, integer *ldwork, integer *
        iwarn, integer *info);

int ib01px(char *job, integer *nobr, integer *n, integer *m,
        integer *l, doublereal *uf, integer *lduf, doublereal *un, integer *
        ldun, doublereal *ul, integer *ldul, doublereal *pgal, integer *
        ldpgal, doublereal *k, integer *ldk, doublereal *r__, integer *ldr,
        doublereal *x, doublereal *b, integer *ldb, doublereal *d__, integer *
        ldd, doublereal *tol, integer *iwork, doublereal *dwork, integer *
        ldwork, integer *iwarn, integer *info);

int ib01py(char *meth, char *job, integer *nobr, integer *n,
        integer *m, integer *l, integer *rankr1, doublereal *ul, integer *
        ldul, doublereal *r1, integer *ldr1, doublereal *tau1, doublereal *
        pgal, integer *ldpgal, doublereal *k, integer *ldk, doublereal *r__,
        integer *ldr, doublereal *h__, integer *ldh, doublereal *b, integer *
        ldb, doublereal *d__, integer *ldd, doublereal *tol, integer *iwork,
        doublereal *dwork, integer *ldwork, integer *iwarn, integer *info);

int ib01qd(char *jobx0, char *job, integer *n, integer *m,
        integer *l, integer *nsmp, doublereal *a, integer *lda, doublereal *
        c__, integer *ldc, doublereal *u, integer *ldu, doublereal *y,
        integer *ldy, doublereal *x0, doublereal *b, integer *ldb, doublereal
        *d__, integer *ldd, doublereal *tol, integer *iwork, doublereal *
        dwork, integer *ldwork, integer *iwarn, integer *info);

int ib01rd(char *job, integer *n, integer *m, integer *l,
        integer *nsmp, doublereal *a, integer *lda, doublereal *b, integer *
        ldb, doublereal *c__, integer *ldc, doublereal *d__, integer *ldd,
        doublereal *u, integer *ldu, doublereal *y, integer *ldy, doublereal *
        x0, doublereal *tol, integer *iwork, doublereal *dwork, integer *
        ldwork, integer *iwarn, integer *info);

int ma02ad(char *job, integer *m, integer *n, doublereal *a,
        integer *lda, doublereal *b, integer *ldb);

int ma02ed(char *uplo, integer *n, doublereal *a, integer *lda);

int ma02fd(doublereal *x1, doublereal *x2, doublereal *c__,
        doublereal *s, integer *info);

int mb01ru(char *uplo, char *trans, integer *m, integer *n,
        doublereal *alpha, doublereal *beta, doublereal *r__, integer *ldr,
        doublereal *a, integer *lda, doublereal *x, integer *ldx, doublereal *
        dwork, integer *ldwork, integer *info);

int mb01rx(char *side, char *uplo, char *trans, integer *m,
        integer *n, doublereal *alpha, doublereal *beta, doublereal *r__,
        integer *ldr, doublereal *a, integer *lda, doublereal *b, integer *
        ldb, integer *info);

int mb01ry(char *side, char *uplo, char *trans, integer *m,
        doublereal *alpha, doublereal *beta, doublereal *r__, integer *ldr,
        doublereal *h__, integer *ldh, doublereal *b, integer *ldb,
        doublereal *dwork, integer *info);

int mb01sd(char *jobs, integer *m, integer *n, doublereal *a,
        integer *lda, doublereal *r__, doublereal *c__);

int mb01td(integer *n, doublereal *a, integer *lda,
        doublereal *b, integer *ldb, doublereal *dwork, integer *info);

int mb01ud(char *side, char *trans, integer *m, integer *n,
        doublereal *alpha, doublereal *h__, integer *ldh, doublereal *a,
        integer *lda, doublereal *b, integer *ldb, integer *info);

int mb01vd(char *trana, char *tranb, integer *ma, integer *
        na, integer *mb, integer *nb, doublereal *alpha, doublereal *beta,
        doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *
        c__, integer *ldc, integer *mc, integer *nc, integer *info);

int mb02pd(char *fact, char *trans, integer *n, integer *
        nrhs, doublereal *a, integer *lda, doublereal *af, integer *ldaf,
        integer *ipiv, char *equed, doublereal *r__, doublereal *c__,
        doublereal *b, integer *ldb, doublereal *x, integer *ldx, doublereal *
        rcond, doublereal *ferr, doublereal *berr, integer *iwork, doublereal
        *dwork, integer *info);

int mb02qy(integer *m, integer *n, integer *nrhs, integer *
        rank, doublereal *a, integer *lda, integer *jpvt, doublereal *b,
        integer *ldb, doublereal *tau, doublereal *dwork, integer *ldwork,
        integer *info);

int mb02rz(char *trans, integer *n, integer *nrhs,
        doublecomplex *h__, integer *ldh, integer *ipiv, doublecomplex *b,
        integer *ldb, integer *info);

int mb02sz(integer *n, doublecomplex *h__, integer *ldh,
        integer *ipiv, integer *info);

int mb02tz(char *norm, integer *n, doublereal *hnorm,
        doublecomplex *h__, integer *ldh, integer *ipiv, doublereal *rcond,
        doublereal *dwork, doublecomplex *zwork, integer *info);

int mb02ud(char *fact, char *side, char *trans, char *jobp,
        integer *m, integer *n, doublereal *alpha, doublereal *rcond, integer
        *rank, doublereal *r__, integer *ldr, doublereal *q, integer *ldq,
        doublereal *sv, doublereal *b, integer *ldb, doublereal *rp, integer *
        ldrp, doublereal *dwork, integer *ldwork, integer *info);

int mb03od(char *jobqr, integer *m, integer *n, doublereal *a,
        integer *lda, integer *jpvt, doublereal *rcond, doublereal *svlmax,
        doublereal *tau, integer *rank, doublereal *sval, doublereal *dwork,
        integer *info);

int mb03ud(char *jobq, char *jobp, integer *n, doublereal *
        a, integer *lda, doublereal *q, integer *ldq, doublereal *sv,
        doublereal *dwork, integer *ldwork, integer *info);

int mb04id(integer *n, integer *m, integer *p, integer *l,
        doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *
        tau, doublereal *dwork, integer *ldwork, integer *info);

int mb04iy(char *side, char *trans, integer *n, integer *m,
        integer *k, integer *p, doublereal *a, integer *lda, doublereal *tau,
        doublereal *c__, integer *ldc, doublereal *dwork, integer *ldwork, integer *info);

int mb04kd(char *uplo, integer *n, integer *m, integer *p,
        doublereal *r__, integer *ldr, doublereal *a, integer *lda,
        doublereal *b, integer *ldb, doublereal *c__, integer *ldc,
        doublereal *tau, doublereal *dwork);

int mb04od(char *uplo, integer *n, integer *m, integer *p,
        doublereal *r__, integer *ldr, doublereal *a, integer *lda,
        doublereal *b, integer *ldb, doublereal *c__, integer *ldc,
        doublereal *tau, doublereal *dwork);

int mb04oy(integer *m, integer *n, doublereal *v,
        doublereal *tau, doublereal *a, integer *lda, doublereal *b, integer *
        ldb, doublereal *dwork);

logical sb02mr(doublereal *reig, doublereal *ieig);
logical sb02ms(doublereal *reig, doublereal *ieig);

int sb02mt(char *jobg, char *jobl, char *fact, char *uplo,
        integer *n, integer *m, doublereal *a, integer *lda, doublereal *b,
        integer *ldb, doublereal *q, integer *ldq, doublereal *r__, integer *
        ldr, doublereal *l, integer *ldl, integer *ipiv, integer *oufact,
        doublereal *g, integer *ldg, integer *iwork, doublereal *dwork,
        integer *ldwork, integer *info);

logical sb02mv(doublereal *reig, doublereal *ieig);
logical sb02mw(doublereal *reig, doublereal *ieig);

int sb02nd(char *dico, char *fact, char *uplo, char *jobl,
        integer *n, integer *m, integer *p, doublereal *a, integer *lda,
        doublereal *b, integer *ldb, doublereal *r__, integer *ldr, integer *
        ipiv, doublereal *l, integer *ldl, doublereal *x, integer *ldx,
        doublereal *rnorm, doublereal *f, integer *ldf, integer *oufact,
        integer *iwork, doublereal *dwork, integer *ldwork, integer *info);

int sb02qd(char *job, char *fact, char *trana, char *uplo,
        char *lyapun, integer *n, doublereal *a, integer *lda, doublereal *t,
        integer *ldt, doublereal *u, integer *ldu, doublereal *g, integer *
        ldg, doublereal *q, integer *ldq, doublereal *x, integer *ldx,
        doublereal *sep, doublereal *rcond, doublereal *ferr, integer *iwork,
        doublereal *dwork, integer *ldwork, integer *info);

int sb02rd(char *job, char *dico, char *hinv, char *trana,
        char *uplo, char *scal, char *sort, char *fact, char *lyapun, integer
        *n, doublereal *a, integer *lda, doublereal *t, integer *ldt,
        doublereal *v, integer *ldv, doublereal *g, integer *ldg, doublereal *
        q, integer *ldq, doublereal *x, integer *ldx, doublereal *sep,
        doublereal *rcond, doublereal *ferr, doublereal *wr, doublereal *wi,
        doublereal *s, integer *lds, integer *iwork, doublereal *dwork,
        integer *ldwork, logical *bwork, integer *info);

int sb02ru(char *dico, char *hinv, char *trana, char *uplo,
        integer *n, doublereal *a, integer *lda, doublereal *g, integer *ldg,
        doublereal *q, integer *ldq, doublereal *s, integer *lds, integer *
        iwork, doublereal *dwork, integer *ldwork, integer *info);

int sb02sd(char *job, char *fact, char *trana, char *uplo,
        char *lyapun, integer *n, doublereal *a, integer *lda, doublereal *t,
        integer *ldt, doublereal *u, integer *ldu, doublereal *g, integer *
        ldg, doublereal *q, integer *ldq, doublereal *x, integer *ldx,
        doublereal *sepd, doublereal *rcond, doublereal *ferr, integer *iwork,
        doublereal *dwork, integer *ldwork, integer *info);

int sb03mv(logical *ltran, logical *lupper, doublereal *t,
        integer *ldt, doublereal *b, integer *ldb, doublereal *scale,
        doublereal *x, integer *ldx, doublereal *xnorm, integer *info);

int sb03mw(logical *ltran, logical *lupper, doublereal *t,
        integer *ldt, doublereal *b, integer *ldb, doublereal *scale,
        doublereal *x, integer *ldx, doublereal *xnorm, integer *info);

int sb03mx(char *trana, integer *n, doublereal *a, integer *
        lda, doublereal *c__, integer *ldc, doublereal *scale, doublereal *
        dwork, integer *info);

int sb03my(char *trana, integer *n, doublereal *a, integer *
        lda, doublereal *c__, integer *ldc, doublereal *scale, integer *info);

int sb03ov(doublereal *a, doublereal *b, doublereal *c__, doublereal *s);

int sb03oy(logical *discr, logical *ltrans, integer *isgn,
        doublereal *s, integer *lds, doublereal *r__, integer *ldr,
        doublereal *a, integer *lda, doublereal *scale, integer *info);

int sb03qx(char *trana, char *uplo, char *lyapun, integer *
        n, doublereal *xanorm, doublereal *t, integer *ldt, doublereal *u,
        integer *ldu, doublereal *r__, integer *ldr, doublereal *ferr,
        integer *iwork, doublereal *dwork, integer *ldwork, integer *info);

int sb03qy(char *job, char *trana, char *lyapun, integer *n,
        doublereal *t, integer *ldt, doublereal *u, integer *ldu, doublereal
        *x, integer *ldx, doublereal *sep, doublereal *thnorm, integer *iwork,
        doublereal *dwork, integer *ldwork, integer *info);

int sb03sx(char *trana, char *uplo, char *lyapun, integer *
        n, doublereal *xanorm, doublereal *t, integer *ldt, doublereal *u,
        integer *ldu, doublereal *r__, integer *ldr, doublereal *ferr,
        integer *iwork, doublereal *dwork, integer *ldwork, integer *info);

int sb03sy(char *job, char *trana, char *lyapun, integer *n,
        doublereal *t, integer *ldt, doublereal *u, integer *ldu, doublereal
        *xa, integer *ldxa, doublereal *sepd, doublereal *thnorm, integer *
        iwork, doublereal *dwork, integer *ldwork, integer *info);

int sb04px(logical *ltranl, logical *ltranr, integer *isgn,
        integer *n1, integer *n2, doublereal *tl, integer *ldtl, doublereal *
        tr, integer *ldtr, doublereal *b, integer *ldb, doublereal *scale,
        doublereal *x, integer *ldx, doublereal *xnorm, integer *info);

int tb01wd(integer *n, integer *m, integer *p, doublereal *
        a, integer *lda, doublereal *b, integer *ldb, doublereal *c__,
        integer *ldc, doublereal *u, integer *ldu, doublereal *wr, doublereal
        *wi, doublereal *dwork, integer *ldwork, integer *info);

int tb05ad(char *baleig, char *inita, integer *n, integer *
        m, integer *p, doublecomplex *freq, doublereal *a, integer *lda,
        doublereal *b, integer *ldb, doublereal *c__, integer *ldc,
        doublereal *rcond, doublecomplex *g, integer *ldg, doublereal *evre,
        doublereal *evim, doublecomplex *hinvb, integer *ldhinv, integer *
        iwork, doublereal *dwork, integer *ldwork, doublecomplex *zwork,
        integer *lzwork, integer *info);

logical select1(doublereal *par1, doublereal *par2);

