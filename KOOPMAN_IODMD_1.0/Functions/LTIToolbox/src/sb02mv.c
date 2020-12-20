/* SB02MV.f -- translated by f2c (version 20041007).
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

logical sb02mv(doublereal *reig, doublereal *ieig)
{
    /* System generated locals */
    logical ret_val;


/*     RELEASE 4.0, WGS COPYRIGHT 1999. */

/*     PURPOSE */

/*     To select the stable eigenvalues for solving the continuous-time */
/*     algebraic Riccati equation. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     REIG    (input) DOUBLE PRECISION */
/*             The real part of the current eigenvalue considered. */

/*     IEIG    (input) DOUBLE PRECISION */
/*             The imaginary part of the current eigenvalue considered. */

/*     METHOD */

/*     The function value SB02MV is set to .TRUE. for a stable eigenvalue */
/*     and to .FALSE., otherwise. */

/*     REFERENCES */

/*     None. */

/*     NUMERICAL ASPECTS */

/*     None. */

/*     CONTRIBUTOR */

/*     V. Sima, Katholieke Univ. Leuven, Belgium, Aug. 1997. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Algebraic Riccati equation, closed loop system, continuous-time */
/*     system, optimal regulator, Schur form. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Executable Statements .. */

    ret_val = *reig < 0.;

    return ret_val;
/* *** Last line of SB02MV *** */
} /* sb02mv_ */

