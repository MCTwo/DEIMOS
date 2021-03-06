*DECK DERF
      DOUBLE PRECISION FUNCTION DERF (X)
C***BEGIN PROLOGUE  DERF
C***PURPOSE  Compute the error function.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C8A, L5A1E
C***TYPE      DOUBLE PRECISION (ERF-S, DERF-D)
C***KEYWORDS  ERF, ERROR FUNCTION, FNLIB, SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C DERF(X) calculates the double precision error function for double
C precision argument X.
C
C Series for ERF        on the interval  0.          to  1.00000E+00
C                                        with weighted error   1.28E-32
C                                         log weighted error  31.89
C                               significant figures required  31.05
C                                    decimal places required  32.55
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH, DCSEVL, DERFC, INITDS
C***REVISION HISTORY  (YYMMDD)
C   770701  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900727  Added EXTERNAL statement.  (WRB)
C   920618  Removed space from variable name.  (RWC, WRB)
C***END PROLOGUE  DERF
      DOUBLE PRECISION X, ERFCS(21), SQEPS, SQRTPI, XBIG, Y, D1MACH,
     1  DCSEVL, DERFC
      LOGICAL FIRST
      EXTERNAL DERFC
      SAVE ERFCS, SQRTPI, NTERF, XBIG, SQEPS, FIRST
      DATA ERFCS(  1) / -.4904612123 4691808039 9845440333 76 D-1     /
      DATA ERFCS(  2) / -.1422612051 0371364237 8247418996 31 D+0     /
      DATA ERFCS(  3) / +.1003558218 7599795575 7546767129 33 D-1     /
      DATA ERFCS(  4) / -.5768764699 7674847650 8270255091 67 D-3     /
      DATA ERFCS(  5) / +.2741993125 2196061034 4221607914 71 D-4     /
      DATA ERFCS(  6) / -.1104317550 7344507604 1353812959 05 D-5     /
      DATA ERFCS(  7) / +.3848875542 0345036949 9613114981 74 D-7     /
      DATA ERFCS(  8) / -.1180858253 3875466969 6317518015 81 D-8     /
      DATA ERFCS(  9) / +.3233421582 6050909646 4029309533 54 D-10    /
      DATA ERFCS( 10) / -.7991015947 0045487581 6073747085 95 D-12    /
      DATA ERFCS( 11) / +.1799072511 3961455611 9672454866 34 D-13    /
      DATA ERFCS( 12) / -.3718635487 8186926382 3168282094 93 D-15    /
      DATA ERFCS( 13) / +.7103599003 7142529711 6899083946 66 D-17    /
      DATA ERFCS( 14) / -.1261245511 9155225832 4954248533 33 D-18    /
      DATA ERFCS( 15) / +.2091640694 1769294369 1705002666 66 D-20    /
      DATA ERFCS( 16) / -.3253973102 9314072982 3641600000 00 D-22    /
      DATA ERFCS( 17) / +.4766867209 7976748332 3733333333 33 D-24    /
      DATA ERFCS( 18) / -.6598012078 2851343155 1999999999 99 D-26    /
      DATA ERFCS( 19) / +.8655011469 9637626197 3333333333 33 D-28    /
      DATA ERFCS( 20) / -.1078892517 7498064213 3333333333 33 D-29    /
      DATA ERFCS( 21) / +.1281188399 3017002666 6666666666 66 D-31    /
      DATA SQRTPI / 1.772453850 9055160272 9816748334 115D0 /
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  DERF
      IF (FIRST) THEN
         NTERF = INITDS (ERFCS, 21, 0.1*REAL(D1MACH(3)))
         XBIG = SQRT(-LOG(SQRTPI*D1MACH(3)))
         SQEPS = SQRT(2.0D0*D1MACH(3))
      ENDIF
      FIRST = .FALSE.
C
      Y = ABS(X)
      IF (Y.GT.1.D0) GO TO 20
C
C ERF(X) = 1.0 - ERFC(X)  FOR  -1.0 .LE. X .LE. 1.0
C
      IF (Y.LE.SQEPS) DERF = 2.0D0*X*X/SQRTPI
      IF (Y.GT.SQEPS) DERF = X*(1.0D0 + DCSEVL (2.D0*X*X-1.D0,
     1  ERFCS, NTERF))
      RETURN
C
C ERF(X) = 1.0 - ERFC(X) FOR ABS(X) .GT. 1.0
C
 20   IF (Y.LE.XBIG) DERF = SIGN (1.0D0-DERFC(Y), X)
      IF (Y.GT.XBIG) DERF = SIGN (1.0D0, X)
C
      RETURN
      END
