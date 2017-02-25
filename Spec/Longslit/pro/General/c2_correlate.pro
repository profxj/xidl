; NAME:
;       C2_CORRELATE
;
; PURPOSE:
;       This function computes the cross correlation Pxy(L) or cross
;       covariance Rxy(L) of two sample populations X and Y as a function
;       of the lag (L).
;
; CATEGORY:
;       Statistics.
;
; CALLING SEQUENCE:
;       Result = c2_correlate(X, Y, Lag)
;
; INPUTS:
;       X:    An n-element vector of type integer, float or double.
;
;       Y:    An n-element vector of type integer, float or double.
;
;     LAG:    A scalar or n-element vector, in the interval [-(n-2), (n-2)],
;             of type integer that specifies the absolute distance(s) between
;             indexed elements of X.
;
; KEYWORD PARAMETERS:
;       COVARIANCE:    If set to a non-zero value, the sample cross
;                      covariance is computed.
;
;       DOUBLE:        If set to a non-zero value, computations are done in
;                      double precision arithmetic.
;
; EXAMPLE
;       Define two n-element sample populations.
;         x = [3.73, 3.67, 3.77, 3.83, 4.67, 5.87, 6.70, 6.97, 6.40, 5.57]
;         y = [2.31, 2.76, 3.02, 3.13, 3.72, 3.88, 3.97, 4.39, 4.34, 3.95]
;
;       Compute the cross correlation of X and Y for LAG = -5, 0, 1, 5, 6, 7
;         lag = [-5, 0, 1, 5, 6, 7]
;         result = bbc_correlate(x, y, lag)
;
;       The result should be:
;         [-0.428246, 0.914755, 0.674547, -0.405140, -0.403100, -0.339685]
;
; PROCEDURE:
;       See computational formula published in IDL manual.
;
; REFERENCE:
;       INTRODUCTION TO STATISTICAL TIME SERIES
;       Wayne A. Fuller
;       ISBN 0-471-28715-6
;
; MODIFICATION HISTORY:
;       Written by:  GGS, RSI, October 1994
;       Modified:    GGS, RSI, August 1995
;                    Corrected a condition which excluded the last term of the
;                    time-series.
;       Modified:    GGS, RSI, April 1996
;                    Simplified CROSS_COV function. Added DOUBLE keyword.
;                    Modified keyword checking and use of double precision.
;       Modified:    W. Biagiotti, July 1997
;                    Moved all constant calculations out of main loop for greatly
;                    reduced processing time.
;       Modified:    W. Biagiotti, Oct 1999
;                    Minor optimimization for yet more speed.


FUNCTION c2_correlate, X, Y, Lag, Covariance = Covariance, Double = Double


  ;Compute the sample cross correlation or cross covariance of
  ;(Xt, Xt+l) and (Yt, Yt+l) as a function of the lag (l).

  ON_ERROR, 2

  TypeX = SIZE(X)
  TypeY = SIZE(Y)
  nX = TypeX[TypeX[0]+2]
  nY = TypeY[TypeY[0]+2]

  if nX ne nY then $
    MESSAGE, "X and Y arrays must have the same number of elements."

  ;Check length.
  if nX lt 2 then $
    MESSAGE, "X and Y arrays must contain 2 or more elements."

  ;If the DOUBLE keyword is not set then the internal precision and
  ;result are identical to the type of input.
  if N_ELEMENTS(Double) eq 0 then $
    Double = (TypeX[TypeX[0]+1] eq 5 or TypeY[TypeY[0]+1] eq 5)

  nLag = N_ELEMENTS(Lag)

  IF (nLag EQ 1) THEN Lag = [Lag] ;Create a 1-element vector.

  IF (Double EQ 0) THEN Cross = FLTARR(nLag) $
  				   ELSE Cross = DBLARR(nLag)

  ; Calculate all constants out of the main processing loop
  Xmean = TOTAL(X, Double = Double) / nX
  Ymean = TOTAL(Y, Double = Double) / nX
  Xdata = X - Xmean
  Ydata = Y - Ymean
  nX1 = nX - 1L
  nY1 = nY - 1L

  FOR k = 0L, nLag-1 DO BEGIN
    IF (Lag[k] GE 0) $
      THEN Cross[k] = TOTAL((Xdata[0L:nX1 - Lag[k]]) * (Ydata[Lag[K]:nX1]), Double = Double) $
      ELSE Cross[k] = TOTAL((Ydata[0L:nY1 - ABS(Lag[k])]) * (Xdata[ABS(Lag[k]):nY1]), Double = Double)
  ENDFOR

  IF (KEYWORD_SET(Covariance) eq 0) THEN BEGIN ; Compute Cross Correlation.

     denom = SQRT((TOTAL((Xdata[0L:nX1])^2, Double = Double)) * $
     			  (TOTAL((Ydata[0L:nY1])^2, Double = Double)))

     Cross = Cross / denom
  ENDIF $
  ELSE Cross = Cross / nX  ; Covariance

  IF (Double EQ 0) THEN RETURN, FLOAT(Cross) $
  				   ELSE RETURN, Cross

END

