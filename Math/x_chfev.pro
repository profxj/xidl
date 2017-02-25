;+ 
; NAME: 
; x_chfev   
;    Version 1.1
;
; PURPOSE:
;!! CHFEV evaluates a cubic polynomial given in Hermite form.
;
; CALLING SEQUENCE:
;   
; INPUTS:
;
; RETURNS:
;  log f(N,X)
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   Oct-2013 Written by JXP
;------------------------------------------------------------------------------
;!*****************************************************************************80
;!
;!! CHFEV evaluates a cubic polynomial given in Hermite form.
;!
;!  Discussion:
;!
;!    This routine evaluates a cubic polynomial given in Hermite form at an
;!    array of points.  While designed for use by SPLINE_PCHIP_VAL, it may
;!    be useful directly as an evaluator for a piecewise cubic
;!    Hermite function in applications, such as graphing, where
;;!    the interval is known in advance.
;!
;!    The cubic polynomial is determined by function values
;!    F1, F2 and derivatives D1, D2 on the interval [X1,X2].
;!
;!  Licensing:
;!
;!    This code is distributed under the GNU LGPL license.
;!
;!  Modified:
;!
;!    25 June 2008
;!
;!  Author:
;!
;!    Original FORTRAN77 version by Fred Fritsch.
;!    FORTRAN90 version by John Burkardt.
;!
;!  Reference:
;;!
;!    Fred Fritsch, Ralph Carlson,
;!    Monotone Piecewise Cubic Interpolation,
;!    SIAM Journal on Numerical Analysis,
;!    Volume 17, Number 2, April 1980, pages 238-246.
;!
;!    David Kahaner, Cleve Moler, Steven Nash,
;!    Numerical Methods and Software,
;!    Prentice Hall, 1989,
;!    ISBN: 0-13-627258-4,
;!    LC: TA345.K34.
;!
;!  Parameters:
;!
;!    Input, real ( kind = 8 ) X1, X2, the endpoints of the interval of
;!    definition of the cubic.  X1 and X2 must be distinct.
;!
;!    Input, real ( kind = 8 ) F1, F2, the values of the function at X1 and
;!    X2, respectively.
;!
;!    Input, real ( kind = 8 ) D1, D2, the derivative values at X1 and
;!    X2, respectively.
;!
;!    Input, integer ( kind = 4 ) NE, the number of evaluation points.
;!
;!    Input, real ( kind = 8 ) XE(NE), the points at which the function is to
;!    be evaluated.  If any of the XE are outside the interval
;!    [X1,X2], a warning error is returned in NEXT.
;!
;!    Output, real ( kind = 8 ) FE(NE), the value of the cubic function
;!    at the points XE.
;!
;!    Output, integer ( kind = 4 ) NEXT(2), indicates the number of
;!    extrapolation points:
;!    NEXT(1) = number of evaluation points to the left of interval.
;!    NEXT(2) = number of evaluation points to the right of interval.
;!
;!    Output, integer ( kind = 4 ) IERR, error flag.
;!    0, no errors.
;!    -1, NE < 1.
;!    -2, X1 == X2.
;!
;-
pro x_chfev, x1, x2, f1, f2, d1, d2, nex, xe, fe, next, ierr 

  compile_opt strictarr

  if ( nex LT 1 ) then begin
     ierr = -1
    ;write ( *, '(a)' ) ' '
    ;write ( *, '(a)' ) 'CHFEV - Fatal error!'
    ;write ( *, '(a)' ) '  Number of evaluation points is less than 1.'
    ;write ( *, '(a,i8)' ) '  NE = ', nex
     stop
  endif

  h = x2 - x1

  if ( h EQ 0.0D+00 ) then begin
     ierr = -2
    ;write ( *, '(a)' ) ' '
    ;write ( *, '(a)' ) 'CHFEV - Fatal error!'
                                ;write ( *, '(a)' ) '  The interval [X1,X2] is of zero length.'
     stop
  endif
;!
;!  Initialize.
;!
  ierr = 0
  next = intarr(2)
  next[0] = 0
  next[1] = 0
  xmi = min([0.0D+00, h], max=xma)
  ;xma = max([0.0D+00, h])
;!
;!  Compute cubic coefficients expanded about X1.
;!
  delta = ( f2 - f1 ) / h
  del1 = ( d1 - delta ) / h
  del2 = ( d2 - delta ) / h
  c2 = -( del1 + del1 + del2 )
  c3 = ( del1 + del2 ) / h
;!
;!  Evaluation loop.
;!
  for i = 0, nex-1 do begin
     x = xe[i] - x1
     fe[i] = f1 + x * ( d1 + x * ( c2 + x * c3 ) )
;!
;!  Count the extrapolation points.
;!
     if ( x LT xmi ) then begin
        next[0] = next[0] + 1
     endif
  
     if ( xma LT x ) then begin
        next[1] = next[1] + 1
     endif

  endfor

  return
end
