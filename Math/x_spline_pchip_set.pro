;+ 
; NAME: 
; x_spline_pchip_set   
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
;!! SPLINE_PCHIP_SET sets derivatives for a piecewise cubic Hermite interpolant.
;!
;!  Discussion:
;!
;!    This routine computes what would normally be called a Hermite
;!    interpolant.  However, the user is only required to supply function
;!    values, not derivative values as well.  This routine computes
;!    "suitable" derivative values, so that the resulting Hermite interpolant
;!    has desirable shape and monotonicity properties.
;!
;!    The interpolant will have an extremum at each point where
;!    monotonicity switches direction.
;!
;!    The resulting piecewise cubic Hermite function may be evaluated
;!    by SPLINE_PCHIP_VAL.
;!
;!    This routine was originally named "PCHIM".
;!
;!  Licensing:
;!
;!    This code is distributed under the GNU LGPL license.
;!
;!  Modified:
;!
;!    15 December 2008
;!
;!  Author:
;!
;!    Original FORTRAN77 version by Fred Fritsch.
;!    FORTRAN90 version by John Burkardt.
;!
;!  Reference:
;!
;!    Fred Fritsch, Ralph Carlson,
;!    Monotone Piecewise Cubic Interpolation,
;!    SIAM Journal on Numerical Analysis,
;!    Volume 17, Number 2, April 1980, pages 238-246.
;!
;!    Fred Fritsch, Judy Butland,
;!    A Method for Constructing Local Monotone Piecewise Cubic Interpolants,
;!    SIAM Journal on Scientific and Statistical Computing,
;!    Volume 5, Number 2, 1984, pages 300-304.
;!
;!  Parameters:
;!
;!    Input, integer ( kind = 4 ) N, the number of data points.  N must be
;!    at least 2.
;!
;!    Input, real ( kind = 8 ) X(N), the strictly increasing independent
;!    variable values.
;!
;!    Input, real ( kind = 8 ) F(N), dependent variable values to be
;!    interpolated.  F(I) is the value corresponding to X(I).
;!    This routine is designed for monotonic data, but it will work for any
;!    F array.  It will force extrema at points where monotonicity switches
;!    direction.
;!
;!    Output, real ( kind = 8 ) D(N), the derivative values at the
;!    data points.  If the data are monotonic, these values will determine
;!    a monotone cubic Hermite function.
;!
;-
function  x_spline_pchip_set, n, x, f

  compile_opt strictarr
  d = dblarr(n)

;!
;!  Check the arguments.
;!
  if ( n LT 2 ) then begin
     ierr = -1
     ;write ( *, '(a)' ) ' '
     ;write ( *, '(a)' ) 'SPLINE_PCHIP_SET - Fatal error!'
     ;write ( *, '(a)' ) '  Number of data points less than 2.'
     stop
  endif

  for i = 1, n-1 do begin
     if ( x[i] LE x[i-1] ) then begin
        ierr = -3
                                ;write ( *, '(a)' ) ' '
                                ;write ( *, '(a)' ) 'SPLINE_PCHIP_SET - Fatal error!'
                                ;write ( *, '(a)' ) '  X array not strictly increasing.'
        stop
     endif
  endfor

  ierr = 0
  nless1 = n - 1
  h1 = x[1] - x[0]
  del1 = ( f[1] - f[0] ) / h1
  dsave = del1
;!
;!  Special case N=2, use linear interpolation.
;!
  if ( n EQ 2 ) then begin
     d[0] = del1
     d[n-1] = del1
     return, d
  endif
;!
;!  Normal case, 3 <= N.
;!
  h2 = x[2] - x[1]
  del2 = ( f[2] - f[1] ) / h2
;!
;!  Set D[0] via non-centered three point formula, adjusted to be
;!  shape preserving.
;!
  hsum = h1 + h2
  w1 = ( h1 + hsum ) / hsum
  w2 = -h1 / hsum
  d[0] = w1 * del1 + w2 * del2

  ;if ( pchst ( d(1), del1 ) <= 0.0D+00 ) then
  if ( (d[0]*del1) LE  0.0D+00 ) then begin
     d[0] = 0.0D+00
  endif else begin
;!
;!  Need do this check only if monotonicity switches.
;!
     ;if ( pchst ( del1, del2 ) < 0.0D+00 ) then
     if ( (del1*del2)  LT 0.0D+00 ) then begin
        dmax = 3.0D+00 * del1

        if ( abs ( dmax ) LT abs ( d[0] ) ) then d[0] = dmax
     endif

  endelse 
;!
;!  Loop through interior points.
;!
  for i = 1, nless1-1 do begin

     if ( i GT 1 ) then begin
        h1 = h2
        h2 = x[i+1] - x[i]
        hsum = h1 + h2
        del1 = del2
        del2 = ( f[i+1] - f[i] ) / h2
     endif
;!
;!  Set D(I)=0 unless data are strictly monotonic.
;!
    d[i] = 0.0D+00

    temp = del1*del2 ; pchst ( del1, del2 )

    if ( temp LT 0.0D+00 ) then begin
       ierr = ierr + 1
       dsave = del2
;!
;!  Count number of changes in direction of monotonicity.
;!
    endif else begin
       if ( temp EQ 0.0D+00 ) then begin

          if ( del2 NE 0.0D+00 ) then begin
             ;if ( pchst ( dsave, del2 ) < 0.0D+00 ) then
             if ( (dsave*del2) LT 0.0D+00 ) then ierr = ierr + 1
             dsave = del2
          endif
;!
;!  Use Brodlie modification of Butland formula.
;!
       endif else begin  ;; TEMP GT 0.

          hsumt3 = 3.0D+00 * hsum
          w1 = ( hsum + h1 ) / hsumt3
          w2 = ( hsum + h2 ) / hsumt3
          dmax = max ([abs ( del1 ), abs ( del2 )])
          dmin = min ([abs ( del1 ), abs ( del2 )])
          drat1 = del1 / dmax
          drat2 = del2 / dmax
          d[i] = dmin / ( w1 * drat1 + w2 * drat2 )
          
       endelse
    endelse

 endfor
;!
;!  Set D(N) via non-centered three point formula, adjusted to be
;!  shape preserving.
;!
  w1 = -h2 / hsum
  w2 = ( h2 + hsum ) / hsum
  d[n-1] = w1 * del1 + w2 * del2
  
  ;if ( pchst ( d(n), del2 ) <= 0.0D+00 ) then
  if (  (d[n-1]*del2) LE 0.0D+00 ) then d[n-1] = 0.0D+00 $
  else begin
     if ( (del1*del2) LT 0.0D+00 ) then begin
                                ;else if ( pchst ( del1, del2 ) < 0.0D+00 ) then
;!
;!  Need do this check only if monotonicity switches.
;!
        dmax = 3.0D+00 * del2
        if ( abs ( dmax ) LT abs ( d[n-1] ) ) then d[n-1] = dmax
     endif
  endelse

  return, d
end
