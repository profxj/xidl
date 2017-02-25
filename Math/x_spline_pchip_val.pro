;+ 
; NAME: 
; x_spline_pchip_val
;    Version 1.0
;
; PURPOSE:
;!    This routine evaluates the cubic Hermite function at the points XE.
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

;!*****************************************************************************80
;!
;!! SPLINE_PCHIP_VAL evaluates a piecewise cubic Hermite function.
;!
;!  Description:
;!
;!    This routine may be used by itself for Hermite interpolation, or as an
;!    evaluator for SPLINE_PCHIP_SET.
;!
;!    This routine evaluates the cubic Hermite function at the points XE.
;!
;!    Most of the coding between the call to CHFEV and the end of
;!    the IR loop could be eliminated if it were permissible to
;!    assume that XE is ordered relative to X.
;!
;!    CHFEV does not assume that X1 is less than X2.  Thus, it would
;!    be possible to write a version of SPLINE_PCHIP_VAL that assumes a strictly
;!    decreasing X array by simply running the IR loop backwards
;!    and reversing the order of appropriate tests.
;!
;!    The present code has a minor bug, which I have decided is not
;!    worth the effort that would be required to fix it.
;!    If XE contains points in [X(N-1),X(N)], followed by points less than
;!    X(N-1), followed by points greater than X(N), the extrapolation points
;!    will be counted (at least) twice in the total returned in IERR.
;!
;!    The evaluation will be most efficient if the elements of XE are
;!    increasing relative to X; that is, for all J <= K,
;!      X(I) <= XE(J)
;!    implies
;!      X(I) <= XE(K).
;!
;!    If any of the XE are outside the interval [X(1),X(N)],
;!    values are extrapolated from the nearest extreme cubic,
;!    and a warning error is returned.
;!
;!    This routine was originally called "PCHFE".
;!
;!  Licensing:
;!
;!    This code is distributed under the GNU LGPL license.
;!
;!  Modified:
;!
;!    14 August 2005
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
;!  Parameters:
;!
;!    Input, integer ( kind = 4 ) N, the number of data points.  N must be
;!    at least 2.
;!
;!    Input, real ( kind = 8 ) X(N), the strictly increasing independent
;!    variable values.
;!
;!    Input, real ( kind = 8 ) F(N), the function values.
;!
;!    Input, real ( kind = 8 ) D(N), the derivative values.
;!
;!    Input, integer ( kind = 4 ) NE, the number of evaluation points.
;!
;!    Input, real ( kind = 8 ) XE(NE), points at which the function is to
;!    be evaluated.
;!
;!    Output, real ( kind = 8 ) FE(NE), the values of the cubic Hermite
;!    function at XE.
;!
;-
function x_spline_pchip_val, n, x, f, d, nex, xe

  compile_opt strictarr
  fe = dblarr(nex)
;!
;!  Check arguments.
;!
  if ( n LT 2 ) then begin
     ierr = -1
                                ;write ( *, '(a)' ) ' '
                                ;write ( *, '(a)' ) 'SPLINE_PCHIP_VAL - Fatal error!'
                                ;write ( *, '(a)' ) '  Number of data points less than 2.'
     stop
  endif

  for i = 1, n-1 do begin
     if ( x[i] LE x[i-1] ) then begin
        ierr = -3
      ;write ( *, '(a)' ) ' '
      ;write ( *, '(a)' ) 'SPLINE_PCHIP_VAL - Fatal error!'
      ;write ( *, '(a)' ) '  X array not strictly increasing.'
        stop
     endif
  endfor

  if ( nex LT 1 ) then begin
     ierr = -4
     ;write ( *, '(a)' ) ' '
    ;write ( *, '(a)' ) 'SPLINE_PCHIP_VAL - Fatal error!'
    ;write ( *, '(a)' ) '  Number of evaluation points less than 1.'
     return, -1
  endif

  ierr = 0
;!
;!  Loop over intervals.
;!  The interval index is IL = IR-1.
;!  The interval is X(IL) <= X < X(IR).
;!
  j_first = 0
  ir = 1

;  do
;!
;!  Skip out of the loop if have processed all evaluation points.
;!
;    if ( nex < j_first ) then
;      exit
;    end if
  flg = 1
  while(flg) do begin
     if ( nex LT j_first ) then begin
        flg = 0
        break
     endif
;!
;!  Locate all points in the interval.
;!
     j_save = nex 
     
     for j = j_first, nex-1 do begin
        if ( x[ir] LE xe[j] ) then begin
           j_save = j
           if ( ir EQ (n-1) ) then begin
              j_save = nex ;+ 1
           endif
           break                ;exit
        endif
     endfor
;!
;!  Have located first point beyond interval.
;!
     j = j_save
     nj = j - j_first
     ;if ir EQ 10 then stop
;!
;!  Skip evaluation if no points in interval.
;!
     if ( nj NE 0 ) then begin
;!
;!  Evaluate cubic at XE(J_FIRST:J-1).
;!
        tmp = fe[j_first:j-1]
        x_chfev, x[ir-1], x[ir], f[ir-1], f[ir], d[ir-1], d[ir], $
                 nj, xe[j_first:j-1], tmp, next, ierc 
        fe[j_first:j-1] = tmp
                                ;call chfev ( x(ir-1), x(ir), f(ir-1), f(ir), d(ir-1), d(ir), &
                                ;  nj, xe(j_first:j-1), fe(j_first:j-1), next, ierc )
        
        if ( ierc LT 0 ) then begin
           ierr = -5
       ;write ( *, '(a)' ) ' '
       ;write ( *, '(a)' ) 'SPLINE_PCHIP_VAL - Fatal error!'
       ;write ( *, '(a)' ) '  Error return from CHFEV.'
           stop
        endif
;!
;!  In the current set of XE points, there are NEXT(2) to the right of X(IR).
;!
        if ( next[1] NE 0 ) then begin
           
           if ( ir LT (n-1) ) then begin
              ierr = -5
          ;write ( *, '(a)' ) ' '
          ;write ( *, '(a)' ) 'SPLINE_PCHIP_VAL - Fatal error!'
          ;write ( *, '(a)' ) '  IR < N.'
              stop
           endif
;!
;!  These are actually extrapolation points.
;!
           ierr = ierr + next[1]
        endif
;!
;!  In the current set of XE points, there are NEXT(1) to the left of X(IR-1).
;!
        if ( next[0] NE 0 ) then begin
;!
;!  These are actually extrapolation points.
;!
           if ( ir LE 2 ) then begin
              ierr = ierr + next[0]
           endif else begin
              j_new = -1
              
              for i = j_first, j - 1 do begin
                 if ( xe[i] LT x[ir-1] ) then begin
                    j_new = i
                    break
                 endif
              endfor
              
              if ( j_new EQ (-1) ) then begin
                 ierr = -5
                ;write ( *, '(a)' ) ' '
                ;write ( *, '(a)' ) 'SPLINE_PCHIP_VAL - Fatal error!'
                ;write ( *, '(a)' ) '  Could not bracket the data point.'
                 stop
              endif
;!
;!  Reset J.  This will be the new J_FIRST.
;!
              j = j_new
;!
;!  Now find out how far to back up in the X array.
;!
              for i = 0, ir-2 do begin
                 if ( xe[j] LT x[i] ) then break       
              endfor
;!
;!  At this point, either XE(J) < X(1) or X(i-1) <= XE(J) < X(I) .
;!
;!  Reset IR, recognizing that it will be incremented before cycling.
;!
              ir = max ( 1, i-1 )
              
           endelse
        endif
        
        j_first = j
        
     endif
     
     ;if ir EQ 10 then stop
     ir = ir + 1
     
     if ( n LT ir ) then flg = 0
     
  endwhile
  
  return, fe
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Tests the monotonic spline
;pro x_tst_monospline
;
;  ;; Test values
;  x = findgen(11)
;  y = [10., 9.5, 9., 8.5, 3.0, 2.0, 1.5, 1.4, 1.2, 1.1, 1.0]
;
;  ;; Evaluation points
;  xev = -1. + findgen(100)/99.*12.
;
;  ;; Traditional spline
;  splin = spl_init(x, y, /double)
;  yev = spl_interp(x, y, splin, xev, /double)
;
;  ;; Mono spline
;  d = x_spline_pchip_set(11, x, y)
;  yev_mono = x_spline_pchip_val(11, x, y, d, 100L, xev)
;
;  ;; Plot
;  x_splot, x, y, xtwo=xev, ytwo=yev, xthr=xev, ythr=yev_mono, $
;           psym1=1, psym3=-3, /bloc
;
;  return
;end
