;+ 
; NAME:
; dla_calcabnd
;  V1.1
;
; PURPOSE:
;    Calculates the abundance relative to the Sun, [X/Y]
;
; CALLING SEQUENCE:
;   dla_calcabnd, dla, nn, X, Y, ans, sig
;
; INPUTS:
;   dla -- DLA structure array
;   nn  -- Index of the structure
;   X   -- Atomic number of first element
;   Y   -- Atomic number of second element
;
; RETURNS:
;
; OUTPUTS:
;  ans --  [X/Y]
;  sig --  error in [X/Y]
;
; OPTIONAL KEYWORDS:
;  /NOSIGY -- Do not include error in Y in the calculation
;     This is only useful when dealing with [X/H]
;  /READ -- Use this when not running from dla_updabd. The reason
;           is that dla_updabd uses non-Log values for X and Y, 
;           and then saves them as Log values. When rerunning this
;           program, you take a log of a log. This flag fixes that. 
;   /ONE -- When called, this ignores Y, and just returns [X] relative
;           to solar.
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   dla_allabd, dla
;
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   06-Oct-2004 Written by JXP
;   05-Aug-2011 Added /READ and /ONE option in case calling after the fact MR & MN
;- 
;------------------------------------------------------------------------------
pro dla_calcabnd, dla, nn, X, Y, ans, sig, NOSIGY=nosigy, READ=read, ONE=one

  if (N_params() LT 4) then begin 
    print,'Syntax - ' + $
             'dla_calcabnd, dla, nn, X, Y, ans, sig, /NOSIGY [v1.1]'
    return
  endif 

  ;; Get abundances
  getabnd, nm, X, Xabnd, flag=1
  getabnd, nm, Y, Yabnd, flag=1

  ;; Calculate [X/Y]

  if not keyword_set(READ) then begin
     if X NE 1 then begin
        logX = alog10(dla[nn].elm[X].clm)
        logXsig = sqrt(((1./(alog(10.0)*dla[nn].elm[X].clm))^2)* $
                       dla[nn].elm[X].sigclm^2)
     endif else begin
        logX = dla[nn].elm[X].clm
        logXsig = dla[nn].elm[X].sigclm
     endelse
     if Y NE 1 then begin
        logY = alog10(dla[nn].elm[Y].clm)
        logYsig = sqrt(((1./(alog(10.0)*dla[nn].elm[Y].clm))^2)* $
                       dla[nn].elm[Y].sigclm^2)
     endif else begin
        logY = dla[nn].elm[Y].clm
        logYsig = dla[nn].elm[Y].sigclm
     endelse
  endif else begin
     logX = dla[nn].elm[X].clm
     logXsig = dla[nn].elm[X].sigclm
     if Y NE 1 then begin
        logY = dla[nn].elm[Y].clm
        logYsig = dla[nn].elm[Y].sigclm
     endif else begin
        logY = dla[nn].nhi
        logYsig = dla[nn].signhi
     endelse
  endelse

  if not keyword_set(ONE) then begin
     ans = logX - logY - Xabnd + Yabnd
  ;;     Error  
     if keyword_set(NOSIGY) then logYsig = 0.
     sig = sqrt(logXsig^2 + logYsig^2)
  endif else begin
     ans = logX - Xabnd
     sig = logXsig
  endelse




  return
end

