;+ 
; NAME:
; x_splncogmax
;   Version 1.0
;
; PURPOSE:
;   
;
; CALLING SEQUENCE:
;   
;   x_splncogmax, ew_red, sig_ew, flambda, Nval, blmt, FNDB=
;
; INPUTS:
;   ew_red  -- Reduced ew
;   sigew  - Error in ew
;   flambda  -- f lambda (assumes A not cm)
;   blmt    --  Limits for b value search (assumes km/s)
;
; RETURNS:
;   ew   - Equivalent width
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
; REVISION HISTORY:
;   15-Sep-2003 Written by JXP
;-
;------------------------------------------------------------------------------

pro x_splncogmax_initcomm

common x_splncogmax_cmm, cog_tau

return
end

;;;
function x_splncogmax_cog, tau

common x_splncogmax_cmm

  ; Common
  cog_tau = tau

  ; Integrate
  ftau = qromo('x_splncogmax_ftau', 0., /double, /midexp)

  ; Other factors
  return, ftau
end
  
;;;;
function x_splncogmax_ftau, x

common x_splncogmax_cmm

  stp1 = -cog_tau * exp(-x^2)
  ; Integrand
  if stp1 LT -80.d then ftauint = 1. else $
    ftauint = (1.d - exp(stp1))
  return, ftauint
end
  

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro x_splncogmax, tmax, tval, tabul, NSTP=nstp, STRCT=strct, $
                  OUTFIL=outfil

  common x_splncogmax_cmm
;
  if  N_params() LT 3  then begin 
    print,'Syntax - ' + $
             'x_splncogmax, tmax, tval, tabul, NSTP= [v1.0]'
    return
  endif 

  x_splncogmax_initcomm

; Optional Keywords

  if not keyword_set(NSTP) then nstp = 1000L
  if not keyword_set(TMIN) then tmin = 1e-4

  tabul = dblarr(nstp) 
  tval = dblarr(nstp) 
  ltmn = alog10(tmin)
  ltmx = alog10(tmax)

; LOOP

  for q=1L,nstp-1 do begin
      expon = ltmn + (ltmx-ltmn)*float(q)/float(nstp-1)
      cog_tau = 10^expon
      tval[q] = cog_tau
      ;; Integrate
      tabul[q] = x_splncogmax_cog(cog_tau)
  endfor

  splint = spl_init(tval, tabul, /double)

  ;; Spline
  xpnt = dblarr(10000L)
  ypnt = dblarr(10000L)
  for q=1L, 10000L-1 do begin
      xpnt[q] = 10^( ltmn + (ltmx-ltmn)*float(q)/float(10000.-1) )
      ypnt[q] = spl_interp(tval, tabul, splint, xpnt[q])
  endfor

  ;; Plot
  clr = getcolor(/load)
  x_splot, xpnt, ypnt, xtwo=tval, ytwo=tabul, /block, psym2=1
  

  strct = { $
            tau: tval, $
            ftau: tabul, $
            splint: splint $
          }
  return

  if keyword_set( OUTFIL ) then mwrfits, strct, outfil, /create

end

 
