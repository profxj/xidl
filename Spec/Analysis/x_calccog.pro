;+ 
; NAME:
; x_calccog
;   Version 1.1
;
; PURPOSE:
;  Perform a COG solution given a list of flambda values, EW
;  measurements and the errors in those measurements.  The user inputs
;  an array of column densities to search over and the code finds the
;  minimum in chi^2 and returns the value.
;
;  This code is currently only used to calculate the best b value for
;  the lines for the given N values.  See fuse_cog for a better COG
;  routine.  [JXP: 5/25/04]
;   
;
; CALLING SEQUENCE:
;   x_calccog, ew_red, sig_ew, flambda, Nval, blmt, FNDB=
;   BANS=, /EXACT
;
; INPUTS:
;   ew_red  -- Reduced ew array
;   sigew   -- Error in ew 
;   flambda  -- f lambda (assumes A not cm)
;   Nval    -- List of Nval to consider
;   blmt    --  Limits for b value search (assumes km/s)
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  /EXACT  -- Calculate the COG exactly  (otherwise use lookup table)
;  /FNDB   -- Find the b-value (only option now)
;
; OPTIONAL OUTPUTS:
;  BANS=  -- The b value which minimizes chi^2
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES/FUNCTIONS CALLED:
;  x_calccog_initcomm
;  x_calccog_cog
;  x_calccog_ftau
;
; REVISION HISTORY:
;   05-Mar-2002 Written by JXP
;-
;------------------------------------------------------------------------------

pro x_calccog_initcomm

common x_calccog_cmm, cog_bval, cog_tau, cog_strct

return
end

;;;
function x_calccog_cog, tau, EXACT=exact

common x_calccog_cmm

  ; Common
  cog_tau = tau

  ;; Integrate (EXACT)
  if keyword_set( EXACT ) then $
    ftau = qromo('x_calccog_ftau', 0., /double, /midexp) $
  else ftau = spl_interp(cog_strct.tau, cog_strct.ftau, cog_strct.splint, $
                         cog_tau)

  ; Other factors
  return, 2.d*cog_bval*ftau/(3.e5)
end
  
;;;;
function x_calccog_ftau, x

common x_calccog_cmm

  stp1 = -cog_tau * exp(-x^2)
  ; Integrand
  if stp1 LT -80.d then ftauint = 1. else $
    ftauint = (1.d - exp(stp1))
  return, ftauint
end
  

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro x_calccog, ew_red, sig_ew, flambda, Nval, blmt, $
               BANS=bans, FNDB=fndb, EXACT=exact, EWONLY=ewonly

  common x_calccog_cmm
;
  if  N_params() LT 5  then begin 
    print,'Syntax - ' + $
             'x_calccog, ew_red, sig_ew, flambda, Nval, blmt, '
    print, '       BANS=, /FNDB  [v1.1]'
    return
  endif 


  x_calccog_initcomm

; Optional Keywords

  if not keyword_set(NCHI) then nchi = 50L
  if not keyword_set(EXACT) then begin
      cog_fil = getenv('XIDL_DIR')+'/Spec/Analysis/cogmax_tab.fits'
      cog_strct = xmrdfits(cog_fil, 1, /silent)
  endif

  npt = n_elements(ew_red)

; BEGIN

  if keyword_set( FNDB ) then begin ; FINDB
      nN = n_elements(Nval)
      bans = fltarr(nN)
      chisq = fltarr(nchi)
      ; Loop on N values
      for qq=0L,nN-1 do begin
          ; Loop on b values
          for ii=0L, nchi-1 do begin
              cog_bval = blmt[0] + (blmt[1]-blmt[0])*float(ii)/float(nchi)
              calcEW = fltarr(npt)
              tau = 1.497e-2*(flambda*1.e-8)*(10^Nval[qq])/(cog_bval*1e5) 
              for jj=0L, npt-1 do $
                calcEW[jj] = x_calccog_cog(tau[jj], EXACT=exact)
              ; Calculate chisq
              chisq[ii] = total( ((calcEW-ew_red)/(sig_ew))^2 )
          endfor
          mn = min(chisq, imn)
          bans[qq] = blmt[0] + (blmt[1]-blmt[0])*float(imn)/float(nchi)
      endfor
  endif else return

end
          
                  
