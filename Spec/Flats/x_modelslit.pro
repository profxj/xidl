;+ 
; NAME:
; x_modelslit
;     Version 1.1
;
; PURPOSE:
;  Model slit fits a full 2-d image of a flat-field frame (non-diffused)
;  It fits one row at a time, so the final image can have structure on a row by
;    row basis
;
; CALLING SEQUENCE:
;   
;  x_modelslit, x, setup, side, [/chk, /clobber]
;
; INPUTS:
;   x     -  MIKE structure
;   setup    -  Setup identifier 
;   [side]   -  Blue (1), Red (2), or both [1,2L]    (Default: [1,2L])
;
; RETURNS:
;
; OUTPUTS:
;   Fills in the profile0 and profile1 tags in the Order structure
;
; OPTIONAL KEYWORDS:
;  /CHK  - Show the profiles and fits order by order
;  /CLOBBER - Clobber an previous work
;  TFLAT_FIL - Name for TFLAT file
;  RESIDUAL_FIL  - Output name for Jacobian matrix
;  PROFILE_FIL - Output name for profile fits
;  DETILT    - Remove a linear tilt from the Cross-section fit
;
;  NXBKT  -  Number of x breakpoints for scattered light fit 
;  NYBKT  -  Number of x breakpoints for scattered light fit 
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   x_slitflat, mike, 1
;
; PROCEDURES/FUNCTIONS CALLED:
;  x_ordermask
;  x_fitgap
;  x_qw
;  x_slitflat_work
;
; REVISION HISTORY:
;   ??-??-2003 Written by SB
;-
;------------------------------------------------------------------------------
pro x_modelslit, tflat, tflativar, ordr_str, focus=focus, $
             scat_model=scat_model, ncoeff=ncoeff, terms=terms, $
             coeff_array=coeff_array

       if NOT keyword_set(ncoeff) then ncoeff = 8
       if NOT keyword_set(focus) then focus = 10.

       ncol = (size(tflat))[1]
       nrow = (size(tflat))[2]
       nord = n_elements(ordr_str)
       model  = tflat*0.0

       x = 2*findgen(ncol)/ncol - 1
       legarr = flegendre(x, ncoeff)
       ordercen =  0.5*(ordr_str.rhedg + ordr_str.lhedg)
       orderwidth = ordr_str.rhedg - ordr_str.lhedg
       xord = findgen(ncol) # replicate(1,nord)        
       r = replicate(1.0, ncol)

       coeff_array = fltarr(nord+ncoeff, nrow)
;       tflativar[420:421,*] = 0.
       

       for i=0L,nrow-1 do begin
          slit_frac = 2.0*(xord - r # (ordercen[i,*]))/(r#orderwidth[i,*])
          slit_basis = 0.5 - errorf(focus*(abs(slit_frac)-1))/2.

          med_width = median(orderwidth[i,*])
          mask = (xord[*,0] GE min(ordr_str.lhedg[i])-1 AND $
                  xord[*,0] LE max(ordr_str.rhedg[i])+1) * $
                  (tflativar[*,i] GT 0)
          good_slit = where(total(slit_basis * $
                  (mask # replicate(1,nord)),1) GT 0.5 * med_width, ngood)
          if ngood LT nord/2 then continue

          ;; Check on ivar
          gd2 = where(tflativar[*,i] GT 0., ngd2)
          if ngd2 LT ncol/2 then continue


          full_basis = [[slit_basis[*,good_slit]],[legarr]]

          res = 0
          for iiter = 1,2 do begin
              sqinv = sqrt(tflativar[*,i] * mask)
              if (total(sqinv GT 0) LT nord+ncoeff) then continue
              
              as =  full_basis * (sqinv # replicate(1,nord+ncoeff))
              alpha = transpose(as) # as
              beta = transpose(as) # (tflat[*,i] * sqinv)

              ;; Singular value? svsol
              svdc, alpha, w, u, v, /double
              res = svsol(u, w, v, beta, /double)
              yfit = full_basis # res
              
              diff = (tflat[*,i] - yfit) * sqinv
              mask = mask * (abs(diff) LT 25)
          endfor

          if keyword_set(res) then begin
           coeff_array[good_slit,i] = res[0:ngood-1]
           coeff_array[nord:*,i] = res[ngood:*]
           model[*,i] = yfit
         endif
   
      endfor

      scat_model = legarr # coeff_array[nord:*,*]

return
end 


