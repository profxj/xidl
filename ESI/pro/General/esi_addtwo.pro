;+ 
; NAME:
; esi_addtwo
;   Version 1.0
;
; PURPOSE:
;    Combines two flats, rejecting Cosmic Rays
;
; CALLING SEQUENCE:
;   
;   img = esi_addtwo(esi, indx, VAR=var)
;
; INPUTS:
;   esi
;   indx
;
; RETURNS:
;   img       - Combine image
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;   VAR       - Variance
;
; COMMENTS:
;
; EXAMPLES:
;   img = esi_addtwo(esi, indx)
;
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   28-Aug-2002 Written by JXP
;-
;------------------------------------------------------------------------------

function esi_addtwo, esi, indx, VAR=var, SCALE=scale, SKY=sky, $
                     SATUR=satur, ARC=arc

;
  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
             'img = esi_addtwo(esi, indx, VAR=, /SCALE) [v1.0]'
    return, -1
  endif 


;  Optional Keywords

; Allow img to be fits file or data

  if keyword_set( SKY ) then begin
      dat1 = xmrdfits(esi[indx[0]].img_final, 2, /silent)
      dat2 = xmrdfits(esi[indx[1]].img_final, 2, /silent)
  endif else begin
      dat1 = xmrdfits(esi[indx[0]].img_final, /silent)
      dat2 = xmrdfits(esi[indx[1]].img_final, /silent)
  endelse

; Variance
  if arg_present(VAR) or keyword_set(SKY) then begin
      var1 = xmrdfits(esi[indx[0]].img_final, 1, /silent)
      var2 = xmrdfits(esi[indx[1]].img_final, 1, /silent)
  endif

  if keyword_set( SCALE ) or keyword_set( SKY ) then begin
      dat2 = dat2 * esi[indx[0]].exp / esi[indx[1]].exp
      var2 = var2 * esi[indx[0]].exp / esi[indx[1]].exp
      scl = esi[indx[0]].exp / esi[indx[1]].exp
  endif else scl = 1.

  if keyword_set( SKY ) then begin

      ;; Mask
      sz = size(dat1, /dimensions)
      msk = bytarr(sz[0],sz[1]) + 1B
      ;; Avoid var<0
      bad1 = where(var1 LE 0. OR dat1 EQ 0.)
      bad2 = where(var2 LE 0. OR dat2 EQ 0.)
      msk[bad1] = 0B
      msk[bad2] = 0B
      gd = where(msk NE 0B)
      ;; Subtract
      sub = dat1[gd] - dat2[gd]

      ;; Stats on the ratio
      djs_iterstat, sub, sigrej=3.0, median=med_rtio, sigma=sig_rtio, maxiter=2

      ;; Find all bad pixels
      bdpix = where(abs(sub-med_rtio) GT 10.*sig_rtio, nbad)
      bdpix = gd[bdpix]

      ;; Average (not add)
      fimg = (dat1+dat2)/2.
      var = (var1+var2)/2.
      
      ;; Take minimum of bad pixels
      if nbad GT 0 then begin
          fimg[bdpix] = dat1[bdpix] < dat2[bdpix]
          if arg_present(VAR) then var[bdpix] = ( var1[bdpix] < var2[bdpix] )
      endif
      
      return, fimg
  endif

  if keyword_set( ARC ) then begin

      if not keyword_set( SATUR ) then satur = 56000
      ;; Mask
      sz = size(dat1, /dimensions)
      msk = bytarr(sz[0],sz[1]) + 1B
      ;; Avoid var<0
      bad1 = where(dat1 EQ 0., nbd1)
      bad2 = where(dat2 EQ 0., nbd2)
      if nbd1 NE 0 then msk[bad1] = 0B
      if nbd2 NE 0 then msk[bad2] = 0B

      gd = where(msk NE 0B)
      ;; Subtract
      sub = dat1[gd] - dat2[gd]

      ;; Stats on the ratio
      djs_iterstat, sub, sigrej=3.0, median=med_rtio, sigma=sig_rtio, maxiter=2

      ;; Find all bad pixels
      bdpix = where(abs(sub-med_rtio) GT 10.*sig_rtio, nbad)
      bdpix = gd[bdpix]

      ;; Average (not add)
      fimg = (dat1+dat2)/2.
      var = (var1+var2)/2.
      
      ;; Take minimum of bad pixels
      if nbad GT 0 then begin
          fimg[bdpix] = dat1[bdpix] < dat2[bdpix]
          if arg_present(VAR) then var[bdpix] = ( var1[bdpix] < var2[bdpix] )
      endif
      
      return, fimg
  endif

end
