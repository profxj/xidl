;+ 
; NAME:
; mike_addtwo
;   Version 1.0
;
; PURPOSE:
;    Combines two flats, rejecting Cosmic Rays
;
; CALLING SEQUENCE:
;   
;   img = mike_addtwo(mike, indx, VAR=var)
;
; INPUTS:
;   mike
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
;   IVAR       - Variance
;
; COMMENTS:
;
; EXAMPLES:
;   img = mike_addtwo(mike, indx)
;
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   30-Apr-2002 Written by JXP
;-
;------------------------------------------------------------------------------

function mike_addtwo, mike, indx, IVAR=ivar, SCALE=scale, SKY=sky, $
                     SATUR=satur, ARC=arc

;
  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
             'img = mike_addtwo(mike, indx, IVAR=, /SCALE) [v1.0]'
    return, -1
  endif 


;  Optional Keywords

; Allow img to be fits file or data

  if keyword_set( SKY ) then begin
      dat1 = xmrdfits(mike[indx[0]].img_final, 2, /silent)
      dat2 = xmrdfits(mike[indx[1]].img_final, 2, /silent)
  endif else begin
      dat1 = xmrdfits(mike[indx[0]].img_final, /silent)
      dat2 = xmrdfits(mike[indx[1]].img_final, /silent)
  endelse

; Variance
  if arg_present(IVAR) or keyword_set(SKY) then begin
      ivar1 = xmrdfits(mike[indx[0]].img_final, 1, /silent)
      ivar2 = xmrdfits(mike[indx[1]].img_final, 1, /silent)
  endif

  if keyword_set( SCALE ) or keyword_set( SKY ) then begin
      dat2 = dat2 * float(mike[indx[0]].exp / mike[indx[1]].exp)
      ivar2 = ivar2 * float(mike[indx[1]].exp / mike[indx[0]].exp)
      scl = float(mike[indx[0]].exp / mike[indx[1]].exp)
  endif else scl = 1.

  if keyword_set( SKY ) then begin

      ;; Mask
      sz = size(dat1, /dimensions)
      msk = bytarr(sz[0],sz[1]) + 1B
      ;; Avoid var<0
      bad1 = where(ivar1 LE 0. OR dat1 EQ 0.)
      bad2 = where(ivar2 LE 0. OR dat2 EQ 0.)
      msk[bad1] = 0B
      msk[bad2] = 0B
      gd = where(msk NE 0B)
      ;; Subtract
      sub = dat1[gd] - dat2[gd]

      ;; Stats on the ratio
      djs_iterstat, sub, sigrej=3.0, median=med_rtio, sigma=sig_rtio, maxiter=5

      ;; Find all bad pixels
      bdpix = where(abs(sub-med_rtio) GT 5.*sig_rtio, nbad)
      bdpix = gd[bdpix]

      ;; Average (not add)
      fimg = (dat1+dat2)/2.
      ivar = 2./(1./ivar1+1./ivar2)
      
      ;; Take minimum of bad pixels
      if nbad GT 0 then begin
          fimg[bdpix] = dat1[bdpix] < dat2[bdpix]
          ivar[bdpix] = ( ivar1[bdpix] < ivar2[bdpix] )
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
      ivar = 2./(1./ivar1+1./ivar2)
      
      ;; Take minimum of bad pixels
      if nbad GT 0 then begin
          fimg[bdpix] = dat1[bdpix] < dat2[bdpix]
          if arg_present(IVAR) then ivar[bdpix] = ( ivar1[bdpix] < ivar2[bdpix] )
      endif
      
      return, fimg
  endif


;;;;;;;;;

  ;; Ratio
  sub = dat1-dat2

  ;; Stats on the ratio
  djs_iterstat, sub, sigrej=3.0, median=med_rtio, sigma=sig_rtio, maxiter=5

  ; Find all bad pixels
  bdpix = where(abs(sub-med_rtio) GT 5.*sig_rtio, nbad)

  ; Average (not add)
  fimg = (dat1+dat2)/2.
  ivar = 2./ (1./ivar1 + 1./ivar2)

  ; Take minimum of bad pixels
  if nbad GT 0 then begin
      fimg[bdpix] = dat1[bdpix] < dat2[bdpix]
      ivar[bdpix] = ivar1[bdpix] < ivar2[bdpix] 
  endif

  return, fimg

end
