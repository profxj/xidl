;+ 
; NAME:
; lrisb_subbias   
;     Version 1.1
;
; PURPOSE:
;    Subtract the overscan region from an LRIS image (blue side)
;
; CALLING SEQUENCE:
;  lrisb_subbias, img, ovimg, /LONG, FITS=, /FULL
;
; INPUTS:
;   img   -- Name of image to ov subtract
;   ovimg -- OV subtracted image
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  /FULL  -- Full readout
;
; OPTIONAL OUTPUTS:
;  FITS  -- Name of fits file to write finimg to
;
; COMMENTS:
;  Currently only good for 1x1 binning
;
; EXAMPLES:
;   lrisb_subbias, 'lblue1010.fits', ovimg
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   30-Apr-2003 Written by JXP
;-
;------------------------------------------------------------------------------
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro lrisr_subbias, img, ovimg, FITS=fits, FULL=full, ONEAMP=oneamp
;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'lrisb_subbias, img, ovimg, /LONG, /FULL, FITS=, [v1.1]'
      return
  endif 
  
;  Optional Keywords
;  if not keyword_set( OVROOT ) then ovroot = 'OV/ov_'
;  if not keyword_set( BIASFIL ) then begin
;      if lrisb[indx[0]].mode EQ 0 then biasfil = 'Bias/BiasI.fits' $
;        else biasfil = 'Bias/BiasS.fits'
;  endif

  

  if keyword_set( FULL ) then begin
      y1 = 615
      y2 = 1445
      x1_n = 40
      x1_x = 1063L
      x2_n = 1064L
      x2_x = 2087L
      ovimg = fltarr(2048L,y2-y1+1)
      ;; Amp1
      fitstr = x_setfitstrct()
      ov1 = djs_median(img[2090:2158L,*],1)
      f1 = x_fitrej(findgen(2048L), ov1, fitstr=fitstr)
      ovimg[0:x1_x-x1_n,*] = $
        img[x1_n:x1_x,y1:y2] -  replicate(1.,x1_x-x1_n+1)#f1[y1:y2]
      ;; Amp2
      fitstr = x_setfitstrct()
      ov2 = djs_median(img[2173:2240L,*],1)
      f2 = x_fitrej(findgen(2048L), ov2, fitstr=fitstr)
      ovimg[x1_x-x1_n+1:*,*] = $
        img[x2_n:x2_x,y1:y2] -  replicate(1.,x2_x-x2_n+1)#f2[y1:y2]
      ;; FITS
      if keyword_set( FITS ) then $
        mwrfits, ovimg, fits, /create, /silent
  endif else begin
      sz = size(img,/dimension)
      y1 = 140
      y2 = 900L < (sz[1]-1)
      x1_n = 42
      x1_x = 1063L
      x2_n = 1064L
      x2_x = 2085L
      nxtot = x2_x - x1_n + 1
      ovimg = fltarr(nxtot,y2-y1+1)
      ;; Amp1
      fitstr = x_setfitstrct()
      ov1 = djs_median(img[2090:2158L,*],1)
      f1 = x_fitrej(findgen(nxtot), ov1, fitstr=fitstr)
      ovimg[0:x1_x-x1_n,*] = $
        img[x1_n:x1_x,y1:y2] -  replicate(1.,x1_x-x1_n+1)#f1[y1:y2]
      ;; Amp2
      fitstr = x_setfitstrct()
      ov2 = djs_median(img[2173:2240L,*],1)
      f2 = x_fitrej(findgen(nxtot), ov2, fitstr=fitstr)
      ovimg[x1_x-x1_n+1:*,*] = $
        img[x2_n:x2_x,y1:y2] -  replicate(1.,x2_x-x2_n+1)#f2[y1:y2]
      ;; FITS
      if keyword_set( FITS ) then $
        mwrfits, ovimg, fits, /create, /silent
  endelse

  print, 'lrisb_subbias: All done!'

  return
end
