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
;  /LONG  -- Long slit mode
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

pro lrisb_subbias, img, ovimg, LONG=long, FITS=fits, FULL=full
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

  
; Open Bias file
;  print, 'lrisb_subbias: Using BIAS file: ', biasfil
;  if x_chkfil(biasfil+'*') EQ 0 then begin
;      print, 'lrisb_subbias: Bias file not found.  Create first!', biasfil
;      stop
;  endif
;  bias = xmrdfits(biasfil, /silent)
;  sz = size(bias, /dimensions)

; Loop

;  for q=0,n_elements(indx)-1 do begin
      ;; Check for output
;      outfil = ovroot+lrisb[indx[q]].img_root
;      a = findfile(outfil+'*', count=na)
;      if na NE 0 and not keyword_set( FORCE ) then begin
;          print, 'lrisb_subbias: File ', outfil, ' found. Not resubtracting'
;          lrisb[indx[q]].img_ov = outfil
;          lrisb[indx[q]].flg_ov = 1
;          continue
;      endif
      ;; Open Raw image
;      raw = xmrdfits(lrisb[indx[q]].rootpth+lris[indx[q]].img_root, 0, head, $
;                    /silent, /fscale)
;      print, 'lrisb_subbias: Subtracting the bias for ', $
;        lrisb[indx[q]].rootpth+lris[indx[q]].img_root

      ;; Allow for mode
  sz = size(img,/dimens)
  if keyword_set( LONG ) and not keyword_set( FULL ) then begin
      ovimg = fltarr(297,sz[1])
      ;; Amp3
      fitstr = x_setfitstrct()
      ov3 = djs_median(img[1081:1110,*],1)
      f3 = x_fitrej(findgen(sz[1]), ov3, fitstr=fitstr)
      ovimg = img[52:348,*] -  replicate(1.,297L)#f3
      ;; FITS
      if keyword_set( FITS ) then $
        mwrfits, ovimg, fits, /create, /silent
  endif

  if keyword_set( LONG ) and keyword_set( FULL ) then begin
      ovimg = fltarr(421,4096L)
      ;; Amp3
      fitstr = x_setfitstrct()
      ov3 = djs_median(img[4470L:4533L,*],1)
      f3 = x_fitrej(findgen(4096L), ov3, fitstr=fitstr)
      ovimg = img[2310L:2730L,*] -  replicate(1.,421)#f3
      ;; FITS
      if keyword_set( FITS ) then $
        mwrfits, ovimg, fits, /create, /silent
  endif

  print, 'lrisb_subbias: All done!'

  return
end
