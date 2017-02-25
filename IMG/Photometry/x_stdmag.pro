;+ 
; NAME:
; x_stdmag   
;   Version 1.1
;
; PURPOSE:
;    Does aperture photometry on the objects in the std_* list files
;
; CALLING SEQUENCE:
;   xdimg_stdmag, struct, SKYR=
;
; INPUTS:
;   img   -- Image
;  strfil -- Star list file
;  outfil -- Output file
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  SKYR =  2 float array describing sky annulus in arcsec
;            (default=15,30)
;  ARCPIX= Arcseconds per pixel
;  SATUR=  Saturation level
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   x_stdmag, 'Final/f_ccd003.fits', 'Photometry/star.list', 'outfil'
;
; PROCEDURES/FUNCTIONS CALLED:
;  X_APER
;  X_CCDINF
;
; REVISION HISTORY:
;   03-Apr-2003 Written by JXP
;-
;------------------------------------------------------------------------------
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro x_stdmag, img, strfil, outfil, MSK=msk, SKYR=skyr, ARCPIX=arcpix, $
              SATUR=satur

;
  if  N_params() LT 1  then begin 
      print, 'Syntax - ' +$
        'x_stdmag, img, strfil, outfil, MSK=, SKYR=, ARCPIX=, SATUR= (v1.1)'
      return
  endif 

  close, /all

;  Optional Keywords

  if not keyword_set( SKYR ) then SKYR = [15., 30.]
  if not keyword_set( ARCPIX ) then arcpix = 0.2
  if not keyword_set( SATUR ) then satur = 25000.

;  Set pixel scale, aperture

  x_ccdinf, struct[stds[0]].ccd, struct[stds[0]].tel, arcpix
  apr = [7.0/arcpix]
  skyrad = [12.0, 30.0]
  skyrad = skyrad/arcpix

; Set files
  datfil = ' '
  listfil = ' '
  nm = ''

  ;; Star file
  a = findfile(strfil, count=cnt)
  if cnt EQ 0 then continue

  ;; Data
  data = xmrdfits(img, /fscale, /silent)

  ;; Mask
  if keyword_set( MSK ) then $
    dmsk = xmrdfits(msk, /fscale, /silent)

  ;; Bad pixels (saturation)
  badpix = [-30000., satur]

  nlin = numlines(strfil)
  nam = strarr(nlin-2)
  xc = fltarr(nlin-2)
  yc = fltarr(nlin-2)
  
  openr, 1, strfil
  readf, 1, datfil
  readf, 1, listfil
  for j=0,nlin-3 do begin
      readf, 1, FORMAT='(a15,f,f)', nm, x, y
      nam[j] = nm
      xc[j] = x
      yc[j] = y
  endfor
          
  close, 1

; Aperture Photometry
  x_aper, data, xc, yc, mags, errap, sky, skyerr, struct[stds[q]].gain, $
    apr, skyrad, badpix, /silent
          
; Deal with exposure time [After to get correct errors]
          
  for j=0,nlin-3 do mags[j] = mags[j] + 2.5*alog10(struct[stds[q]].exp)
      
; Output
  openw, 1, outfil
  printf, 1, struct[stds[q]].filter
  printf, 1, listfil
  printf, 1, struct[stds[q]].AM
  
  for j=0,nlin-3 do printf, 1, FORMAT='(a15,f,f)', nam[j], mags[j], errap[j]
  close, 1
      
  print, 'All done with photometry!'
  
end
