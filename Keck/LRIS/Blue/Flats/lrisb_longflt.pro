;+ 
; NAME:
; lrisb_longflt   
;     Version 1.1
;
; PURPOSE:
;    Median combine all ZRO frames (darks)
;
; CALLING SEQUENCE:
;   
;  lrisb_subbias, lris, indx
;
; INPUTS:
;   lrisb   -  ESI structure
;   indx  -  Index numbers of frame to subtract (default output is OV)
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  /NOFIT -- Recommended for the 1200 grism
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;  Currently only good for 1x1 binning
;
; EXAMPLES:
;   lrisb_subbias, lris, [47L,48L,49L]
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   30-Apr-2003 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro lrisb_longflt, flatfil, nrmflat, FITS=fits, GAIN=gain, $
                   NOFIT=nofit

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'lrisb_longflt, img, newflt, FITS= '
      return
  endif 

; Optional keywords
  if not keyword_set( GAIN ) then gain = 1.6

  readcol, flatfil, flats, FORMAT='A'

  nflt = n_elements(flats)


; Bias subtract

;  raw = xmrdfits(flatfil,/fscale)
  for i=0L,nflt-1 do begin
      ;; 
      raw = xmrdfits(flats[i],/fscale)
      ;; Outfil
      outnm = 'ov_'+flats[i]
      ;; Subtract
      lrisb_subbias, raw, flat, /LONG, FITS=outnm
  endfor

; COMBINE
  xcombine, 'ov_'+flats, cmbflat, FCOMB=2, SCALE='MED', gain=gain, RN=3.

; Delete OV image
  for i=0L,nflt-1 do spawn, '\rm '+'ov_'+flats[i]
      
  sz = size(cmbflat,/dimensions)

; Now normalize the whole flat

  medflt = djs_median(cmbflat[sz[0]/2 - 15: sz[0]/2 + 15,*],1)
  if not keyword_set( NOFIT ) then $
    bset = bspline_iterfit(findgen(n_elements(medflt)),medflt, yfit=yfit,$
                           everyn=15) $
  else yfit = medflt
  nrmflat = cmbflat / ( replicate(1., sz[0]) # yfit)
  
; Output
  if keyword_set( FITS ) then mwrfits, nrmflat, fits, /create

;  Optional Keywords
  print, 'lrisb_longflt: All done!'

  return
end
