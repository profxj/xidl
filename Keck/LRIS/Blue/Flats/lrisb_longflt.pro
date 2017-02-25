;+ 
; NAME:
; lrisb_longflt   
;     Version 1.1
;
; PURPOSE:
;    Read in a list of flat files, bias subtract them, combine
;  and then normalize by fitting a bspline down the middle.
;
; CALLING SEQUENCE:
;  lrisb_longflt, flatfil, nrmflat, FITS=, GAIN=, /NOFIT, /FULL
;
; INPUTS:
;   flatfil -- List of flats
;
; RETURNS:
;   nrmflat -- Normalized flat image
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  /NOFIT -- Do not fit a bspline to the median.
;          (Recommended for the 1200 grism)
;  GAIN=  -- Gain [default: 1.6]
;  /FULL  -- Keyword to lrisb_subbias
;
; OPTIONAL OUTPUTS:
;  FITS=  -- File to write nrmflat to 
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
                   NOFIT=nofit, FULL=full, CHK=chk

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
      spos = strpos(flats[i], 'w/')
      if spos EQ -1 then spos = 0 else spos = spos + 2
      ;; Outfil
      outnm = 'OV/ov_'+strmid(flats[i],spos)
      ;; Subtract
      lrisb_subbias, raw, flat, /LONG, FITS=outnm, FULL=full
      if i NE 0 then all_out = [all_out, outnm] else all_out = [outnm]
  endfor

; COMBINE
  xcombine, all_out, cmbflat, FCOMB=2, SCALE='MED', gain=gain, RN=3.

; Delete OV image
  for i=0L,nflt-1 do spawn, '\rm '+'ov_'+flats[i]
      
  sz = size(cmbflat,/dimensions)

; Now normalize the whole flat

  medflt = djs_median(cmbflat[sz[0]/2 - 15: sz[0]/2 + 15,*],1)
  if not keyword_set( NOFIT ) then $
    bset = bspline_iterfit(findgen(n_elements(medflt)),medflt, yfit=yfit,$
                           everyn=15) $
  else yfit = medflt
  if keyword_set(CHK) then x_splot, medflt, ytwo=yfit, /bloc
  nrmflat = cmbflat / ( replicate(1., sz[0]) # yfit)
  
; Output
  if keyword_set( FITS ) then mwrfits, nrmflat, fits, /create

  if not keyword_set( SVOV ) then begin
      for i=0L,nflt-1 do begin
          spawn, '\rm '+all_out[i]
      endfor
  endif
;  Optional Keywords
  print, 'lrisb_longflt: All done!'

  return
end
