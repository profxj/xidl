;+ 
; NAME:
; x_extract   
;   Version 1.1
;
; PURPOSE:
;    Extracts a 1D spectrum from a 2D image given an aperture.  This
;  routine is superseded by x_extobjbox
;
; CALLING SEQUENCE:
;   spec = x_extract(img, aper, trace, [error, sky], CAPER=, GAIN=,
;                    RN=, /OPTIMAL)
;
; INPUTS:
;   img       - 2D Image
;   aper      - Aperture to extract
;   trace     - Trace of the spectrum
;   [sky]     - Optional input required for error output
;
; RETURNS:
;   spec      - 1D Spectrum
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   CAPER      - Line chosen for aperture: Require with TRACE
;   GAIN       - CCD gain
;   RN       - Read noise
;
; OPTIONAL OUTPUTS:
;   var      - Error array
;
; COMMENTS:
;
; EXAMPLES:
;   spec = x_extract(img, aper)
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   05-Dec-2001 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function x_extract, img, aper, trace, svar, sky, CAPER=caper, GAIN=gain, $
                    RN=rn, OPTIMAL=optimal, SKYRMS=skyrms, VAR=var

;
  if  N_params() LT 3  then begin 
    print,'Syntax - ' + $
             'spec = x_extract(img, aper, trace, [svar, sky], CAPER=, GAIN=,'
    print, '        RN=, /OPTIMAL, SKYRMS=, VAR=var) [V1.0]'
    return, -1
  endif 

  sz = size(img, /dimensions)
  npix = sz[0] < n_elements(trace)

;  Optional Keywords

  if keyword_set( TRACE ) AND not keyword_set( CAPER ) then begin
      print, 'Assuming aperture comes from the center line'
      caper = (sz[0]-1)/2
  endif

  if not keyword_set( RN ) then rn = 5.
  if not keyword_set( GAIN ) then gain = 1.
  if arg_present( svar ) then begin
      if not keyword_set( SKY ) then begin
          print, 'x_extract: Sky not set -- taking sky = 0. everywhere'
          sky = fltarr(sz[0],sz[1])
      endif
      if not keyword_set( VAR ) then stop
      if not keyword_set( SKYRMS ) then begin
          print, 'x_extract: SKYRMS not set -- taking skyrms = 0. everywhere'
          skyrms = fltarr(sz[0])
      endif
      svar = fltarr(npix)
  endif

; Define spec

  spec = fltarr(npix)

; EXTRACT

  if not keyword_set( OPTIMAL ) then begin  ; BOXCAR
      center = trace[long(caper)]
      for q=0L, npix-1 do begin
          offset = trace[q] - center
          spec[q] = gain*$
            total( img[q,round(aper[0]+offset):round(aper[1]+offset)], 2)
          if arg_present(svar) then begin
              svar[q] = $
                total( var[q,round(aper[0]+offset):round(aper[1]+offset)], 2)
          endif
      endfor
  endif else begin ; OPTIMAL
      print, 'Not set for optimal!'
      return, -1
  endelse

; RETURN

  return, spec
end
