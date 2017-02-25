;+ 
; NAME:
; x_skyshift
;     Version 1.1
;
; PURPOSE:
;    Identify a shift in the sky spectrum between the one extracted 
;    for a given frame and the UVES archived spectrum.
;
; CALLING SEQUENCE:
;   
;  x_skyshift, bset_prof, uveslog, uvesflux, xwave, THRESH=thresh
;
; INPUTS:
;  bset_prof - The bspline info from the sky fit
;  uveslog   - alog10(sky wavelengths) from UVES line list
;  uvesflux  - Flux of the sky lines
;  xwave     - Wavelengths of the sky
;
; RETURNS:
;   The sky shift in fractional pixels
;
; OUTPUTS:
;   A sky subtracted image added as extension #2 to the fits file.
;   A file in the directory Sky/ describing the fit
;
;
; OPTIONAL KEYWORDS:
;   THRESH  -  Minimum threshold to include line in the analysis
;             (Default: 0.4)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;          sky_shift = x_skyshift(bset_prof, uveslog, uvesflux, $
;                                    objstr[mm].sky_wv[0:nrow-1] ) 
;
; PROCEDURES/FUNCTIONS CALLED:
;  ladfit
;  bspline_valu()
;
; REVISION HISTORY:
;   ??-2003     Written by SMB
;-
;------------------------------------------------------------------------------
function x_skyshift, bset_prof, uveslog, uvesflux, xwave, thresh=thresh

  if  N_params() LT 4  then begin 
      print,'Syntax - ' + $
        'shift = x_skyshift(bset_prof, uveslog, uvesflux, xwave, THRESH=)' + $
        ' [v1.1]'
      return, -1
  endif 

  if NOT ARG_PRESENT(thresh) then thresh=0.4

  nord = bset_prof.nord
  goodsky = where(bset_prof.bkmask GT 0 AND bset_prof.icoeff GT 0.2)
  
  maxwv = max(bset_prof.fullbkpt[nord+goodsky] , min=minwv)

  uvesin = where(uveslog LT maxwv AND uveslog GT minwv $
                 AND uvesflux GT thresh, nuvesin)
  
  if nuvesin EQ 0 then begin
      print, 'X_SKYSHIFT: No skylines found, no shift calculated'
      return, 0.
  endif
  
  xuse = where(xwave GT minwv AND xwave LT maxwv, nuse)
  
  dispcoeff = ladfit(0.5*(xwave[xuse]+xwave[xuse[1:*]]), $
                     xwave[xuse[1:*]]-xwave[xuse])
  
  uvesdisp = poly(uveslog[uvesin], dispcoeff)
  
  template = fltarr(nuse)
  for i=0, nuvesin-1 do begin
      diff = abs(xwave[xuse] - uveslog[uvesin[i]])/uvesdisp[i]
      diffin = where(diff LT 10.0, nin)
      if nin GT 10 then $
        template[diffin] = template[diffin] + $
        exp(-0.5*diff[diffin]^2) * uvesflux[uvesin[i]]
  endfor
  
  lag = lindgen(11)-5
  skyflux = bspline_valu(xwave[xuse], bset_prof)
  corr = c_correlate(skyflux, template, lag)
  
  totalcorr = total(corr)
  if totalcorr LT 0.5 then begin
      print, 'X_SKYSHIFT: found', nuvesin,' lines', format='(a,i4,a,$)'
      print, ', but low correlation: ', totalcorr, $
        format='(a,f8.3)'
      return, 0.0
  endif
  
  shift_pix = total(lag*corr)/totalcorr
  print, 'X_SKYSHIFT: found ', nuvesin, ' lines.  Correlation/Shift is ', $
    totalcorr, shift_pix, format='(a,i3,a,f8.3, f8.3)'
  
  return, shift_pix
end


     
