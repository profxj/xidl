;+ 
; NAME:
; x_rebin2dspec
;   Version 1.1
;
; PURPOSE:
;    Rebin a single data set to a new wavlength scale
;      Simple adding (no weighting by S/N)
;
; CALLING SEQUENCE:
;  x_rebin2dspec, wv, fx, nwwv, nwfx, VAR=, NWVAR=, GDPIX=,
;                  /SILENT, /REDBLUE, /CR, /REBINC
;
; INPUTS:
;   wv   -- Wavelength array
;   fx   -- Flux array
;   newwv -- New wavelength array to rebin to
;
; RETURNS:
;
; OUTPUTS:
;   newfx -- New flux array
;
; OPTIONAL KEYWORDS:
;  VAR=     -- Data is float
;  /REBINC  -- Rebin using a C program
;  /REDBLUE -- Data runs from red to blue
;
; OPTIONAL OUTPUTS:
;  NEWVAR=  - New variance array
;  /CR      - VAR = -1 flag CR and eliminate the pix
;
; COMMENTS:
;
; EXAMPLES:
;   x_rebin2dspec, gdpix, orig_wv, orig_fx, newwv, newfx
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   12-Apr-2002 Written by JXP
;   24-Aug-2002 Added C program
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro x_rebin2dspec, wv, fx, nwwv, nwfx, VAR=var, NWVAR=nwvar, GDPIX=gdpix,$
                   SILENT=silent, REDBLUE=redblue, CR=cr, REBINC=rebinc

;
  if  N_params() LT 4  then begin 
    print,'Syntax - ' + $
      'x_rebin2dspec, wv, fx, nwwv, nwfx, VAR=, NWVAR=, ' + $
      'GDPIX=, /CR, /REBINC [V1.1]'
    return
  endif 


;  Optional Keywords
  sz = size(wv, /dimensions)
  if not keyword_set(VAR) then var = fltarr(sz[0],sz[1]) + 1.
  if not keyword_set(GDPIX) then gdpix = where(fltarr(sz[0],sz[1]) EQ 0.)

  nwpix = n_elements(nwwv)

; Calculate wavelength endpoints for each pixel

  ; wvl = Lower endpoint
  ; wvh = Upper endpoint
  if keyword_set(REDBLUE) then begin  ; Wavelengths go red to blue on the CCD
      wvh = (wv + shift(wv, 1, 0))/2.
      wvl = (wv + shift(wv, -1, 0))/2.
      wvh[0,*] = wv[0,*] + (wv[0,*] - wv[1,*])/2.
      wvl[sz[0]-1,*] = wv[sz[0]-1,*] + (wv[sz[0]-2,*] - wv[sz[0]-1,*])/2.
  endif else begin
      wvl = (wv + shift(wv, 1, 0))/2.
      wvh = (wv + shift(wv, -1, 0))/2.
      wvl[0,*] = wv[0,*] - (wv[1,*] - wv[0,*])/2.
      wvh[sz[0]-1,*] = wv[sz[0]-1,*] + (wv[sz[0]-1,*] - wv[sz[0]-2,*])/2.
  endelse

; Refine gdpix to avoid wave=0. values
  a = where(wv[(gdpix-1)>0] * wv[(gdpix+1)<(sz[0]*sz[1]-1)] EQ 0, na)
  if na NE 0 then gdpix[a] = -1
  ;; Also avoid wvh > wvl
  a = where(wvh[gdpix] LE wvl[gdpix], na)
  if na NE 0 then begin
;      print, 'x_rebin2dspec: Warning!  Wavelength solution not monotonic'
;      print, 'x_rebin2dspec: Rejecting bad pixels..'
      gdpix[a] = -1
  endif
  

; Calculate endpoints of the final array
  bwv = dblarr(nwpix+1)
  bwv[0:nwpix-1] = (nwwv + shift(nwwv,1))/2.
  bwv[0] = nwwv[0] - (nwwv[1]-nwwv[0])/2.
  bwv[nwpix] = nwwv[nwpix-1] + (nwwv[nwpix-1]-nwwv[nwpix-2])/2.

; Create tmp arrays for final array

  nwfx  = fltarr(nwpix)
  nwvar = dblarr(nwpix)

;  Loop!

;  if not keyword_set(SILENT) then print, 'x_rebin2dspec: Big Loop!'
  ngd = n_elements(gdpix)
  if not keyword_set( REBINC ) then begin
      ;; Loop on pixels
      for i=0L, ngd-1 do begin
          q = gdpix[i]
                                ;; Rejected Good pixels
          if q EQ -1 then continue
          
          ;; Fill up the Final array
          ;; No overlap
          if wvh[q] LE bwv[0] OR wvl[q] GE bwv[nwpix] then continue
          ;; No zero's
          if (wvh[q] * wvl[q]) EQ 0. then continue
          
          ;; Find pixel that bw is within
          if wvl[q] LT bwv[0] then i1 = [0L] else $
            i1 = where(wvl[q] LE shift(bwv,-1) AND wvl[q] GT bwv)
          ;; Same for hw
          if wvh[q] GT bwv[nwpix] then i2 = [nwpix-1L] else $
            i2 = where(wvh[q] LE shift(bwv,-1) AND wvh[q] GT bwv)
          j1 = i1[0]
          j2 = i2[0]
          
          ;; Now Sum up
          for kk=j1,j2 do begin
              ;; Rejected pixels do not get added in
              if var[q] GT 0. then begin
                  frac = ( (wvh[q] < bwv[kk+1]) - (wvl[q] > bwv[kk]) ) $
                    / (wvh[q]-wvl[q])
                  nwfx[kk] = nwfx[kk] + frac * fx[q]
              endif
              ;; Variance
              if nwvar[kk] NE -1 then begin ; CR junk
                  if var[q] LE 0. then begin
                      if keyword_set( CR ) then begin
                          if var[q] EQ -1. then nwvar[kk] = -1.
                      endif 
                  endif else nwvar[kk] = nwvar[kk] + frac * var[q]
              endif 
          endfor
      endfor
  endif else begin
      ndim = 2
      dim = lonarr(10)
      dim[0] = ngd
      dim[1] = nwpix
      if keyword_set( CR ) then dim[2] = 1L
      
      soname = filepath('libxmath.' + idlutils_so_ext(), $
                        root_dir=getenv('XIDL_DIR'), subdirectory='/lib')
      retval = call_external(soname, 'rebin2dspec', $
                             ndim, dim, gdpix, double(wvl), double(wvh),$
                             bwv, float(fx), double(var), nwfx, nwvar)
  endelse

  return
end
