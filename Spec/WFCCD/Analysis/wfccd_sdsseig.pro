;+ 
; NAME:
; wfccd_sdsseig
;    Version 1.0
;
; PURPOSE:
;   Converts SDSS eigenfunctions into WFCCD
;
; CALLING SEQUENCE:
;   
;   wfccd_sdsseig, wfccd, maskid, expsr, XSIZE=, YSIZE=
;
; INPUTS:
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   XSIZE      - Size of gui in screen x-pixels (default = 1000)
;   YSIZE      - Size of gui in screen y-pixels (default = 600)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   wfccd_sdsseig, wfccd, maskid, expsr
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   19-Apr-2002 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro wfccd_sdsseig, sdss_fil, wfccd_fil, CRVAL1=crval1, CDELT=cdelt, NPIX=npix

;
  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
      'wfccd_sdsseig, sdss_fil, wfccd_fil, CRVAL1=, CDELT=, NPIX= [v1.0]'
    return
  endif 

;  Optional Keywords

  if not keyword_set( WV0 ) then crval1 = alog10(3200.d)
  if not keyword_set( WVDELT ) then cdelt = 1.4486623d-4
  if not keyword_set( NPIX ) then npix = 3200L

;  Input

  sdss = mrdfits(sdss_fil, 0, hsdss)
  sz_sdss = size(sdss, /dimensions)
  wav0 = sxpar(hsdss, 'COEFF0')
  wav1 = sxpar(hsdss, 'COEFF1')

  ; SDSS Wavelengths
  wv_sdss = 10^(wav0 + findgen(sz_sdss[0])*wav1)

; New wavelength solution
  wv_wfccd = 10^(crval1 + findgen(npix)*cdelt)

; Create final array
  wfccd = fltarr(npix, sz_sdss[1])
  
; Rebin

  for q=0L,sz_sdss[1]-1 do begin
      stop
      x_specrebin, wv_sdss, sdss[*,q], wv_wfccd, nwfx
      wfccd[*,q] = nwfx
      ; Plot
;      clr = getcolor(/load)
;      plot, wv_wfccd, nwfx
;      oplot, wv_sdss, sdss[*,q], color=clr.red
  endfor

  mkhdr, head, 4, [npix,sz_sdss[1]]
  sxaddpar, head, 'COEFF0', alog10(3200.d)
  sxaddpar, head, 'COEFF1', 100.d/2.9979e5/alog(10.d)   ; 100km/s pixels

; Write
  mwrfits, wfccd, wfccd_fil, head, /create

  return
end
