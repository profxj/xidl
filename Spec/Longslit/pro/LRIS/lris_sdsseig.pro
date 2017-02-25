;+ 
; NAME:
; lris_sdsseig
;    Version 1.0
;
; PURPOSE:
;   Converts SDSS eigenfunctions into LRIS. Defaults to parameters
;   best for D560 600/10000 (blue) and 600/4000 (red)
;
; CALLING SEQUENCE:
;   
;   lris_sdsseig, sdss_fil, lris_fil[],
;
; INPUTS:
;   sdss_fil -- SDSS (from $IDLSPEC2d/templates/ typically) to
;               rebin to desired LRIS blue and red spectra
;
; RETURNS:
;
; OUTPUTS:
;   lris_fil -- Blue and red spectra from rebinned SDSS tamplates
;
; OPTIONAL KEYWORDS:
;   CRVAL1 -- alog10(min wavelength) for both blue and red spectra
;   CDELT1 -- alog10(1+ <dwv/wv>) for both blue and red spectra
;   NPIX -- number of pixels in new blue and red spectra
;   /VIEW -- plot each rebinned eigenspectrum over original
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   lris_sdsseig, lris, maskid, expsr
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   19-Apr-2002 Written by JXP
;    5-Jun-2008 Adopted from wfccd_sdsseig, KLC
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro lris_sdsseig, sdss_fil, lris_fil, CRVAL1=crval1, CDELT1=cdelt1, $
                  NPIX=npix, VIEW=view, single=single
;
if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
      'lris_sdsseig, sdss_fil, lris_fil[], CRVAL1=[], CDELT1=[], '
    print,'        NPIX=[] [v1.0]'
    return
endif 

;  Optional Keywords

if not keyword_set( CRVAL1 ) then crval1 = [alog10(3200.d),alog10(5747.d)]
if not keyword_set( CDELT1 ) then cdelt1 = [0.000125556,7.8156869e-05]
if not keyword_set( NPIX ) then npix = [2050,2060]

if not keyword_set(single) and $
   n_elements(lris_fil) ne 2 and n_elements(crval1) ne 2 and $
   n_elements(cdelt1) ne 2 and n_elements(npix) ne 2 then $
  stop,'lris_sdsseig: need blue and red values'

;  Input

sdss = mrdfits(sdss_fil, 0, hsdss)
sz_sdss = size(sdss, /dimensions)
wav0 = sxpar(hsdss, 'COEFF0')
wav1 = sxpar(hsdss, 'COEFF1')

;; SDSS Wavelengths
wv_sdss = 10^(wav0 + findgen(sz_sdss[0])*wav1)

;; Loop blue and red
for ii=0,n_elements(crval1)-1 do begin
    ;; New wavelength solution
    wv_lris = 10^(crval1[ii] + findgen(npix[ii])*cdelt1[ii])

    ;; Create final array
    lris = fltarr(npix[ii], sz_sdss[1])
    
    ;; Rebin
    for q=0L,sz_sdss[1]-1 do begin
        x_specrebin, wv_sdss, sdss[*,q], wv_lris, nwfx
        lris[*,q] = nwfx
        if keyword_set(view) then $
          x_splot,wv_sdss,sdss[*,q],xtwo=wv_lris,ytwo=lris[*,q],$
          title=strtrim(lris_fil[ii],2)+' '+strtrim(q,2),/block
    endfor

    ;; Modify header
    head = hsdss
;    mkhdr, head, 4, [npix[ii],sz_sdss[1]]
    sxaddpar, head, 'NAXIS1', npix[ii]
    sxaddpar, head, 'COEFF0', crval1[ii]
    sxaddpar, head, 'COEFF1', cdelt1[ii] ; I hope

    ;; Write
    mwrfits, lris, lris_fil[ii], head, /create
    print,'lris_sdsseig: created ',lris_fil[ii]
endfor                          ;loop n_elements(crval1)

return
end
