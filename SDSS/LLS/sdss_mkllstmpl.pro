;+ 
; NAME:
; sdss_mkllstmpl
;
; PURPOSE:
;    Routine to make a series of LLS absorption templates by only
;    varying NHI
;
; CALLING SEQUENCE:
;
; INPUTS:
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   27-Feb-2004 Written by JXP
;-
;------------------------------------------------------------------------------
pro sdss_mkllstmpl, tautrial, outfil

IF NOT keyword_set(tautrial) THEN begin
    nhv = 16. + findgen(20)*0.2
    tautrial=10^(nhv-20.2)
endif
;IF NOT keyword_set(btrial) THEN btrial=[10,15,20,25,30,40,50,80,150]
IF NOT keyword_set(outfil) THEN outfil = 'nhi16_19b30.fits'
IF NOT keyword_set(npoly) THEN npoly=2
IF NOT keyword_set(taucut) THEN taucut=0.0

;open template spectrum
;template=mrdfits(QSOtemplate,/silent)
;tlogwav = findgen(n_elements(template))*0.0001 + QSOtempwav0

;open optical depth spectrum, default NHI = 20.2, b=30.0 from Burles
;inwav = mrdfits('tau_lyman_series.fits',0)
intau = xmrdfits('tau_lyman_series.fits',1)
npin = n_elements(intau)

inwav = 10^(2.7d + dindgen(npin)*0.0001/50.)
sdss_wav = 10^(2.7d + dindgen(npin/50L)*0.0001)
npix = n_elements(sdss_wav)
dwv = sdss_wav - shift(sdss_wav,1)
dwv[0] = dwv[1]

fwhm = 2.
for qq=0L,n_elements(tautrial)-1 do begin
    ;; Update tau
    newtau = intau * tautrial[qq]
    subfx = exp(-newtau)

    fx = fltarr(npix)
    ;; 
    soname = filepath('libxmath.so', $
                      root_dir=getenv('XIDL_DIR'),  $
                      subdirectory='/lib')
    retval = call_external(soname, 'x_gsmooth', $
                           npix, npin, double(sdss_wav), double(inwav), $
                           double(dwv), float(subfx), float(fx), $
                           float(FWHM))
    ;; Write
    mwrfits, fx, outfil, create=(qq EQ 0)

endfor 

 return

END
