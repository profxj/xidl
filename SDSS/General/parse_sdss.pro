;+ 
; NAME:
; parse_sdss   
;   Version 1.1
;
; PURPOSE:
;    Given an SDSS spectrum, return flux, wavelength, etc.
;
; CALLING SEQUENCE:
; parse_sdss, fil, flux, wave, conti, SIG=, IVAR=, ZQSO=, $
;                  NPIX=, HEAD=, CFIL=, CDIR=
;
; INPUTS:
;  fil -- Name of SDSS file OR [PLATE, FIBER]
;
; RETURNS:
;  flux -- Spectral flux array
;  wave -- Wavelength array
;  conti -- Continuum array  (requires CFIL)
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  CFIL= -- Name of continuum file
;  CDIR= -- Name of directory to continuum files
;
; OPTIONAL OUTPUTS:
;  SIG=  -- Sigma array
;  IVAR= -- Inverse variance array
;  ZQSO= -- Redshift of the quasar (or galaxy)
;  NPIX= -- Number of pixels in the spectrum
;  HEAD= -- Full header of the SDSS spectrum
;
; COMMENTS:
;
; EXAMPLES:
;   parse_sdss, fil
;
;
; PROCEDURES/FUNCTIONS CALLED:
; REVISION HISTORY:
;   05-Sep-2003 Written by JXP
;   24-Feb-2014 Updated by KLC ot handle new SDSS spectum format
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro parse_sdss, ifil, flux, wave, conti, SIG=sig, IVAR=ivar, ZQSO=zqso, $
                   NPIX=npix, HEAD=head, CFIL=cfil, CDIR=cdir, DATFIL=fil

;
  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
          'parse_sdss, fil ([PLATE,FIBER]), flux, wave, ' + $
          'conti, SIG=, IVAR=, ZQSO=, NPIX= '
    print, '   HEAD=,  CFIL=, CDIR=  [v1.1]'
    return
  endif 

; Optional Keywords

  if size(ifil,/type) NE 7 then begin
      if n_elements(ifil) NE 2 then begin
          print, 'parse_sdss: Wrong format!'
          return
      endif
      sdss_objinf, ifil, filnm=fil
  endif else fil = ifil

;  Open the file
  dat = xmrdfits(strtrim(fil,2), 0, head, /silent)
  if size(dat,/n_dim) eq 0 then begin
     dat = xmrdfits(strtrim(fil,2), 1, head, /silent)
     if size(dat,/type) ne 8 then $
        stop,'parse_sdss stop: input spectrum of unknown format'

     ;; Assume new SDSS format
     flux = dat.flux
     npix = n_elements(flux)
     wave = 10.^dat.loglam
     ivar = dat.ivar
     sig = flux * 0
     gds = where(ivar gt 0.)
     sig[gds] = sqrt(1./ivar[gds])
  endif else begin
     flux = dat[*,0]
     npix = n_elements(flux)
     
     ;; Wave
     iwave = sxpar(head, 'CRVAL1')
     cdelt = sxpar(head,'CD1_1')
     
;     wave= 10.^(iwave + findgen(npix)*(0.0001));-10^(iwave)+10^(iwave)
     wave= 10.^(iwave + dindgen(npix)*cdelt)
     if abs(cdelt-0.0001d) gt 1.e-7 then $
        print,'parse_sdss: pixel scale not the expected size.',cdelt
     
     ;; SIG
     sig = dat[*,2]

     ;; IVAR
     if arg_present( IVAR ) then begin
        gds = where(sig GT 0.)
        ivar = fltarr(npix)
        ivar[gds] = 1. / sig[gds]^2
     endif

  endelse 

  ;; Zqso
  if arg_present( ZQSO ) then zqso = sxpar(head, 'Z')
  
  
  ;; Continuum
  if keyword_set( CDIR ) then begin
     sdss_name = strtrim( sxpar(head, 'NAME'), 2)+'-'+ $
                 strtrim( sxpar(head, 'FIBERID'), 2)
     cfil = CDIR+sdss_name+'.fits'
  endif
  if keyword_set( CFIL ) then conti = xmrdfits(cfil, /silent)
  
  return
end

