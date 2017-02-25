;+ 
; NAME:
; boss_readspec   
;   Version 1.1
;
; PURPOSE:
;    Given an BOSS FITS file (KG format), 
;           return flux, wavelength, etc.
;
; CALLING SEQUENCE:
; boss_readspec, fil, flux, wave, conti, SIG=, IVAR=, ZQSO=, $
;                  NPIX=, HEAD=, CFIL=, CDIR=
;
; INPUTS:
;  fil -- Name of BOSS file OR [PLATE, FIBER]
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
;   boss_readspec, fil
;
;
; PROCEDURES/FUNCTIONS CALLED:
; REVISION HISTORY:
;   05-Sep-2003 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro boss_readspec, ifil, flux, wave, conti, SIG=sig, IVAR=ivar, ZQSO=zqso, $
                   NPIX=npix, HEAD=head, CFIL=cfil, CDIR=cdir, DATFIL=fil

;
  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
          'boss_readspec, fil ([PLATE,FIBER]), flux, wave, ' + $
          'conti, SIG=, IVAR=, ZQSO=, NPIX= '
    print, '   HEAD=,  CFIL=, CDIR=  [v1.1]'
    return
  endif 

; Optional Keywords

  if size(ifil,/type) NE 7 then begin
      if n_elements(ifil) NE 2 then begin
          print, 'boss_readspec: Wrong format!'
          return
      endif
      boss_objinf, ifil, filnm=fil
  endif else fil = ifil

;  Open the file
  dat = xmrdfits(strtrim(fil,2), 1, head, /silent)
  flux = dat.flux
  npix = n_elements(flux)
  sig = fltarr(npix)

; Wave
  wave= 10.^dat.loglam

; SIG
  ivar = dat.ivar
  good = where(ivar GT 0., complement=bad, ncomplem=nbad)
  if nbad GT 0 then sig[bad] = -1.
  sig[good] = 1./sqrt(ivar[good])

; Zqso
  if arg_present( ZQSO ) then stop

; Continuum (Use KG code)
  if keyword_set( CDIR ) then begin
     stop
     boss_name = strtrim( sxpar(head, 'NAME'), 2)+'-'+ $
                 strtrim( sxpar(head, 'FIBERID'), 2)
      cfil = CDIR+boss_name+'.fits'
  endif
  if keyword_set( CFIL ) then conti = xmrdfits(cfil, /silent)

  return
end

