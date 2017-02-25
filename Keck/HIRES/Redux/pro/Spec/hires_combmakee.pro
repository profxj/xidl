;+ 
; NAME:
; hires_combmakee
;    Version 1.1
;
; PURPOSE:
;   Combines multiple exposures of MAKEE output [old CCD only] and shoves the 2D
;   spectrum into a format for continuum fitting.
;
;  The main call is to long_combspec
;
; CALLING SEQUENCE:
;   hires_combspec, hires, setup, obj_id, chip, exp_id
;
; INPUTS:
;   files   -- MAKEE Flux- files assumed
;
; RETURNS:
;
; OUTPUTS:
;   echfspec   -  HIRES fspec structure (fits file; FSpec/name.fits)
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS: ;
; EXAMPLES:
;   hires_combmakee, ['Flux-079.fits', 'Flux-080.fits']
;
; PROCEDURES/FUNCTIONS CALLED:
;  long_combspec
;
; REVISION HISTORY:
;   Sep-2011 Written by JXP
;-
;------------------------------------------------------------------------------
;------------------------------------------------------------------------------
;------------------------------------------------------------------------------

pro hires_combmakee, files, outfil 

;
  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
      'hires_combmakee, files, outfil [v1.0]'
    return
  endif 

;  Optional Keywords
  if not keyword_set(KEYINDX) then keyindx = 0L
  if not keyword_set(MEDINDX) then medindx = 100
  if size(ORDNM,/type) EQ 0 then ordnm = 1
  if not keyword_set(SIGREJ) then sigrej=4.
  if not keyword_set(REJSIG) then rejsig=4.
  if not keyword_set(MINPIX1) then minpix1 = 100
  if not keyword_set(OBJ_NM) then obj_nm = 'a'
  if not keyword_set( SNRMIN ) then snrmin = 2.
  if keyword_set(NOREJ) then REJSIG = -1.


  ;; Read in the Spectra and pack for long_combspec
  nfil = n_elements(files)
  f0 = xmrdfits(files[0])
  sz0 = size(f0, /dimen)
  npix = sz0[0]
  norder = sz0[1]

  a_influx = fltarr(npix, norder, nfil)
  a_invar = fltarr(npix, norder, nfil)
  a_loglam = dblarr(npix, norder, nfil)

  for qq=0L,nfil-1 do begin
     ;; Read flux
     flux = xmrdfits(files[qq],0,head)
     np = n_elements(flux[*,0])
     if np NE npix then stop  ;; If you have to pad, don't pass zero value wavelengths!
     a_influx[0:np-1, *, qq] = flux

     ;; Error array
     ipos = strpos(files[qq], 'Flux')
     if ipos LT 0 then stop
     sfile = strmid(files[qq],0,ipos)+'Var'+strmid(files[qq],ipos+4)
     var = xmrdfits(sfile,0)
     a_invar[0:np-1, *, qq] = var > 0.

     ;; Wavelengths
     hires_rdmakeewave, head, wva
     a_loglam[0:np-1, *, qq] = alog10(wva)
  endfor

  ;; Generate final wavelength array
  ;;  2.0 km/s pixels
  gdw = where(a_loglam GT 0.)
  mnwv = min(a_loglam[gdw], max=mxwv)

  CD1_1 =  2.897300918105d-06
  crval1 = mnwv
  final_npix = (mxwv - mnwv) / cd1_1 + 2
  final_wv = 10.d^(crval1 + (cd1_1 * dindgen(final_npix) ))


  ;; Final Output
  echfspec = { hiresfspecstrct }
  echfspec.nexp = nfil
  echfspec.phys_ordr = lindgen(norder) + 1  ;; Not really physical

  for qq=0L,norder-1 do begin

     print, 'qq: ', qq
     ;; Pack
     influx = reform(a_influx[*,qq,*], npix, nfil)

     inivar = reform(a_invar[*,qq,*], npix, nfil)
     gdv = where(inivar GT 0.)
     inivar[gdv] = 1./inivar[gdv]

     inloglam = reform(a_loglam[*,qq,*], npix, nfil)

     ;; Final wave
     a = where(final_wv LE 10.d^min(inloglam[0,*]),na)
     b = where(final_wv GE 10.d^max(inloglam),nb)
     f_npix = b[0]-na+2
     newloglam = alog10( final_wv[a[na-1]+lindgen(f_npix)] )  ; 2 is for padding

     long_combspec, influx, inivar, inloglam, insky = insky, innivar = innivar $
                    , newloglam = newloglam, newflux = newflux $
                    , newivar = newivar, newnivar = newnivar $
                    , newmask = newmask, newsky = newsky $
                    , iref = iref, SIGREJ = SIGREJ, CHECK = CHECK $
                    , NOSHIFT = NOSHIFT, NOSHARP = NOSHARP, NOREJ = NOREJ $
                    ,MEDSCALE=MEDSCALE , DV=500. $
                    , SN2 = SN2, YMULT = YMULT, DEBUG=DEBUG 
;     x_splot, 10.d^newloglam, newflux, /bloc, xtwo=10.d^inloglam[*,0], ytwo=influx[*,0]
;     x_splot, 10.d^newloglam, newivar, /bloc, xtwo=10.d^inloglam[*,0], ytwo=inivar[*,0]

     ;; Save in structure
     echfspec.npix[qq] = f_npix
     echfspec.wave[0:f_npix-1,qq] = 10.d^newloglam
     echfspec.fx[0:f_npix-1,qq] = newflux
     gd = where(newivar GT 0.)
     echfspec.var[gd, qq]  = 1./newivar[gd]
     print, 'Next order now..'
  endfor

  ;;;; OUTPUT  ;;;;
  hires_wrfspec, echfspec, outfil 

;  close, 56

  print, 'hires_combspec:  Output is in ', outfil
  print, 'hires_combspec:  All done!'


  return
end
  

