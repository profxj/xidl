;+ 
; NAME:
; sdss_qsoflux
;
; PURPOSE:
;    Measure the QSO fluxes for the QAL structure and fill up the
;    appropriate tag.  Uses SDSS photometry
;
; CALLING SEQUENCE:
;   sdss_qsoflux, qalfil
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
;   04-Jan-2005 Written by JXP
;-
;------------------------------------------------------------------------------
pro sdss_qsoflux, qalfil, drsumm, ISTRT=istrt, CON_DIR=con_dir

  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
             'sdss_qsoflux, qalfil, drsumm, [V1.0]'
    return
  endif 

  if not keyword_set(ISTRT) then istrt = 0L
  if not keyword_set(CON_DIR) then con_dir = 'ABSLIN/'

  ;;Read in QAL file
  qalstr=xmrdfits(qalfil,1,/silent)
  if not tag_exist(qalstr, 'QSO_FLUX') then begin
      dum = replicate({qalcharstrct},n_elements(qalstr))
      copy_struct, qalstr, dum
      qalstr = temporary(dum)
  endif
  nqal = n_elements(qalstr)

  ;; Read in photometry
  dr = xmrdfits(drsumm, 1, /silent)

  ;; Read in Filter curves and Spline
  readcol, getenv('XIDL_DIR')+'SDSS/Photo/u.dat', uwv, utran
  uspl = spl_init(uwv, utran, /double)
  readcol, getenv('XIDL_DIR')+'SDSS/Photo/g.dat', gwv, gtran
  gspl = spl_init(gwv, gtran, /double)
  readcol, getenv('XIDL_DIR')+'SDSS/Photo/r.dat', rwv, rtran
  rspl = spl_init(rwv, rtran, /double)
  readcol, getenv('XIDL_DIR')+'SDSS/Photo/i.dat', iwv, itran
  ispl = spl_init(iwv, itran, /double)

  cnst = x_constants()
  ;; Start loopin on QAL struct
  for i=istrt,nqal-1 do begin
      if i MOD 100 EQ 0 then print, 'i = ', i
      if qalstr[i].z_qso LT 2.2 then continue

      ;; Read data
      parse_sdss, getenv('SDSSPATH')+qalstr[i].file_name, flux, wave, $
        conti, CDIR=con_dir, SIG=sig, NPIX=npix
      dwv = wave - shift(wave,1)
      dwv[0] = dwv[1]
      wrest = wave / (1.+qalstr[i].z_qso)

      ;; Mask sky Lines
      smed = median(sig, 40L)
      gd = where(smed NE 0.)
      nrm_s = fltarr(npix)
      nrm_s[gd] = sig[gd] / smed[gd]
      djs_iterstat, nrm_s, sigma=3., mask=mask
      repl = where(mask EQ 0, nrep)
      if nrep NE 0 then flux[repl] = conti[repl]

      ;; Convert flux to fnu
      fnu = flux * wave^2  / cnst.c * 1e-8
      
      ;; Identify filter
      cen = 1325.*(1.+qalstr[i].z_qso)

      ;; Match with DR
      mtch = where(qalstr[i].plate EQ dr.plate AND $
                   qalstr[i].fiberid EQ dr.fiberid, nm)
      if nm NE 1 then stop
      ;; Take product
      if cen LT 4900. then begin ;; g filter
          freg = where(wave GE 3850. AND wave LE 5500.)
          crv = spl_interp(gwv, gtran, gspl, wave[freg], /double) 
          mag = dr[mtch].psf_g
      endif
      if cen GE 4900. AND cen LT 6825. then begin ;; r filter
          freg = where(wave GE 5400. AND wave LE 7000.)
          crv = spl_interp(rwv, rtran, rspl, wave[freg], /double) 
          mag = dr[mtch].psf_r
      endif
      if cen GT 6825. then begin ;; i filter
          freg = where(wave GE 6600. AND wave LE 8500.)
          crv = spl_interp(iwv, itran, ispl, wave[freg], /double) 
          mag = dr[mtch].psf_i
      endif
      ;; Find fudge factor
      totf = total(crv * fnu[freg] * dwv[freg])  ; This is not quite kosher
      denom = total(crv * dwv[freg])
      avfn = totf / denom
      magf = 10^(-1.*(mag + 48.6)/2.5)
      fact = avfn / magf 
      ;; Find median flamb valu
      reg = where(wrest GE 1300. AND wrest LE 1350., nreg)
      if nreg EQ 0 then stop
      medf = median(flux[reg])
      ;; Save
      qalstr[i].qso_flux = medf / fact * 1e17
  endfor
      
  ;; Write
  mwrfits,qalstr,qalfil,/create
  spawn, 'gzip -f '+qalfil
  ;; All done
  print, 'sdss_qsoflux: All done!'

end

