;+ 
; NAME:
; deimos_sensstd
;     Version 1.0
;
; PURPOSE:
;    Process a standard star and create the sensitivity file
;
; CALLING SEQUENCE:
;   
;
; INPUTS:
;   esi     -  ESI structure
;
; RETURNS:
;
; OUTPUTS:
;  Image with rdxstdered light removed
;
; OPTIONAL KEYWORDS:
;   DFLAT      - Use Dome flats where possible
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   esi_echrdxstd, esi
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   12-Jan-2008 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro deimos_sensstd, std_fil, CLOBBER=clobber, CONDITIONS=conditions, $
                    NOPROC=noproc, NO_UPD_WEB=no_upd_web, _EXTRA=extra, EXPTIME=exptime

;
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'deimos_sensstd, std_fil [v1.0]'
      return
  endif 

  if not keyword_set(WVTIME) then wvtime = [5000., 6200., 8000, 9000]
  if not keyword_set(EXPTIME) then exptime = 1. 


  if x_chkfil(STD_FIL+'*') NE 1 then begin
      print, 'No DEIMOS file '+std_fil
      print, 'Returning...'
      return
  endif

  ;; Grab meta data
  t_meta = xmrdfits(std_fil,1,/silen)
  meta = {thrustrct}
  copy_struct, t_meta, meta
  meta.grating = strtrim(meta.grating,2)
  mgrat = meta.grating

  outfil = getenv('DEIMOS_CALIBS') + '/'+mgrat+'/sens_'+mgrat+ $
           '_'+meta.date+'_'+ x_padstr(meta.frame,4,'0', /TRIM,/REVERSE)+'.fits'

  ;; Create the directory if necessary
  a = findfile(getenv('DEIMOS_CALIBS')+'/'+mgrat+'/..', count=count)
  if count EQ 0 then file_mkdir, getenv('DEIMOS_CALIBS')+'/'+mgrat
  if count EQ 0 then file_mkdir, getenv('DEIMOS_CALIBS')+'/'+mgrat+'/THRU_DIR'

  
  ;; Check
  if x_chkfil(outfil+'*') and not keyword_set(CLOBBER) then begin
      print, 'esi_echthruput: File '+outfil+' exists.  Use /CLOBB to overwrite.'
      if not keyword_set(NO_UPD_WEB) then begin
          files = findfile(getenv('DEIMOS_CALIBS')+'/'+mgrat+'/sens_*')
          mkhtml_specthru, files, TITLE='DEIMOS '+mgrat+' Throughput Measurements', $
                           WVTIME=WVTIME, $
                           OUTPTH = getenv('DEIMOS_CALIBS')+'/'+mgrat+'/'
      endif
      return
  endif 

  if keyword_set(CONDITIONS) then meta.conditions = conditions $
  else meta.conditions = ' '

  ;; Get 'approved' standard file
  x_stdcalibfil, meta.ra, meta.dec, calibfil, std_name
  if strlen(calibfil) EQ 0 then begin
      print, 'deimos_senstd: Get a better standard, please'
      return
  endif
  meta.std_name = std_name

  ;; Read in data
  spectrum = xmrdfits(std_fil, 2, /silen)
  box_wave = spectrum.wavelength
  box_fx = spectrum.counts / exptime ;; electrons/s/pix

  ;; Get sensitivity function
  cal_fil = getenv('XIDL_DIR')+$
            '/Spec/Longslit/calib/standards/calspec/'+calibfil
  std = xmrdfits(cal_fil,1)
  std_wv = std.wavelength
  c = x_constants()
  std_fx = std.flux * (std_wv)^2 * 1d-8 / c.c ;; fnu
  std_ab = -2.5*alog10(std_fx) - 48.6  ;; AB mag

  ;; Fit smooth function to the sensitivity function

  fitstr = x_setfitstrct(FUNC='LEGEND', NORD=4L, LSIG=3., HSIG=3., $
                         FLGREJ=1)

  rej_wv = [ $
;        [4460, 4524.], $
;        [4853, 4872.], $
;        [5014, 5048.], $
;        [6273, 6301.], $
        [6884, 6991.], $
        [7584, 7767.] $
        ]
  sz_rej = size(rej_wv, /dimen)
        
  npix = n_elements(box_wave)
  sv_wv = dblarr(npix)
  sv_sens = fltarr(npix)
  sv_sens2 = fltarr(npix)

  ;; Get the standard values
  linterp, std_wv, std_fx, box_wave, fx
  box_fx = box_fx / fx
  stop

  ;; Mask 
  msk = replicate(1B,npix)
  for ss=0L,sz_rej[1]-1 do begin
      bd = where( box_wave GT rej_wv[0,ss] and $
                  box_wave LT rej_wv[1,ss], nbd)
     if nbd NE 0 then msk[bd] = 0B
  endfor

;  gd = where(msk, nwav, complement = bad)
  ;; Fit sensitivity function
;  fit = x_fitrej(box_wave, box_fx, fitstr=fitstr, rejpt=rejpt)
;  stop

  ;; Smooth
  fit = smooth(box_fx, 11)

  ;; Save
  allwv = box_wave
  dwv = allwv - shift(allwv,1) 
  dwv[0] = dwv[1]
  allsens = 1./fit
  allsens2 = 1./(fit / dwv)  ;;   1/[electrons/Ang/flux] (Fnu)

  ;; Convert to AB
  allsens = -2.5*alog10(allsens) - 48.6
  allsens2 = -2.5*alog10(allsens2) - 48.6

  ;; Chop down to 1000 evaluations
  npix = n_elements(allwv)
  step = npix/1000.
  idx = round(lindgen(999)*step)
  allsens = allsens[idx]
  allsens2 = allsens2[idx]
  allwv = allwv[idx]

  if keyword_set(PLOT) then x_splot, allwv, allsens, /blo
  if keyword_set(PLOT) then x_splot, allwv, allsens2, /blo

  ;; Extinction
  mkhdr, head, allwv
  sxaddpar, head, 'TELESCOP', 'Keck'
  sxaddpar, head, 'ELAPTIME', 1.
  sxaddpar, head, 'AIRMASS', meta.airmass
  extinct = long_extinct(allwv, head, /NOTIME)

  ;; Efficiency
  x_initkeck, telescope 
  linterp, std_wv, std_fx, allwv, all_std  ;; fnu (erg/s/cm^2/Hz)
  nphot = telescope.area * all_std * 1. / (c.h * allwv)         ;; Trust me
  nobs = all_std / 10.^(-1.*(allsens2+48.6)/2.5) 
  eff = (nobs*extinct)/nphot
  if keyword_set(PLOT) then x_splot, allwv, (nobs/nphot), /bloc

  ;; Output
  out_str = { $
            wav: allwv, $
            zp_pix: allsens, $
            zp_ang: allsens2, $
            eff: eff $
            }
  mwrfits, meta, outfil, /create
  mwrfits, out_str, outfil

  ;; Update the HTML files
  if not keyword_set(NO_UPD_WEB) then begin
      files = findfile(getenv('DEIMOS_CALIBS')+'/'+mgrat+'/sens_*')
      mkhtml_specthru, files, TITLE='DEIMOS '+mgrat+' Throughput Measurements', $
                       WVTIME=WVTIME, $
                       OUTPTH = getenv('DEIMOS_CALIBS')+'/'+mgrat+'/'
  endif
  
end
              
