;+ 
; NAME:
; lris_sensstd
;     Version 1.1
;
; PURPOSE:
;    Process a standard star and create the sensitivity file for LRIS
;
; CALLING SEQUENCE:
;   
;
; INPUTS:
;   std_fil -- A processed standard star using the LowRedux pipeline
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

pro lris_sensstd, std_fil, CLOBBER=clobber, CONDITIONS=conditions, $
                  NOPROC=noproc, ALL_GRATING=all_grating, STD_OBJ=std_obj, $
                  CALIBFIL=calibfil

;
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'lris_sensstd, std_fil [v1.0]'
;      return
  endif 


  ;; Reduce + Extract standard

  if x_chkfil(STD_FIL+'*') NE 1 then begin
      print, 'No Standard file '+std_fil
      return
  endif
  ;; Parse header to create output file name
  head = xheadfits(std_fil)
  if strmatch(sxpar(head,'INSTRUME'),'LRISBLUE') then cside = 'b' else cside='r'
  date = sxpar(head, 'DATE')
;  stop
  date = strmid(date,8,2)+$
         x_getmonth(long(strmid(date,5,2)))+$
         strmid(date,0,4)

  ;; Meta data
  meta = {thrustrct}
  meta.observatory = 'Keck'
  meta.instrument = sxpar(head,'INSTRUME')
  if strmatch(meta.instrument,'LRISBLUE',/fold_case) then $
    meta.grating = sxpar(head,'GRISNAME') else $
    meta.grating = sxpar(head,'GRANAME')

  prs = strsplit(strtrim(meta.grating,2),'/', /extract)
  lbl = strtrim(strjoin([cside,prs[0], '_', prs[1]]),2)

  frame = sxpar(head,'FRAMENO')
  outfil = getenv('LRIS_THRU') + '/'+lbl+ '/sens_LRIS'+lbl+'_'+date+'_'+ $
           x_padstr(frame,4,'0', /TRIM,/REVERSE)+'.fits'

  ;; Create the directory if necessary
  a = findfile(getenv('LRIS_THRU')+'/'+lbl+'/..', count=count)
  if count EQ 0 then file_mkdir, getenv('LRIS_THRU')+'/'+lbl
  if count EQ 0 then file_mkdir, getenv('LRIS_THRU')+'/'+lbl+'/THRU_DIR'
  
  ;; Check
  if x_chkfil(outfil+'*') and not keyword_set(CLOBBER) then begin
      print, 'lris_sensstd: File '+outfil+' exists.  Use /CLOBB to overwrite.'
      return
  endif 

  ;; RA/DEC
  ras = sxpar(head, 'RA')
  decs = sxpar(head, 'DEC')
  x_radec, ras, decs, rad, decd


  ;; Airmass
  hangl = float(strsplit(sxpar(head,'HA'),':',/extrac))
  if hangl[0] LT 0. then hangl = hangl[0]-hangl[1]/60. $
  else hangl = hangl[0]+hangl[1]/60. 
  meta.airmass = airmass(40., decd, hangl)

  meta.date = date
  meta.ra = ras                   ;; J2000
  meta.dec = decs                 ;; J2000
  mask = sxpar(head,'SLITNAME')
  prs = strsplit(mask, '_', /extra)
  np = n_elements(prs)
  if np NE 2 then begin
      if strmatch(strtrim(mask,2), 'direct') then meta.slit_width=9.999
  endif else meta.slit_width = float(prs[1])
  meta.central_wave = sxpar(head,'WAVELEN')
  meta.blocking = strtrim(sxpar(head,'DICHNAME'),2)
  binning = long(strsplit(sxpar(head, 'BINNING'), ',', /extract))
  meta.spec_bin = binning[1]
  meta.filename = ' '

  if keyword_set(CONDITIONS) then meta.conditions = conditions $
  else meta.conditions = ' '

  ;; Get 'approved' standard file
  x_stdcalibfil, meta.ra, meta.dec, icalibfil, std_name, TOLER=1000.
  if not keyword_set(CALIBFIL) then calibfil=icalibfil
  if strlen(calibfil) EQ 0 then begin
      print, 'kast_echsenstd: Get a better standard, please'
      return
  endif
  meta.std_name = std_name

  ;; Read in data
  objstr = xmrdfits(std_fil, 5, /sile)
  mx = max(objstr.peakflux,pimx)
  if not keyword_set(STD_OBJ) then imx = pimx else imx = STD_OBJ
  print, 'lris_sensstd: Using object -- ', imx
  objstr=objstr[imx]
  ;; Reverse if need be
  if objstr.wave_box[0] GT objstr.wave_box[1] then begin
      objstr.wave_box = reverse(objstr.wave_box)
      objstr.flux_box = reverse(objstr.flux_box)
      objstr.ivar_box = reverse(objstr.ivar_box)
  endif
  ;;
  exp = sxpar(head,'EXPOSURE') > sxpar(head,'EXPTIME')
  exp = sxpar(head,'TTIME') > exp
  objstr.flux_box = objstr.flux_box / exp  ;; electrons per second
  npix = n_elements(objstr.flux_box)

  ;; Get sensitivity function
  print, 'lris_sensstd: Using the calibration file: ', calibfil
  cal_fil = getenv('XIDL_DIR')+$
            '/Spec/Longslit/calib/standards/calspec/'+calibfil
  std = xmrdfits(cal_fil,1)
  std_wv = std.wavelength
  c = x_constants()
  std_fx = std.flux * (std_wv)^2 * 1d-8 / c.c ;; fnu
  std_ab = -2.5*alog10(std_fx) - 48.6  ;; AB mag

  ;; 

  fitstr = x_setfitstrct(FUNC='LEGEND', NORD=11L, LSIG=3., HSIG=3., $
                         FLGREJ=1)

;  rej_wv = [ $
;        [4460, 4524.], $
;        [4853, 4872.], $
;        [5014, 5048.], $
;        [6273, 6301.], $
;        [6864, 6961.], $
;        [7589, 7737.] $
;        ]
;  sz_rej = size(rej_wv, /dimen)
;  STOP

  sv_wv = dblarr(npix)
  sv_sens = fltarr(npix)
  sv_sens2 = fltarr(npix)

  ;; Get the standard values
  gdw = where(objstr.wave_box GT 0., npix)
  linterp, std_wv, std_fx, objstr.wave_box[gdw], fx
  objstr.flux_opt = 0.
  objstr.flux_opt[gdw] = objstr.flux_box[gdw] / fx


  ;; Mask/ 
  msk = replicate(1B,n_elements(objstr.wave_box))
;  for ss=0L,sz_rej[1]-1 do begin
;      bd = where( objstr.box_wv GT rej_wv[0,ss] and $
;                  objstr.box_wv LT rej_wv[1,ss], nbd)
;      if nbd NE 0 then msk[bd] = 0B
;  endfor

  ;; Fit
  gd = where(objstr.ivar_box GT 0. and msk, nwav, complement = bad)
;  x_splot, objstr.wave_box[gd], objstr.flux_opt[gd], /blo

  sig_fit  = sqrt(1./objstr.ivar_box[gd])
  fit = x_fitrej(objstr.wave_box[gd], $
                 objstr.flux_opt[gd], $
                 sig = sig_fit[gd], $
                 fitstr=fitstr, rejpt=rejpt)
;      x_splot, objstr.box_wv[gdw], x_calcfit(objstr.box_wv[gdw], $
;                                                fitstr=fitstr), $
;               ytwo=objstr.flux[gdw], /bloc

  ;; Save
  allwv = where(objstr.wave_box GT objstr.wave_box[gd[0]] and $
                objstr.wave_box LT objstr.wave_box[gd[nwav-1]], nsv)
  sv_wv[0:nsv-1] = objstr.wave_box[allwv]
  dwv = sv_wv[0:nsv-1] - shift(sv_wv[0:nsv-1],1) 
  dwv[0] = dwv[1]
  sv_sens[0:nsv-1] = x_calcfit(objstr.wave_box[allwv], fitstr=fitstr) ;/ $
  sv_sens2[0:nsv-1] = x_calcfit(objstr.wave_box[allwv], fitstr=fitstr) / $
                      dwv  ;;   electrons/Ang/flux (Fnu)

  ;; Eliminate overlap
  iwv = where(sv_wv GT 0., nwv) 
  allwv = sv_wv[ iwv ]
  allsens = sv_sens[iwv]
  allsens2 = sv_sens2[iwv]

  srt = sort(allwv)
  allwv = allwv[srt]
  allsens = allsens[srt]
  allsens2 = allsens2[srt]


  msk = replicate(1B,nwv)
  for ii=0L,nwv-2 do begin
      if abs(allwv[ii]-allwv[ii+1]) LT 1e-3 then begin
          allsens[ii:ii+1] = total([allsens[ii],allsens[ii+1]])
          allsens2[ii:ii+1] = total([allsens2[ii],allsens2[ii+1]])
          msk[ii+1] = 0B
      endif
  endfor

  bad = where(allsens < 0)
  MSK[bad] = 0
  idx = where(MSK, nwv) 
  allwv = allwv[idx]
  allsens = 1./allsens[idx]
  allsens2 = 1./allsens2[idx]

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


  
end
              
