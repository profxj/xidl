;+ 
; NAME:
; esi_sensstd
;     Version 1.0
;
; PURPOSE:
;    Process a standard star and create the sensitivity file
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

pro esi_echsensstd, std_fil, CLOBBER=clobber, CONDITIONS=conditions, $
                    NOOSCAN=nooscan, NOPROC=noproc, NO_UPD_WEB=no_upd_web, _EXTRA=extra

;
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'esi_echsensstd, std_fil [v1.0]'
;      return
  endif 

  if not keyword_set(WVTIME) then wvtime = [4500., 6200., 7700, 9000]

  ;; Reduce + Extract standard

  if x_chkfil(STD_FIL+'*') NE 1 then begin
      print, 'No ESI file '+std_fil
      return
  endif
  ;; Parse header to create output file name
  head = xheadfits(std_fil)
  date = sxpar(head, 'DATE-OBS')
  date = strmid(date,8,2)+$
         x_getmonth(long(strmid(date,5,2)))+$
         strmid(date,0,4)

  frame = sxpar(head,'FRAMENO')
  outfil = getenv('ESI_CALIBS') + '/STD/sens_ECH_'+date+'_'+ $
           x_padstr(frame,4,'0', /TRIM,/REVERSE)+'.fits'
  
  ;; Check
  if x_chkfil(outfil+'*') and not keyword_set(CLOBBER) then begin
      print, 'esi_echthruput: File '+outfil+' exists.  Use /CLOBB to overwrite.'
      if not keyword_set(NO_UPD_WEB) then begin
          files = findfile(getenv('ESI_CALIBS')+'/STD/sens_ECH*')
          mkhtml_specthru, files, TITLE='ESI ECH Throughput Measurements', $
                           WVTIME=WVTIME, $
                           OUTPTH = getenv('ESI_CALIBS')+'/STD/'
      endif
      return
  endif 

  ;; Process
  if not keyword_set(NOPROC) then esi_echqckrdx, std_fil, /STD, NFIND=1L, NOOSCAN=nooscan;, /NOSKY, /CHK

  ;; Meta data
  esi_strct, esi, IMG=[std_fil], /NOEDIT, /NOFILE
  
  ;; Precess RA/DEC
  x_radec, esi[0].ra, esi[0].dec, rad, decd
  precess, rad, decd, sxpar(head, 'EQUINOX'), 2000.
  x_radec, ras, decs, rad, decd, /flip

  ;;

  meta = {thrustrct}
  meta.observatory = 'Keck'
  meta.instrument = 'ESI'
  meta.grating = 'ECH'
  meta.blocking = ' '
  meta.filename = std_fil
  meta.airmass = esi[0].AM
  meta.date = date
  meta.ra = ras                   ;; J2000
  meta.dec = decs                 ;; J2000
  meta.slit_width = esi[0].slit
  meta.spec_bin = esi[0].rbin

  if keyword_set(CONDITIONS) then meta.conditions = conditions $
  else meta.conditions = ' '

  ;; Get 'approved' standard file
  x_stdcalibfil, meta.ra, meta.dec, calibfil, std_name
  if strlen(calibfil) EQ 0 then begin
      print, 'esi_echsenstd: Get a better standard, please'
      return
  endif
  meta.std_name = std_name

  ;; Read in data
  objfil = 'Extract/Obj_'+esi[0].img_root
  objstr = xmrdfits(objfil, 1, /sile)
  objstr.box_fx = objstr.box_fx / esi[0].exp  ;; electrons per second

  ;; Get sensitivity function
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

  rej_wv = [ $
        [4460, 4524.], $
        [4853, 4872.], $
        [5014, 5048.], $
        [6273, 6301.], $
        [6864, 6961.], $
        [7589, 7737.] $
        ]
  sz_rej = size(rej_wv, /dimen)
        

  sv_wv = dblarr(10,5000L)
  sv_sens = fltarr(10, 5000L)
  sv_sens2 = fltarr(10, 5000L)

  for qq=0L,9 do begin

      ;; Get the standard values
      gdw = where(objstr[qq].box_wv GT 0., npix)
      linterp, std_wv, std_fx, objstr[qq].box_wv[gdw], fx
      objstr[qq].flux[gdw] = objstr[qq].box_fx[gdw] / fx


      ;; Mask/ 
      msk = replicate(1B,n_elements(objstr[qq].box_wv))
      for ss=0L,sz_rej[1]-1 do begin
          bd = where( objstr[qq].box_wv GT rej_wv[0,ss] and $
                      objstr[qq].box_wv LT rej_wv[1,ss], nbd)
          if nbd NE 0 then msk[bd] = 0B
      endfor

      ;; Fit
      gd = where(objstr[qq].box_var GT 0. and msk, nwav, complement = bad)
;      x_splot, objstr[qq].box_wv[gd], objstr[qq].flux[gd], /blo

      sig_fit  = sqrt(objstr[qq].box_var[gd])
      fit = x_fitrej(objstr[qq].box_wv[gd], $
                     objstr[qq].flux[gd], $
                     sig = sig_fit[gd], $
                     fitstr=fitstr, rejpt=rejpt)
;      x_splot, objstr[qq].box_wv[gdw], x_calcfit(objstr[qq].box_wv[gdw], $
;                                                fitstr=fitstr), $
;               ytwo=objstr[qq].flux[gdw], /bloc

      ;; Save
      allwv = where(objstr[qq].box_wv GT objstr[qq].box_wv[gd[0]] and $
                    objstr[qq].box_wv LT objstr[qq].box_wv[gd[nwav-1]], nsv)
      sv_wv[qq,0:nsv-1] = objstr[qq].box_wv[allwv]
      dwv = sv_wv[qq,0:nsv-1] - shift(sv_wv[qq,0:nsv-1],1) 
      dwv[0] = dwv[1]
      sv_sens[qq,0:nsv-1] = x_calcfit(objstr[qq].box_wv[allwv], fitstr=fitstr) ;/ $
      sv_sens2[qq,0:nsv-1] = x_calcfit(objstr[qq].box_wv[allwv], fitstr=fitstr) / $
        dwv  ;;   electrons/Ang/flux (Fnu)
  endfor

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
  allsens = allsens[idx]  > 10.
  allsens2 = allsens2[idx] > 10.
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


  ;; Update the HTML files
  if not keyword_set(NO_UPD_WEB) then begin
      files = findfile(getenv('ESI_CALIBS')+'/STD/sens_ECH*')
      mkhtml_specthru, files, TITLE='ESI ECH Throughput Measurements', $
                       WVTIME=WVTIME, $
                       OUTPTH = getenv('ESI_CALIBS')+'/STD/'
  endif
  
end
              
