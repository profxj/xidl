pro esi_calcs2n_wrapper,wave,s2n,wvmn=wvmn,wvmx=wvmx,slitwidth=slitwidth, $ 
                        binning=binning,mag=mag,airmass=airmass,seeing=seeing,exptime=exptime, $
                        ffilter=filter,template=template,redshift=redshift,mtype=mtype

  x_initesi, str_instr,  INFIL=infil, STR_TEL=str_tel, SLIT=slitwidth
  str_obs = x_obsinit(infil)
;  if (keyword_set(GRATING)) then grating = grating
  str_obs.seeing = 0.7  ; 
  if (keyword_set(seeing)) then str_obs.seeing = seeing
  if (keyword_set(airmass)) then str_obs.airmass = airmass
  if (keyword_set(exptime)) then str_obs.exptime = exptime
  if (keyword_set(mag)) then str_obs.mstar = mag
  if (keyword_set(redshift)) then str_obs.redshift = redshift
  if (keyword_set(mtype)) then str_obs.mtype = mtype
  if (keyword_set(template)) then begin 
     str_obs.template = getenv('TEMPLATE_DIR') + template
     str_obs.vega_template = getenv('TEMPLATE_DIR') + "alpha_lyr_stis_005.fits"
  endif
  if (keyword_set(filter)) then str_obs.filter = getenv('FILTER_DIR') +filter

  state = {             $
          nwv: 0L, $
          dwv: 10., $
          wave: fltarr(1000L), $
          s2n: fltarr(1000L), $
          wvmn: 4000., $
          wvmx: 10000., $
          infil: '', $
          slits: ['0.3','0.5', '0.75', '1.0'], $
          binning: ['1x1','2x1', '3x1', '2x2'], $
          deckidx: 0L, $
          instr: ['ESI'], $
          pixel: 0., $
          flg_plot: 0, $
          side: 0, $
          str_instr: str_instr, $ 
            str_tel: str_tel, $
            str_obs: str_obs $
  }

  if (keyword_set(wvmn)) then state.wvmn = wvmn
  if (keyword_set(wvmx)) then state.wvmx = wvmx

  state.nwv = long((state.wvmx-state.wvmn)/state.dwv) + 1
  state.wave[0:state.nwv-1] = state.wvmn + findgen(state.nwv)*state.dwv

;  stop
  if (keyword_set(binning)) then begin
;     for in_ind = 0,1 do begin
     state.str_instr.bins = long(strmid(binning,0,1))
     state.str_instr.bind = long(strmid(binning,2,1))
;     endfor
  endif

  keck_calcs2n, state.wave[0:state.nwv-1], 3, $
                STATE=state, /nopr, S2N=s2n, FSTRCT=fstrct

;  state.s2n[0:state.nwv-1] = s2n
;  STOP
  wave = state.wave[0:state.nwv-1]
  good = where(finite(s2n))
  ssize = size(s2n,/n_elements)
  obj = make_array(ssize)
  sky = make_array(ssize)
  noise = make_array(ssize)
  size = size(fstrct.sky,/n_elements)

  obj[0:ssize-1] += fstrct.star
  sky[0:ssize-1] += fstrct.sky
  noise[0:ssize-1] += fstrct.noise
  
  print,"wave:",wave[good]
  print,"s2n:",s2n[good]
  print,"obj:",obj[good]
  print,"noise:",noise[good]
  print,"sky:",sky[good]
  return
end
