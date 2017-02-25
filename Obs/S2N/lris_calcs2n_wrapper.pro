pro lris_calcs2n_wrapper,wave,s2n,wvmn=wvmn,wvmx=wvmx,slitwidth=slitwidth,dichroic=dichroic, $ 
                         binning=binning,grism=grism,grating=grating,mtype=mtype,$ 
                         mag=mag,airmass=airmass,seeing=seeing,exptime=exptime,$
                         ffilter=filter,template=template,redshift=redshift

  x_initlris, str_instr, INFIL=infil, STR_TEL=str_tel, SLIT=slitwidth, GRISM=grism, GRATING=grating
  str_obs = x_obsinit(infil)
  str_obs.seeing = 1.5  ; 
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
          wave: fltarr(10000L), $
          s2n: fltarr(10000L), $
          wvmn: 3500., $
          wvmx: 8000., $
          flg_wave: 0, $ ; 1=Blue,  2=Red
          infil: '', $
          bdispersers: ['B300','B400','B600','B1200'], $
          rdispersers: ['600/7500','600/10000','600/5000','400/8500','300/5000','150/7500','831/8200','900/5500','1200/7500','1200/9000'], $
          dichroics: ['d46','d56','d50','d68'], $
          slits: ['1.0', '1.5', '2.0', '1.0'], $
          binning: ['1x1', '2x1', '3x1', '2x2'], $
          deckidx: 0L, $
          instr: ['Lris'], $
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

  state.str_instr[0].dichroic = 'd46'
  if (keyword_set(dichroic)) then  state.str_instr[0].dichroic = dichroic

;  state.str_instr[0].binning = '1x1'
;  state.str_instr[1].binning = '1x1'
;  if (keyword_set(binning)) then  state.str_instr[0].binning = binning
;  if (keyword_set(binning)) then  state.str_instr[1].binning = binning

  if (keyword_set(binning)) then begin
     for in_ind = 0,1 do begin
        state.str_instr[in_ind].bins = long(strmid(binning,0,1))
        state.str_instr[in_ind].bind = long(strmid(binning,2,1))
     endfor
  endif

  state.str_instr[0].grating = 'G2'
  if (keyword_set(grism)) then  state.str_instr[0].grating = grism
  state.str_instr[1].grating = '600/7500'
  if (keyword_set(grating)) then  state.str_instr[1].grating = grating

  lris_calcs2n, state.wave[0:state.nwv-1], 3, $
                STATE=state, /nopr, S2N=s2n,  B_FSTRCT=b_fstrct, $
                R_FSTRCT=r_fstrct, $
                BLUE_IDX=blue_idx, RED_IDX=red_idx
  state.s2n[0:state.nwv-1] = s2n
;  STOP
  wave = state.wave[0:state.nwv-1]
  good = where(finite(s2n))
  ssize = size(s2n,/n_elements)
  obj = make_array(ssize)
  sky = make_array(ssize)
  noise = make_array(ssize)
  bsize = size(b_fstrct.sky,/n_elements)

  obj[0:bsize-1] += b_fstrct.star
  obj[bsize:ssize-1] += r_fstrct.star
  sky[0:bsize-1] += b_fstrct.sky
  sky[bsize:ssize-1] += r_fstrct.sky
  noise[0:bsize-1] += b_fstrct.noise
  noise[bsize:ssize-1] += r_fstrct.noise
  
  print,"wave:",wave[good]
  print,"s2n:",s2n[good]
  print,"obj:",obj[good]
  print,"noise:",noise[good]
  print,"sky:",sky[good]
  return
end
