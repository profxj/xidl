;+ 
; NAME:
; hires_calcs2n
;    Version 1.1
;
; PURPOSE:
;     This program computes count rates and expected S/N for
;     the HIRES spectrograph.
;
; CALLING SEQUENCE:
;  x_fndchrt, targlist, OUTDIR=, IMSIZE=, SURVEY=
;
; INPUTS:
;  targlist  -- ASCII file containing  (QSO,  RA,  DEC)
;  h0 -- H0 in units of 100 km/s/Mpc
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   imsize - Arcmin of image (default is 5')
;
; OPTIONAL OUTPUTS:
;  OUTDIR=  -- Name of output directory
;
; COMMENTS:
;
; EXAMPLES:
;   x_fndchrt, 'targets.list'
;
; PROCEDURES/FUNCTIONS CALLED:
;  showfits
;  querydss
;
; REVISION HISTORY:
;     HISTORY
;     ??/??/?? Written by G. Donald Penrod 
;     11/13/89 Modified by S. Vogt       - made inputs less confusing
;     06/08/92 Modified by M. Keane      - changed from Hamilton to HIRES 
;     02/07/96 Modified by C. Churchill  - structured queries by function
;                                        - set defaults for Decker C1
;                                        - added comments
;     20-Oct-2005 Ported to IDL by JXP
;-
;------------------------------------------------------------------------------
pro hires_calcs2n, iwv, flg, INFIL=infil, NOPRINT=noprint, PLOT=plot, $
                   PSFILE=psfile, STATE=state, S2N=sn, IORDER=iorder, $
                   BLAZE=blaze, PIXEL=pixel, FSTRCT=fstrct

  if  N_params() LT 1 then begin 
      print,'Syntax - ' + $
        'hires_calcs2n, wave, [flg], GUI=, [v1.0]'
      return
  endif 

  ;; Parse infil
  if not keyword_set(CSZ) then csz = 1.5
  if not keyword_set(FLG) then flg = 1


  if not keyword_set(STATE) then begin
      ;; Instrument + Telescope
      x_inithires, str_instr, flg, STR_TEL=str_tel, $
        DECKER=hires_decker, INFIL=infil ;; OLD HIRES
      ;; Observing conditions
      str_obs = x_obsinit( infil, /MAUNAKEA )
  endif else begin
      str_obs = state.str_obs
      str_instr = state.str_instr
      str_tel = state.str_tel
  endelse
 
 nwv = n_elements(iwv)
 wave = iwv[sort(iwv)]

 ;; Spectral stuff
 m      = long(str_instr.MLAMBDA/wave)
 center = str_instr.MLAMBDA/m
 fsr = center/m
 
 low = where(center LT 3800.,nlow,complement=high,ncomplement=nhigh)
 sep = fltarr(nwv)
 iorder = lonarr(nwv)
 if nlow NE 0 then begin
     sep[low]    = 2.*str_instr.DELY*(str_instr.MLAMBDA/m[low] - $
                                 str_instr.MLAMBDA/(m[low]+1)) 
     iorder[low] = 2
 endif
 if nhigh NE 0 then begin
     sep[high]  = str_instr.DELY*(str_instr.MLAMBDA/m[high] $
                                  - str_instr.MLAMBDA/(m[high]+1)) 
     iorder[high] = 1
 endif
 
  height = str_instr.sheight
  width = str_instr.swidth

  dark = str_instr.dark

  binc = str_instr.bind
  binr = str_instr.bins

  seeing = str_obs.seeing
  phase = str_obs.mphase
  air = str_obs.airmass
  time = str_obs.exptime

  ;; Print
  if not keyword_set(NOPRINT) then begin
      print, 'DECKER  = ', hires_decker
      print, 'BINC    = ', binc
      print, 'BINR    = ', binr
      print, 'SEEING  = ', seeing
      print, 'MOON    = ', phase
      print, 'AIRMASS = ', air
      print, 'EXPTIME = ', time
  endif

  if str_obs.mstar LT 1 then mstar = 17. else mstar = str_obs.mstar
  if str_obs.mtype EQ 0 then mtype = 1 else mtype = str_obs.mtype


;  if mtype EQ 1  then n0 = nj0 else n0 = nab0 
  case mtype of 
      1: begin
          nj0     = x_fluxjohnson(wave)
          n0 = nj0
      end
      2: begin
;          nab0    = 3.54e-9 * wave *1e-8 / c.c / c.h  ;; photons/s/cm^2/Ang
          nab0    = 10.d^(-0.4*48.6) / 6.626d-27 / wave  ;; photons/s/cm^2/Ang
          n0 = nab0
      end
      else: stop
  endcase
  if str_obs.template and str_obs.filter then begin
     temp = xmrdfits(str_obs.template,1)
                                ; this is how single_spec2mag wants
                                ; the flux units, 1e-17 ergs/cm/cm/s/A
     temp.wavelength *= (1+str_obs.redshift) ; note, we are doing NO bounds checking.
     temp_mstar = single_spec2mag(temp.wavelength,temp.flux*1e17,str_obs.filter)
     ; this is an AB magnitude in the filter
     ; we need the AB offset
     temp.flux *= 10^(0.4*temp_mstar)  ; force to have 0 mag in filter
     temp.flux *= temp.wavelength /(6.626d-27 *2.99792e18) ; convert to photons
     n0 = interpol(temp.flux,temp.wavelength,wave,/SPLINE)
  endif
      
  slit0   = x_gaussslit(width/seeing, height/seeing, 0., 0.)
  rows    = long(3*seeing/str_instr.SCALE_PERP+0.999)
  columns = 2. > ( (width < 3.*seeing)/str_instr.SCALE_PARA )
  nsky    = float((long(height/str_instr.SCALE_PERP+0.999) - rows))/rows     

  ;; Readno
  read = str_instr.readno * sqrt(rows / float(binr))

  if not keyword_set(NOPRINT) then $
    print, 'Slit width projects to ', columns, ' pixels'
 
 ;  grab through put from within subroutine and communicate
  case flg of
      1: thru = hires_thruput(wave, center, iorder, fsr, BLAZE=blaze)  ; Old HIRES
      2: thru = hires_thruput(wave, center, iorder, fsr, flg=1) ; New HIRES
      3: thru = mthr_thruput(wave) ; MTHR
      else: stop
  endcase

;
  pixel    = binc*(wave/str_instr.R)
  slarea   = 3*seeing*width 

  ;; Grab extinction
  extinct  = maunakea_trans(wave)
  slit1    = slit0*10^(-0.4*extinct*air)
;
  star     = n0*(10^(-0.4*mstar))*str_tel.area*thru*slit1*pixel*time
  dstar    = sqrt(star)
  noise    = read
;  noise    = read/sqrt(float(binr*binc))
  magsky   = maunakea_sky(wave, phase, /NOEMPIR)
  projslit = columns*(wave/str_instr.R)

  ;;  compute the sky counts and noise
  sky   = n0*(10^(-0.4*magsky))*slarea*str_tel.area*thru*pixel*time
  dsky  = sqrt(sky)
;
;  compute the full signal to noise
  ndark  = binc*dark*rows*time/3600.
  ddark  = sqrt(ndark)
  tnoise = sqrt(star+(1.+1./nsky)*(noise*noise+sky+ndark))
  sn     = star/tnoise
;
  fstrct = { $
             sn: sn, $
             star: star, $
             tnoise: tnoise, $
             extinct: extinct, $
             noise: noise, $
             ndark: ndark, $
             slit0: slit0, $
             thru: thru, $
             sky: sky $
           }

  ;; Loop on wave
  if not keyword_set(NOPRINT) then begin
      for qq=0L,nwv-1 do begin
          print,  'Wavelength    =', wave[qq]
          print,  'Object counts =', star[qq],dstar[qq], $
            format=('(11x,a,f8.0,1x,f7.1)')
          print,  'Sky counts    =', sky[qq],dsky[qq], $
            format=('(11x,a,f8.0,3x,f5.1)')
          print,  'Dark counts   =', ndark,ddark, format=('(11x,a,f8.0,3x,f5.1)')
          print,  'Readout       =', noise, format=('(11x,a,f5.1)')
          print,  'Net Object    =', star[qq],tnoise[qq], $
            format=('(11x,a,f8.0,1x,f7.1)')
          print,  'Net S/N: ', sn[qq], ' per ', 1000.*pixel[qq], 'mA pixel ', $
            format=('(1x, a, f6.0, a, f4.0, a)')
          print, sn[qq]*sqrt(columns), ' per ', 1000.*projslit[qq], $
            'mA res element', $
            format=('(10x,f6.0,a,f4.0,a)')
      endfor
  endif


  return
end
