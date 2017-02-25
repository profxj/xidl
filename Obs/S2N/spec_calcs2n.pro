;+ 
; NAME:
; spec_calcs2n
;    Version 1.1
;
; PURPOSE:
;     This program computes count rates and expected S/N for
;     a generic slit spectrograph
;
; CALLING SEQUENCE:
;  spec_calcs2n, iwv, flg
;
; INPUTS:
; iwv -- Input wavelengths
; flg -- Flag for the instrument (1=Old HIRES, 2=New HIRES, 3=ESI)
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
;     23-Mar-2011 Generatlized for "all" Keck spectrometers
;     23-Aug-2011 Generatlized for "all" spectrometers
;-
;------------------------------------------------------------------------------

function generate_template, str_obs, wave

  temp = xmrdfits(str_obs.template,1)
                                ; this is how single_spec2mag wants
                                ; the flux units, 1e-17 ergs/cm/cm/s/A
  temp.wavelength *= (1+str_obs.redshift) ; note, we are doing NO bounds checking.
  temp_mstar = single_spec2mag(temp.wavelength,temp.flux*1e17,str_obs.filter)
                                ; this is an AB magnitude in the filter
                                ; we need the AB offset

  if str_obs.mtype EQ 1 then begin
     vega = xmrdfits(str_obs.vega_template,1) 
     vega_mstar = single_spec2mag(vega.wavelength,vega.flux*1e17,str_obs.filter)
     temp_mstar = temp_mstar - vega_mstar
  endif 

  temp.flux *= 10^(0.4*temp_mstar)                         ; force to have 0 mag in filter
  temp.flux *= temp.wavelength /(6.626d-27 *2.99792e18)    ; convert to photons
  n0 = interpol(temp.flux,temp.wavelength,wave,/SPLINE)

  ; minimalist bounds checking
  
  max_index = where(wave GT max(temp.wavelength),nmax)
  if nmax GT 0 then begin
     n0[max_index] = 0
     print,"flagging nmax pixels",nmax
  endif
  min_index = where(wave LT min(temp.wavelength),nmin)
  if nmin GT 0 then begin
     n0[min_index] = 0
     print,"flagging nmin pixels",nmin
  endif

  return,n0
end


function spec_calcs2n, wave, thru, str_tel, str_instr, str_obs

  if  N_params() LT 5 then begin 
      print,'Syntax - ' + $
        'spec_calcs2n, wave, str_tel, str_instr, str_obs, [v1.1]'
      return, -1
  endif 

  ;; More generic Instrument Stuff
  height = str_instr.sheight
  width = str_instr.swidth
  dark = str_instr.dark
  
  binc = str_instr.bind
  binr = str_instr.bins
  
  seeing = str_obs.seeing
  phase = str_obs.mphase
  air = str_obs.airmass
  time = str_obs.exptime

  ;; SED
  if str_obs.mstar LT 1 then mstar = 17. else mstar = str_obs.mstar
  if str_obs.mtype EQ 0 then mtype = 2 else mtype = str_obs.mtype
  
  
  case mtype of 
     1: begin
        nj0     = x_fluxjohnson(wave) ;; photons/s/cm^2/Ang
        n0 = nj0
     end
     2: begin
;          nab0    = 3.54e-9 * wave *1e-8 / c.c / c.h  ;; photons/s/cm^2/Ang
        nab0    = 10.d^(-0.4*48.6) / 6.626d-27 / wave ;; photons/s/cm^2/Ang
        n0 = nab0
     end
     else: stop
  endcase
      
  ;; Spreading out the light
  slit0   = x_gaussslit(width/seeing, height/seeing, 0., 0.)
  rows    = long(3*seeing/str_instr.SCALE_PERP+0.999)
  columns = 2. > ( (width < 3.*seeing)/str_instr.SCALE_PARA )
  nsky    = float((long(height/str_instr.SCALE_PERP+0.999) - rows))/rows     
  
  ;; Readno
  read = str_instr.readno * sqrt(rows / float(binr))
  
  if not keyword_set(NOPRINT) then $
     print, 'Slit width projects to ', columns, ' pixels'
  
 ;;;;;
  ;; S/N
  pixel    = binc*(wave/str_instr.R)  ;; Ang
  slarea   = 3*seeing*width 
  
  ;; Grab extinction
  case strtrim(str_tel.name,2) of 
     'KeckI': extinct  = maunakea_trans(wave)
     'KeckII': extinct  = maunakea_trans(wave)
     'Lick-3m': extinct  = mtham_trans(wave)
     'APF': extinct  = mtham_trans(wave)
     else: stop
  endcase
  slit1    = slit0*10^(-0.4*extinct*air)

  ;; Sky
  case strtrim(str_tel.name,2) of 
     'KeckI': magsky   = maunakea_sky(wave, phase, /NOEMPIR)
     'KeckII': begin
        case str_instr.name of
           'DEIMOS': begin
              phase = 0L ;; Only New Moon so far
              if str_instr.grating EQ '1200' then flg_sky = 1L 
              magsky = maunakea_sky(wave, phase, flg_sky=flg_sky)
           end
           'ESI': begin
              phase = 0L ;; Only New Moon so far
              magsky = maunakea_sky(wave, phase)
           end
           'LRIS': begin
              phase = 0L ;; Only New Moon so far
              magsky = maunakea_sky(wave, phase)
           end
           else: magsky   = maunakea_sky(wave, phase, /NOEMPIR)
        endcase
     end
     'Lick-3m': magsky   = mtham_sky(wave, phase)
     'APF': magsky   = mtham_sky(wave, phase)
     else: stop
  endcase

  if str_obs.template and str_obs.filter then  n0 = generate_template(str_obs,wave)

  ;; Object
  star     = n0*(10^(-0.4*mstar))*str_tel.area*thru*slit1*pixel*time
  dstar    = sqrt(star)
  noise    = read
;  noise    = read/sqrt(float(binr*binc))
  projslit = columns*(wave/str_instr.R)
;  stop
  
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
  ;; Final structure
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

  return, fstrct
end
