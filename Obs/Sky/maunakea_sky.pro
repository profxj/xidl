;+ 
; NAME:
; manuakea_sky
;    Version 1.1
;
; PURPOSE:
;    Gives an estimate of the sky brightness for Mauna Kea
;
; CALLING SEQUENCE:
;  mag = maunakea_sky( wave, phase )
;
; INPUTS:
;  wave=  -- Wavelength array
;  phase= -- Phase of the moon [days]
;
; RETURNS:
;  mag=  -- Sky brightness in AB mags
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
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   27-Oct-2005 Written by JXP based on HIRES S2N code
;-
;------------------------------------------------------------------------------
function maunakea_sky, wave, iphase, NOEMPIRI=noempiri, FLG_SKY=flg_sky

common maunakea_sky_cmmn, mkea_swv, mkea_sfx, mkea_sfil, mkea_wvmnx
;     Sky brightness at Mauna Kea. 
;..............................................................................
;
;  NWAVE=NPHASE 
;  6=5

  nwv = n_elements(wave)
  msky = fltarr(nwv)
  if not keyword_set(FLG_SKY) then flg_sky = 0

  if not keyword_set(NOEMPIRIC) then begin
     ;; Read in Empirical sky model (will need to vary with moon phase)
      dir = getenv('XIDL_DIR')
      if dir eq '' then dir = '/local/home/randyc/idl/xidl/' ; for alamoana
      ;;
      case flg_sky of
         0: tfil = dir+'Obs/Sky/Empirical/mkea_sky_newmoon_DEIMOS_600_2011oct.fits'
         1: tfil = dir+'Obs/Sky/Empirical/mkea_sky_newmoon_DEIMOS_1200_2011oct.fits'
         else: stop
      endcase
      if not keyword_set(MKEA_SFIL) then mkea_sfil = ''
      if not strmatch(mkea_sfil,tfil) then begin
          mkea_swv = xmrdfits(tfil, 0, /silent)
          mkea_sfx = xmrdfits(tfil, 1, /silent)
          mkea_sfil = tfil
          mkea_wvmnx = [min(mkea_swv, max=mx), mx]
      endif

     emp_wave = where(wave GE mkea_wvmnx[0] and wave LE mkea_wvmnx[1], ngd, $
                      complement=extrap, ncomplement =nextrap)
  endif else begin  ;; Don't use the Empirical model
     ngd = 0L
     extrap = lindgen(nwv)
     nextrap = nwv
  endelse

  ;; Empirical?
  if ngd GT 0 then begin
     ;; Interpolate
     flam = interpol(mkea_sfx, mkea_swv, wave[emp_wave])
     ;a = where(wave LT mkea_swv, na)
     ;if na GT 0 then flam[a] = mkea_sfx[0]
     
     fnu = flam / (3e10) * wave[emp_wave] * (wave[emp_wave] * 1e-8)
     msky[emp_wave] = -1. * alog10(fnu) / 0.4 - 48.6
  endif

  if nextrap GT 0 then begin 

     ;; Max phase
     phase = iphase < 13.9
;
     xwave = [3500, 4200, 5500, 6700, 7800, 22000.]
     xphase = [0., 3., 7., 10., 14.]     
     xsky = [ [22.4, 23.0, 21.9, 21.2, 19.9, 12.0], $
              [21.5, 22.4, 21.7, 20.8, 19.9, 12.0], $ 
              [19.9, 21.6, 21.4, 20.6, 19.7, 12.0], $
              [18.5, 20.7, 20.7, 20.3, 19.5, 12.0], $
              [17.0, 19.5, 20.0, 19.9, 19.2, 12.0] ] ;; Vega values
     
;
     ;; Convert to AB
     ;; AB  U_AB = U + 0.71;  B_AB = B-0.11; V_AB = V, R_AB = R+0.199,
     ;; I_AB=I+0.454  [Probably optimal for galaxies]
     xsky[0,*] = xsky[0,*] + 0.71
     xsky[1,*] = xsky[1,*] - 0.11
     xsky[3,*] = xsky[3,*] + 0.199
     xsky[4,*] = xsky[4,*] + 0.454
;  stop  ;; Might be a bug here -- Fixed by JXP on 5/25/2011 (I hope)
     
     for ii=0L,nextrap-1 do begin
        qq = extrap[ii]
        gd = where(xwave LE wave[qq], ngd)
        if ngd LE 1 then j = 0 else begin
           mn = min(wave[qq] - xwave[gd], mnj)
           j = gd[mnj]
        endelse
        
        gd = where(xphase LE phase)
        if ngd LE 1 then k = 0 else begin
           mn = min(phase - xphase[gd], mnk)
           k = gd[mnk]
        endelse 
;
;
        t = (wave[qq] - xwave[j])/(xwave[j+1]-xwave[j])
        u = (phase - xphase[k])/(xphase[k+1]-xphase[k])
;       
        msky[qq] = (1.-t)*(1.-u)*xsky[j,k] + t*(1-u)*xsky[j+1,k] $
                   + t*u*xsky[j+1,k+1] +(1-t)*u*xsky[j,k+1]
     endfor
  endif 
;        
  return, msky
end
