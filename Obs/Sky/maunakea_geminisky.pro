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


function pivot_wave, wave, trans

  norm = int_tabulated(wave, trans * wave)
  pivot = sqrt(norm)                                  
  pivot = pivot/sqrt(int_tabulated(wave,trans/wave))

  return, pivot


end

pro filtermag,nfilt,filt,wave,flux,filtmag,pivotwave

  filtfile =  getenv('XIDL_DIR') + "/Obs/Sky/subaru" 
  if nfilt eq 0 then begin
     filtfile =  getenv('XIDL_DIR') + "/Obs/Sky/Buser_" 
  endif
  filtfile = filtfile + filt + ".dat"
  readcol, filtfile, filtwave,filttrans,/SILENT
  pivotwave += pivot_wave(filtwave,filttrans)

  filtmag = single_spec2mag( wave, flux,filtfile, careful=0.8)

  return


end

pro compute_correction,finalmags,xwave,xmag,skywave,dmags,skycorrect

  skycorrect = spline(xwave,dmags,skywave)
;  order = 3
;  skymodel = poly_fit(xwave,dmags,order)  
;  skycorrect = dblarr(n_elements(skywave))
;  for i = 0,order do begin
;     skycorrect += skymodel[i]*skywave^i
;  endfor

  bigcorrect = where(skycorrect GT 2,nbig)
  smallcorrect = where(skycorrect LT 0.5,nsmall)

  if nbig gt 0 then skycorrect[bigcorrect] = 2
  if nsmall gt 0 then skycorrect[smallcorrect] = 0.5

  return
end


;------------------------------------------------------------------------------
function maunakea_geminisky, wave, iphase, dispersion=dispersion

;     Sky brightness at Mauna Kea. 
;..............................................................................
;
;  NWAVE=NPHASE 
;  6=5


  if n_elements(wave) eq 0 then begin
     print,"needs an input wavelength array in Angstroms"
     return,-1
  endif

  clr = getcolor(/load)

  disp = 0
  if keyword_set(dispersion) then begin
     if dispersion gt 0 then disp = dispersion
  endif

  filename = getenv('XIDL_DIR') + "/Obs/Sky/skybg_50_10.dat"
  readcol, filename, skywave, skybright, /silent

  planck = 6.62606896e-27; woohoo, ergs!
  c =  2.99798e18;  Angstroms
  skywave *= 10; Angstroms!
  skybright = 1e12 * skybright * planck * c / skywave ; convert to ergs
; flux needs to be in 1e17 ergs for spec2mag


  ;; Max phase
  phase = iphase < 13.9
  ;
  
  
  filtlist = ['U','B','V','R','I','z']
  xphase = [0., 3., 7., 10., 14.]     
  xsky = [ [ 22.4, 21.5, 19.9, 18.5, 17.0], $
           [ 23.0, 22.4, 21.6, 20.7, 19.5], $
           [ 21.9, 21.7, 21.4, 20.7, 20.0], $
           [ 21.2, 20.8, 20.6, 20.3, 19.9], $
           [ 19.9, 19.9, 19.7, 19.5, 19.2], $
           [ 18.3, 18.3, 18.1, 17.8, 17.6] ]
; Vega values execpt for z
; z comes from a comb of two places
;
; http://www.naoj.org/Observing/Instruments/SCam/exptime.html
; for data at three days
;
; http://www.eso.org/observing/etc/doc/ut1/fors/helpfors.html#atmosphere
; for variation with moon
; apparently Paranal is darker than Mauna Kea in the z (as is CTIO).
; I am not sure I believe these, but whatever
;
  ;; Convert to AB
  ;; AB  U_AB = U + 0.71;  B_AB = B-0.11; V_AB = V, R_AB = R+0.199,
   ;; I_AB=I+0.454  [Probably optimal for galaxies]
  
  xsky[*,0] += 0.71
  xsky[*,1] -= 0.11
  xsky[*,3] += 0.199
;  xsky[*,4] += 0.454
  xsky[*,4] += 0.0
  xsky[*,5] += 0.0

  finalmags =  dblarr(n_elements(filtlist))
  xwave = dblarr(n_elements(filtlist))
  xmag = dblarr(n_elements(filtlist))
  for nfilt=0,n_elements(filtlist)-1 do begin
     twave = 0.0D
     tmag = 0.0D
     filtermag,nfilt,filtlist[nfilt],skywave,skybright,tmag,twave
     xwave[nfilt] += twave
     xmag[nfilt] += tmag
     linterp,xphase,xsky[*,nfilt],phase,curmag
     finalmags[nfilt] += curmag
;     print,nfilt," ",filtlist[nfilt],xwave[nfilt],xmag[nfilt],finalmags[nfilt]
  endfor
;  plot,skywave,skybright,color=clr.blue

  speccorrect = dblarr(n_elements(skybright))
  dmags = dblarr(n_elements(filtlist))
  dmags += 2
  iter = 0
  while (max(abs(dmags-1)) GT .05 and iter le 9) do begin

     compute_correction,finalmags,xwave,xmag,skywave,dmags,speccorrect
     skybright /= speccorrect
;     oplot,skywave,skybright

     for nfilt=0,n_elements(filtlist)-1 do begin
        tmag = 0.0D
        twave = 0.0D
        filtermag,nfilt,filtlist[nfilt],skywave,skybright,tmag,twave
        xmag[nfilt] = tmag
;        print,nfilt," ",filtlist[nfilt],xmag[nfilt],finalmags[nfilt]
     endfor
     dmags = xmag-finalmags
     dmags = 10^(-0.4*dmags)

;     print ,iter," ",avg(dmags)-1
     iter += 1
     
  endwhile

  msky= dblarr(n_elements(wave))
  linterp,skywave,skybright,wave,msky
        
  return, msky
end


