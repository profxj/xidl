;+
; NAME:
; cldy_wrin
;  V1.1
;
; PURPOSE:
;    Create a CLOUDY input file
;   
; INPUT
;   params:  array of parameters in the following order:  z, logN_HI, metallicity, logU, logn_H
;   fluxfil:  text file of rydberg/flux pairs
; 
; OPTIONAL
;   v0602a: new syntax for Cloudy v06.02
;   cie: constant temperature, mmodel collisional ionization equilibrium 
;        (exclude ionization parameter and fluxfil)
;     /NOU:  Do not vary U.  Good for varying n_H only
;
; CALLING SEQUENCE:
;
;   cldy_wrin, params, fluxfil, OUTFIL=outfil
;
; REVISION HISTORY:
;   10-Sep-2005 Written by GEP
;   13-Apr-2006 added v0602a and cie, KLC
;-
;------------------------------------------------------------------------------

pro cldy_wrin, params, fluxfil, OUTFIL=outfil, TITLE=title, $
               OLD_VER=old_ver,cie=cie, NODOUB=nodoub, SHAPECMD=shapecmd, $
               FNU=FNU, FLG_NH=flg_NH, HEAT=heat

  if (N_params() LT 1) then begin
      print, 'Syntax - cldy_wrin, params, fluxfil'
      return
  endif

if not keyword_set(OUTFIL) then outfil = 'cloudy.in'
if not keyword_set(TITLE) then title = 'model'

openw, 1, outfil

printf, 1, 'title "'+title+'"'

printf, 1, 'radius 30'          ;log inner radius (cm)
printf, 1, 'print last'         ;Don't print results from iterations
if not keyword_set(OLD_VER) then printf,1,'CMB redshift '+strtrim(string(params[0]),2) $
else printf, 1, 'fireball [z='+strtrim(string(params[0]),2)+']'

  if not keyword_set(NODOUB) then $
    printf, 1, 'double optical depths'

if not keyword_set(CIE) then begin
    if keyword_set(v0602a) then printf, 1, 'stop temperature 10K linear' $
    else printf, 1, 'stop temperature = 10K [linear]'
endif

printf, 1, 'hden '+strtrim(string(params[4]),2)

;; metallicity scale factor (if metals < 10, assumed to be log scale factor)
printf, 1, 'metals '+strtrim(string(params[2]),2) 
case flg_NH of 
   0: printf, 1, 'stop neutral column density '+strtrim(string(params[1]),2) ;log N(HI) (cm^-2)
   1: printf, 1, 'stop column density '+strtrim(string(params[1]),2) ;log N(H) (cm^-2)
   else: stop
endcase
if keyword_set(v0602a) then printf, 1, 'set trimming -20' $
else printf, 1, 'set trim -20'

flg_shape = 0
if keyword_set(cie) then begin
    ;;Collisional ionization equilibrium (if cie < 10, assumed to be logT)
    printf,1,'constant temperature, t='+strtrim(string(cie),2)+'K'
    if not keyword_set(FLUXFIL) then flg_shape =1
endif

if keyword_set(SHAPECMD) then begin
    printf, 1, shapecmd
    printf, 1, 'ionization parameter = ' + strtrim(string(params[3]),2)
    flg_shape = 1
endif

if keyword_set(HEAT) then $
   printf, 1, 'hextra  ' + strtrim(string(HEAT),2) 

if flg_shape NE 1 then begin
    if not keyword_set(FLUXFIL) then stop
                                ;printf,1,'table HM05 redshift '+strtrim(string(params[0]),2) + $
                                ;       ' quasar'$
    ;;log(U) (since hden (n(h)) set, this changes flux of ionizing photons)
    ;; Normalize
    if not keyword_set(FNU) then $
       printf, 1, 'ionization parameter = ' + strtrim(string(params[3]),2) $
    else printf, 1, 'f(nu) = ' + string(fnu, '(f7.3)')

    ;;File of format:
    ;;interpolate (nu(Ryd)  log(f_nu))
    ;;continuee (nu(Ryd)  log(f_nu))...
    readcol, fluxfil, first, ryd, flux, FORMAT='A,A,A', /silent
    
    nflux = n_elements(flux)
    if nflux EQ 0 then stop, 'No flux values!'

    for i = 0, nflux-1 do begin
        printf, 1, first[i]+' '+ryd[i]+' '+flux[i]
    endfor
endif

printf, 1, 'iterate to convergence'
close, 1

end
