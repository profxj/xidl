;+ 
; NAME:
; x_fndflux   
;    Version 1.1
;
; PURPOSE:
;  OBSOLETE:  Use x_calibstd
;
; CALLING SEQUENCE:
;   x_fndflux, wv, fx, star, fitstr
;
; INPUTS:
;   wv   - Wavelength array
;   fx   - Stellar flux
;   star - Specphot standard ('LTT7379')
;   fitstr - Fit of the star
;
; RETURNS:
;
; OUTPUTS:
;   fitstr    - Fit structure
;
; OPTIONAL KEYWORDS:
;   INTER - Interactive
;   
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   x_fndflux, wv, fx, star, fitstr
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   15-Apr-2002 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro x_fndflux, wv, fx, star, fitstr, INTER=inter, LSTROOT=lstroot, EXTRAP=extrap
;
  if  N_params() LT 4  then begin 
    print,'Syntax - ' + $
      'x_fndflux, wv, fx, star, fitstr, /INTER, LSTROOT=, /EXTRAP [v1.0]'
    return
  endif 

;  Optional Keywords

  if not keyword_set( LSTROOT ) then $
    lstroot = getenv('XIDL_DIR')+'/Spec/Flux/Lists/'

; Read in Star file

  case star of 
      'LTT7379': strfil = lstroot+'LTT7379.lst'
      'BD+284211': strfil = lstroot+'BD284211.lst'
      'G191B2B': strfil = lstroot+'G191B2B.lst'
      else : begin
          print, 'x_fndflux: Star file not found'
          return
      end
  endcase

  readcol, strfil, str_wv, str_dwv, str_mag
  nstr = n_elements(str_wv)

; Fit structure
  if not keyword_set( FITSTR ) then begin
      fitstr = { fitstrct }
      fitstr.func = 'BSPLIN'
      fitstr.flg_rej = 0 
      fitstr.nord = (nstr/2) < 30L
  endif

; Extapolate
  if keyword_set( EXTRAP ) then begin
      add_wv = (lindgen(extrap)+1)*100. + str_wv[nstr-1]
      add_dwv = replicate(str_dwv[nstr-1], extrap)
      add_mag = replicate(str_mag[nstr-1], extrap)
      nstr = nstr + extrap
      str_wv = [str_wv, add_wv]
      str_dwv = [str_dwv, add_dwv]
      str_mag = [str_mag, add_mag]
  endif

; Calculate observed fluxes

  ; Add up e^-/Ang

  obj_fx = fltarr(nstr)
  for j=0L,nstr-1 do begin
      sum = 0.
      ;; Dlambda
      dlamb = str_dwv[j]
      ;; Find left edge
      lbin = str_wv[j]-dlamb/2.
      lft = (where(wv LT lbin AND shift(wv,-1) GE lbin))[0]
      if lft EQ -1 then continue  ; Bin extends beyond data
      cpix = total(wv[lft:lft+1])/2.
      if lbin GT cpix then begin
          sum = sum + ((1.-((lbin-cpix)/(wv[lft+1]-cpix)))/2. + 0.5)*fx[lft+1]
          lft = lft + 1
      endif else sum = sum + (cpix - lbin)/(cpix-wv[lft])*fx[lft]/2.
      lft = lft+1  ; Update left for whole pixels
      ;; Right edge
      rbin = str_wv[j]+dlamb/2.
      rgt = (where(wv LT rbin AND shift(wv,-1) GE rbin))[0]
      if rgt EQ -1 then continue  ; Bin extends beyond data
      cpix = total(wv[rgt:rgt+1])/2.
      if rbin LT cpix then begin
          sum = sum + ((1.-((cpix-rbin)/(cpix-wv[rgt])))/2. + 0.5)*fx[rgt]
          rgt = rgt - 1
      endif else sum = sum + (rbin-cpix)/(wv[rgt+1]-cpix)*fx[rgt+1]/2.
      rgt = rgt-1  ; Update left for whole pixels
      ;; Whole pixels
      sum = sum + total(fx[lft:rgt])
      ;; Save
      obj_fx[j] = sum/dlamb  ; e- per Ang
  endfor

  ;; Good values
  gdval = where(obj_fx GT 0., ngd)

  ;; Now fnu
  str_fnu = 10^(-(str_mag[gdval] + 48.595)/2.5)  ; Standard formula
  obj_fnu = obj_fx[gdval] * str_wv[gdval]^2 / (3e18)  ; c in Ang (dlamb already done)

  ;; Take ratio
  rto = str_fnu / obj_fnu

  ;; Reset nord
  fitstr.nord = fitstr.nord < (ngd-3)
  ;; Fit log
  fit = x1dfit(str_wv[gdval], alog10(rto), FITSTR=fitstr, /INTER)


;  conv = x_calcfit(wv, FITSTR=fitstr)
;  tru = 10.^conv

;  dwv = (shift(wv,-1)-shift(wv, 1))/2.
;  dwv[0] = wv[1]-wv[0]
;  npix = n_elements(wv)
;  dwv[npix-1] = wv[npix-1]-wv[npix-2]

;  new_fnu = (fx/dwv) * wv^2 / (3e18) * tru
;
;  plot, wv, new_fnu
;  oplot, str_wv[gdval], str_fnu, psym=1
;  stop

  return
end
  
