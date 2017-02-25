;+ 
; NAME:
; x_pltvpfit
;  V1.2
;
; PURPOSE:
;    Given a ISM structre write the files out
;
; CALLING SEQUENCE:
;  flux = x_pltvpfit(wrest, wv, vpfils, FWHM=)
;
; INPUTS:
;  wrest -- Rest wavelength [Ang]
;  wv    -- Wavelength array
;  vpfils -- Structure with VPFIT info
;
; RETURNS:
;
; OUTPUTS:
;  Ion structure.  Column densities have linear values
;
; OPTIONAL KEYWORDS:
;  FWHM=  -- FWHM in pixels of the spectral resolution
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   05-Apr-2007 Written by JXP 
;-
;------------------------------------------------------------------------------
function x_pltvpfit, wrest, wv, vpfils, FWHM=fwhm, LINES=sv_lin

; parse_ismlst -- Reads in DLA data to a structure

  if (N_params() LT 2) then begin 
    print,'Syntax - ' + $
             'pltvp = ism_pltvpfit(wrest, wv, vpfils) [v1.0]'
    return, -1
  endif 
  if not keyword_set(FWHM) then fwhm = 4.

  ;; Setup
  npix = n_elements(wv)
  fin_fx = replicate(1.,npix)
  wmn = min(wv,max=wmx)


  ;; Bookkeeping
  getion, wrest, ion, FNM=nm
  lin = x_setline(wrest)

  nfil = n_elements(vpfils)
  flgsv = 0
  for kk=0L,nfil-1 do begin
      ;; Parse the VPFIT file
      prs_vpfit, vpfils[kk], strct, /silent
      ncomp = n_elements(strct)

      mt = where((1+strct.zabs)*wrest GT wmn and (1+strct.zabs)*wrest LT wmx AND $
                 strmatch(strct.ion,nm) and strlen(strct.ion) EQ strlen(nm), nmt)
      if nmt EQ 0 then continue

      ;; Create lines
      all_lin = replicate(lin,nmt)
      
      ;; Fill it up
      for qq=0L,nmt-1 do begin
          ;; Values
          all_lin[qq].zabs = strct[mt[qq]].zabs
          all_lin[qq].b = strct[mt[qq]].b
          all_lin[qq].N = strct[mt[qq]].N
      endfor

      ;; Create
      fx = x_voigt(wv, all_lin, FWHM=fwhm)
      fin_fx = fin_fx * fx
      if flgsv EQ 0 then begin
          sv_lin = all_lin 
          flgsv = 1
      endif else sv_lin = [sv_lin,all_lin]
  endfor

  return, fin_fx
end
