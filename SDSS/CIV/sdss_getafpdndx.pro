;+ 
; NAME:
; sdss_getafpdndx.pro
;    Version 1.0
;
; PURPOSE:
;   Explicitly query user bias file to calculate accepted
;   false-positive dN/dX (as in sdss_calcdndx() structure. Typically
;   called by sdss_functions.pro but couldn't live there
;   because calls sdss_functions.pro which in turn call it.
;
; CALLING SEQUENCE:
;   
;   afpdndx = sdss_getafpdndx(zlim,ewlim)
;
; INPUTS: 
;
; RETURNS: 
;   dN/dX structure for accepted false-positives, from sdss_calcdndx()
;
; OUTPUTS: 
;
; OPTIONAL KEYWORDS:
;   _extra= gets passed to sdss_calcdndx(), which passes it to sdss_getdxzw()
;
; OPTIONAL OUTPUTS:
;
; COMMENTS: 
;    This is unfortunately a bit hard-coded.
;
; EXAMPLES:
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   11-Apr-2012  created by KLC
;-
;---------------------------------------------------------------------------  
function sdss_getafpdndx, zlim, ewlim, biasuser_fil=biasuser_fil, $
                          user_indx=user_indx, _extra=extra
  if n_params() ne 2 then begin
     print,'Syntax - sdss_getafpdndx(zlim, ewlim, [biasuser_fil=, user_indx=, _extra=])'
     return,-1
  endif

  sdssdir = sdss_getsdssdir()

  if not keyword_set(biasuser_fil) then $
     biasuser_fil = sdssdir+'mcfkiv/cmpltrec_userbias_noBAL.fit' 
  if n_elements(ewlim) eq 1 then ewrng = [ewlim,!values.f_infinity] $
  else ewrng = ewlim

  userbiasstr = xmrdfits(biasuser_fil,1,/silent) 
  fpstr = xmrdfits(biasuser_fil,3,/silent)

  ;; Find the best binning to use for the given redshift range
  user_indx = sdss_getuserbiasindex(userbiasstr,zlim)

  ;; Must work self-consistently with the redshift bin
  sub = where(fpstr[userbiasstr.afp_indx].zabs_orig[0] ge $
              userbiasstr.zlim[user_indx,0] $
              and fpstr[userbiasstr.afp_indx].zabs_orig[0] lt $
              userbiasstr.zlim[user_indx,1],nsub)
  sub = userbiasstr.afp_indx[sub] 

  ;; output will be for only one bin b/c of how completeness
  ;; structure designed
  ;; _extra includes /dz, and stuff for sdss_getdxw()
  afpdndx = $
     sdss_calcdndx(fpstr[sub], sdssdir+userbiasstr.cmplt_fil[user_indx],$
                   ewrng, dblt_name=userbiasstr.dblt_name,ewmax=3.5,/silent,$
                   final=0, afp_flg=0, _extra=extra) 
  print,'sdss_getafpdndx():',zlim[0],zlim[1],ewrng[0],ewrng[1],$
        user_indx,afpdndx.dndx[0],$
        format='(a,1x,2(f5.2,1x),1x,2(f6.2,1x),1x,i1,2x,f6.4)'
  return, afpdndx
end                             ; sdss_getafpdndx()


