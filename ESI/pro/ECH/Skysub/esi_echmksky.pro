;+ 
; NAME:
; esi_echmksky   
;     Version 1.0
;
; PURPOSE:
;    Create a Sky array from the sky fits
;
; CALLING SEQUENCE:
;   
;  esi_echmksky, esi, obj_id
;
; INPUTS:
;   esi   -  ESI structure
;   indx  -  Indices of objects to process
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   FLAT  - Flat file
;   BIAS  - Bias frame
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   esi_echmksky, esi 
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   19-Nov-2002 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro esi_echmksky, esi, obj_id, exp, BORDR=bordr, ALL_MNXWV=all_mnxwv

;
  if  N_params() LT 3  then begin 
      print,'Syntax - ' + $
        'esi_echmksky, esi, obj_id, exp, [v1.0]'
      return
  endif 
  
;  Optional Keywords
  if not keyword_set(BORDR) then bordr = 5L  ;; Order to begin bspline
  if not keyword_set(CRVAL1) then crval1 = 3900.d
  if not keyword_set(CDELT) then cdelt = 0.00001447624d
  if not keyword_set(NPIX) then npix = 28800L

  if not keyword_set( all_mnxwv ) then begin
      all_mnxwv = dblarr(10,2)
      all_mnxwv[0,*] = [4000., 4340.]
      all_mnxwv[1,*] = [4340., 4680.]
      all_mnxwv[2,*] = [4680., 5020.]
      all_mnxwv[3,*] = [5020., 5400.]
      all_mnxwv[4,*] = [5400., 5900.]
      all_mnxwv[5,*] = [5900., 6400.]
      all_mnxwv[6,*] = [6400., 7250.]
      all_mnxwv[7,*] = [7250., 8150.]
      all_mnxwv[8,*] = [8150., 9350.]
      all_mnxwv[9,*] = [9350., 10200.]
  endif

; Sky file
  indx = where(esi.flg_anly NE 0 AND esi.mode EQ 2 AND $
               esi.obj_id EQ obj_id AND $
               strtrim(esi.type,2) EQ 'OBJ', nindx)
  if nindx EQ 0 or exp GE nindx then begin
      print, 'esi_echmksky: No image to deal with!', obj_id, exp
      return
  endif
  skyfil = 'Sky/sky_'+esi[indx[exp]].img_root
  a = findfile(skyfil+'*',count=na)
  if na EQ 0 then begin
      print, 'esi_echmksky: No Sky file!  Returning...'
      return
  endif

; Outfil
  posi = strpos(esi[indx[exp]].img_root, '.fits')
  outfil = 'Sky/Sky_'+strmid(esi[indx[exp]].img_root,0,posi)+'F.fits'

; Create Final arrays
  fin_flux = fltarr(npix)
  fin_wv = 10^(alog10(CRVAL1) + dindgen(npix)*cdelt)
  ordr_i = lonarr(npix)

; Loop on pix grabbing order
  for i=0L,npix-1 do begin
      a = where(fin_wv[i] GT all_mnxwv[*,0] AND fin_wv[i] LE all_mnxwv[*,1],na)
      if na EQ 0 then ordr_i[i] = -1 else ordr_i[i] = a[0]
  endfor

; Loop on fits
  for qq=0L,9 do begin
      sky_fit = xmrdfits(skyfil, qq+1, /silent)
      gd = where(ordr_i EQ qq, ngd)
      if ngd EQ 0 then continue
      if qq LT BORDR then begin ;; Interpolate on POLY
          gdsky = where(sky_fit[*,0] GT 0. AND sky_fit[*,0] GT all_mnxwv[qq,0])
          linterp, sky_fit[gdsky,0], sky_fit[gdsky,1], fin_wv[gd], tmp
          fin_flux[gd] = tmp
      endif else begin
          fin_flux[gd] = bspline_valu(fin_wv[gd], sky_fit)
      endelse
  endfor

; Output
  print, 'esi_echmksky: Writing to '+outfil
  mkhdr, head, fin_flux
  sxaddpar, head, 'CRVAL1', alog10(crval1)
  sxaddpar, head, 'CDELT1', cdelt
  sxaddpar, head, 'CRPIX1', 1
  sxaddpar, head, 'CTYPE1', 'LINEAR'
  sxaddpar, head, 'DC-FLAG', 1
  mwrfits, fin_flux, outfil, head, /create, /silent
  
;  DONE
  print, 'esi_echmksky: All done! '
  return
end
