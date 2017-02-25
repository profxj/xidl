;+ 
; NAME:
; mike_coadd  
;     Version 1.1
;
; PURPOSE:
; Combine the individual 1D spectra from each order to create 1
; continuous spectrum.  The code attempts to match flux at the order
; ends and it is recommended that you have fluxed the obj prior to
; this step.
;
; CALLING SEQUENCE:
;   
;  mike_coadd, mike, obj_id, CRVAL1=, CDELT=, NPIX=, SPECFIL=, /STD
;
; INPUTS:
;   mike   -  MIKE structure
;   obj_id   -  Object ID  (e.g. 0L, 1L, etc)
;
; RETURNS:
;
; OUTPUTS:
;  One flux and one error array in 1D
;
; OPTIONAL KEYWORDS:
;    ORDRS=    - Orders to coadd (default: [0L,9L])
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   mike_coadd, mike, 2, CRVAL1=4000.0d, NPIX=50000, specfil='file.fits.gz'
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   13-Oct-2002 Written by JXP
;   04-Feb-2003 Polished (JXP)
;   12-Apr-2004 Updated for MIKE (RAS)
;-
;------------------------------------------------------------------------------

pro mike_coadd, mike, obj_id, side, CRVAL1=crval1, CDELT=cdelt, NPIX=npix,$
                     STD=std, ORDRS=ordrs, SPECFIL=specfil

  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'mike_coadd, mike, obj_id, side, CRVAL1=, CDELT=, NPIX=, SPECFIL=, STD=, '
      print, '    ORDRS=   [v1.1]'
      return
  endif

;  Optional Keywords

  if not keyword_set(OBJ_NM) then obj_nm = 'a'
  if not keyword_set(CRVAL1) then crval1 = 3000.d
  if not keyword_set(NPIX) then npix = 50000L

; Grab exposure
  if not keyword_set( STD ) then begin
      allexp = where(mike.type EQ 'OBJ' AND mike.flg_anly NE 0 AND $
                     mike.obj_id EQ obj_id, nexp)
  endif else begin
      allexp = obj_id[0]
      nexp = 1
  endelse

  if not keyword_set ( SIDE ) then begin
      side = mike[allexp[0]].side
  endif

  if not keyword_set (SPECFIL) then begin
      if (side EQ 1) then begin
          specfil = 'FSpec/'+ $
            strcompress(strtrim(mike[allexp[0]].Obj,2),/remove_all)+ $
            obj_nm+'_b.fits'
      endif else if (side EQ 2) then begin
          specfil = 'FSpec/'+ $
            strcompress(strtrim(mike[allexp[0]].Obj,2),/remove_all)+ $
            obj_nm+'_r.fits'
      endif else stop
  endif
  stop

  mike_wrfspec, spec2d, specfil, /READ
  if spec2d.flg_flux NE 1 then begin
      print, 'mike_coadd: Spectrum not fluxed! Be aware...', specfil
;      return
  endif

; Final wave array
; Enforce that wavelength array matches the output from mike_box,
; which has wavelength array starting at 3000A with 1.5 km/s pixels in
; the blue, and 2.1 km/s pixels in the red.  This way, the data will
; only have been rebinned once.

  if not keyword_set (crval1) then crval1=3000.0d
  velpix = (side EQ 1 ? 1.50d : 2.10d) * 2
  cdelt  = alog10(1.0d + velpix / 299792.458d)
  wave0  = alog10(3000.0d)

  ioffset = floor((alog10(crval1)-alog10(3000.0d))/(cdelt))
  tot_wave = 10^(wave0 + (dindgen(npix)+ioffset)*cdelt)

  weight = fltarr(npix)
  tot_flux = fltarr(npix)
  sig = replicate(-1.,npix)

  order_list = where(spec2d[0].phys_ordr NE 0, nord)
  ordrs[0] = (order_list[0] LT order_list[nord-1]) ? order_list[0] $
    : order_list[nord-1]

  ordrs[1] = (order_list[0] LT order_list[nord-1]) ? order_list[nord-1] $
    : order_list[0]

; Loop

  for qq=ordrs[0],ordrs[1] do begin

      ;; Find first index
      if (side EQ 1) then begin
          indx = where(abs((tot_wave-spec2d.wave[0,qq])/spec2d.wave[0,qq]*299792.458) LT 0.1d, nindx)
      endif else begin
          indx = where(abs((tot_wave-spec2d.wave[0,qq])/spec2d.wave[0,qq]*299792.458) LT 0.1d, nindx)
      endelse

      if (nindx NE 1) then stop
      indx = indx[0]

      ;; Good points
      a = where(spec2d.fx[*,qq] NE 0.00000, ngd)

      ;; Add em in
      
      wtmp = spec2d.var[a,qq]
      tmp = tot_flux[indx+a] + spec2d.fx[a,qq] / wtmp
      for i=0L,ngd-1 do begin
          tot_flux[indx+a[i]] = tmp[i]
          weight[indx+a[i]] = weight[indx+a[i]] + 1./wtmp[i]
      endfor
  endfor
  a = where(weight NE 0.)
  tot_flux[a] = tot_flux[a] / weight[a]
  spec1d = temporary(tot_flux)
  if arg_present(WAV) then wav = temporary(tot_wave)
  sig[a] = sqrt(1./weight[a])

  ;; Header
;  mkhdr, head, spec1d
  head = xheadfits(mike[allexp[0]].img_final)
  sxaddpar, head, 'CRVAL1', alog10(tot_wave[0])  
  sxaddpar, head, 'CDELT1', cdelt
  sxaddpar, head, 'CRPIX1', 1
  sxaddpar, head, 'CTYPE1', 'LINEAR'
  sxaddpar, head, 'DC-FLAG', 1
  sxaddpar, head, 'BITPIX', -32
  sxaddpar, head, 'NAXIS', 1
  sxaddpar, head, 'NAXIS1', n_elements(spec1d)
  sxdelpar, head, 'NAXIS2'

  ;; Output
  if (side EQ 1) then begin
      outfil = 'FSpec/'+strtrim(mike[allexp[0]].Obj,2)+obj_nm+'_Fb.fits'
      sigfil = 'FSpec/'+strtrim(mike[allexp[0]].Obj,2)+obj_nm+'_Eb.fits'
  endif else begin
      outfil = 'FSpec/'+strtrim(mike[allexp[0]].Obj,2)+obj_nm+'_Fr.fits'
      sigfil = 'FSpec/'+strtrim(mike[allexp[0]].Obj,2)+obj_nm+'_Er.fits'
  endelse

  mwrfits, spec1d, outfil, head, /create, /silent
  mwrfits, sig, sigfil, head, /create, /silent

  print, 'mike_coadd: All done!'
  return
end
