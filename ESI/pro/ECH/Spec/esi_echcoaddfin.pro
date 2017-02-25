;+ 
; NAME:
; esi_echcoaddfin  
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
;  esi_echcoaddfin, esi, obj_id, CRVAL1=, CDELT=, NPIX=, /STD
;
; INPUTS:
;   esi   -  ESI structure
;   obj_id   -  Object ID  (e.g. 0L, 1L, etc)
;
; RETURNS:
;
; OUTPUTS:
;  One flux and one error array in 1D
;
; OPTIONAL KEYWORDS:
;   ORDRS=    - Orders to coadd (default: [0L,9L])
;   SPECFIL=  - Useful for fluxing files which resulted from a
;               combination of multiple nights
;   OUTROOT=  -- Output name
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   tst_coadd, esi, indx
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   13-Oct-2002 Written by JXP
;   04-Feb-2003 Polished (JXP)
;-
;------------------------------------------------------------------------------

pro esi_echcoaddfin, esi, obj_id, CRVAL1=crval1, CDELT=cdelt, NPIX=npix $
                     , STD=std, ORDRS=ordrs, SPECFIL=specfil, OUTROOT=outroot $
                     , OUTNM = outnm, obj_nm = obj_nm, OLDMSK = OLDMSK $
                     , NOVAR = NOVAR, SKY = SKY, CBIN = CBIN

  if  N_params() LT 2  AND NOT keyword_set( SPECFIL ) then begin 
      print,'Syntax - ' + $
        'esi_echcoaddfin, esi, obj_id, CRVAL1=, CDELT=, NPIX=, STD=, '
      print, '    ORDRS=   [v1.1]'
      return
  endif

;  Optional Keywords

  if not keyword_set(OBJ_NM) then obj_nm = 'a'
  if not keyword_set(CRVAL1) then crval1 = 3900.d
  if not keyword_set(CDELT) then cdelt = 0.00001447624d
  if not keyword_set(NPIX) then npix = 33000L
  if not keyword_set(ORDRS) then ordrs = [0L, 9L]

  if keyword_set( SPECFIL ) and not keyword_set( OUTROOT ) then begin
      print, 'esi_echcoaddfin: You must set OUTROOT too!!'
      return
  endif

; Grab exposure
  if keyword_set( ESI ) and not keyword_set( SPECFIL ) then begin
      if not keyword_set( STD ) then begin
          allexp = where(esi.type EQ 'OBJ' AND esi.flg_anly NE 0 AND $
                         esi.mode EQ 2 AND esi.obj_id EQ obj_id, nexp)
      endif else begin
          allexp = obj_id[0]
          nexp = 1
      endelse
      specfil = 'FSpec/'+strtrim(esi[allexp[0]].Obj,2)+obj_nm+'_ech.fits'
  endif else begin
      if keyword_set(obj_id) then allexp = obj_id[0] else allexp = 0
  endelse
  
  IF NOT KEYWORD_SET(CBIN) THEN BEGIN
     IF KEYWORD_SET(ESI) THEN BEGIN 
        binvec = esi[allexp].CBIN
        CBIN = total(binvec)/double(nexp)
     ENDIF ELSE CBIN = 1
  ENDIF
    
  x_wrechfspec, spec2d, specfil, /READ
  if spec2d.flg_flux EQ 0 then begin
      print, 'esi_echcoaddfin: Spectrum not fluxed! Returning...', specfil
      return
  endif
  bad = fltarr(10, 10, 2)
  IF KEYWORD_SET(OLDMSK) THEN BEGIN
      bad[0, 0, *] = [4287.4, 4292.6]
      bad[0, 1, *] = [4325, 4347.]
      bad[1, 0, *] = [4190, 4230.]
      bad[1, 1, *] = [4463, 4518.]
      bad[1, 2, *] = [4651, 4675.5]
      bad[2, 0, *] = [4440, 4480.]
      bad[2, 1, *] = [5015, 5044.]
      bad[3, 0, *] = [4711, 4800.5]
      bad[4, 0, *] = [5000, 5180.5]
      bad[4, 1, *] = [5180.5, 5450.5] ; Bad wave soln
      bad[5, 0, *] = [6560, 6600]
;  bad[5,1,*] = [6225, 6350]
      bad[7, 0, *] = [8160, 8300]
      bad[8, 0, *] = [8000, 8100]
  ENDIF ELSE BEGIN
      CASE CBIN OF 
          1.0: BEGIN            
              bad[0, 0, *] = [4330, 4400] ; UPPER EDGE
              bad[1, 0, *] = [4190, 4230] ; LOWER EDGE
              bad[1, 1, *] = [4645, 4685] ; HOTSPOT
              bad[2, 0, *] = [4440, 4500] ; LOWER EDGE
              bad[2, 1, *] = [5028, 5065] ; UPPER EDGE= + BAD COLUMN
              bad[3, 0, *] = [4675, 4860] ; LOWER EDGE
              bad[4, 0, *] = [5050, 5250] ; LOWER EDGE
              bad[5, 0, *] = [5600, 5710] ; LOWER EDGE
              bad[6, 0, *] = [6235, 6280] ; LOWER EDGE
          END
          2.0: BEGIN
              bad[0, 0, *] = [4325, 4347]      ; UPPER EDGE
              bad[1, 0, *] = [4190, 4230]      ; LOWER EDGE
              bad[1, 1, *] = [4651, 4676]      ; HOTSPOT
              bad[2, 0, *] = [4440, 4480]      ; LOWER EDGE
              bad[2, 1, *] = [5044, 5060] ; UPPER EDGE= + BAD COLUMN
              bad[3, 0, *] = [4675, 4831] ; LOWER EDGE
              bad[4, 0, *] = [5000, 5220] ; LOWER EDGE
              bad[5, 0, *] = [5600, 5700] ; LOWER EDGE
              bad[6, 0, *] = [6225, 6270.0] ; LOWER EDGE
          END
          ELSE: message, 'Problem with binning'
      ENDCASE
  ENDELSE
  
  tot_wave = 10^(alog10(CRVAL1) + dindgen(npix)*cdelt)
  weight = dblarr(npix)
  tot_flux = dblarr(npix)
  tot_novar = dblarr(npix)
  tot_sky   = dblarr(npix)
  sig = replicate(-1.D, npix)
  nosig = replicate(-1.D, npix)
; Loop

  for qq=ordrs[0],ordrs[1] do begin

      ;; Find first index
      indx = where(abs(tot_wave - spec2d.wave[0,qq]) LT 0.005, nindx)
      if nindx EQ 0 then stop
      indx = indx[0]

      ;; Good var
      a = where(spec2d.var[*,qq] GT 0, ngd)

      ;; Bad values
      if bad[qq,0,0] NE 0. then begin
          msk = replicate(1,ngd)
          for i=ordrs[0],ordrs[1] do begin
              b = where(spec2d.wave[a,qq] GT bad[qq,i,0] AND $
                        spec2d.wave[a,qq] LE bad[qq,i,1], nb)
              if nb NE 0 then msk[b] = 0
          endfor
          a = a[where(msk EQ 1, ngd)]
      endif

      ;; Add em in
      
      wtmp   = double(spec2d.var[a, qq])
      tmp    = tot_flux[indx+a]  + double(spec2d.fx[a, qq])/double(wtmp)
      novtmp = tot_novar[indx+a] + double(spec2d.novar[a, qq])/double(wtmp^2)
      stmp   = tot_sky[indx+a]   + double(spec2d.sky[a, qq])/double(wtmp)
      for i = 0L, ngd-1 do begin
          tot_flux[indx+a[i]] = tmp[i]
          tot_novar[indx+a[i]] = novtmp[i]
          tot_sky[indx + a[i]] = stmp[i]
          weight[indx+a[i]] = weight[indx+a[i]] + 1.D/double(wtmp[i])
      endfor
  endfor
  a = where(weight NE 0.)
  tot_flux[a] = tot_flux[a] / weight[a]
  tot_sky[a]  = tot_sky[a] / weight[a]
  tot_novar[a] = tot_novar[a]/weight[a]^2
  spec1d = tot_flux
  if arg_present(WAV) then wav = tot_wave
  sig[a] = sqrt(1./weight[a])
  nosig[a] = sqrt(tot_novar[a])
  ;; Header
  if keyword_set( ESI ) then $
    head = xheadfits(esi[allexp[0]].rootpth+esi[allexp[0]].img_root) $
  else mkhdr, head, spec1d

  sxaddpar, head, 'CRVAL1', alog10(crval1)
  sxaddpar, head, 'CDELT1', cdelt
  sxaddpar, head, 'CRPIX1', 1
  sxaddpar, head, 'CTYPE1', 'LINEAR'
  sxaddpar, head, 'DC-FLAG', 1
  sxaddpar, head, 'BITPIX', -32
  sxaddpar, head, 'NAXIS', 1
  sxaddpar, head, 'NAXIS1', n_elements(spec1d)
  sxdelpar, head, 'NAXIS2'
  sxdelpar, head, 'BZERO'
  sxdelpar, head, 'BSCALE'

  ;; Output
  if not keyword_set( OUTROOT ) then $
    outroot = 'FSpec/'+strtrim(esi[allexp[0]].Obj,2)+obj_nm

  outfil   = outroot + '_F.fits'
  sigfil   = outroot + '_E.fits'
  nosigfil = outroot + '_N.fits'
  skyfil   = outroot + '_S.fits'
;  outfil = 'FSpec/'+strtrim(esi[allexp[0]].Obj,2)+obj_nm+'_F.fits'
;  sigfil = 'FSpec/'+strtrim(esi[allexp[0]].Obj,2)+obj_nm+'_E.fits'
  print, 'esi_echcoaddfin: Writing ', outfil, ' ', sigfil
  mwrfits, spec1d, outfil, head, /create, /silent
  mwrfits, sig, sigfil, head, /create, /silent
  IF KEYWORD_SET(NOVAR) THEN mwrfits, nosig, nosigfil, head, /create, /silent
  IF KEYWORD_SET(SKY)   THEN mwrfits, tot_sky, skyfil, head, /create, /silent
  ;;
  print, 'esi_echcoaddfin: All done'
  return
end
