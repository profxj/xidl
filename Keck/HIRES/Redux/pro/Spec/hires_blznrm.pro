;+ 
; NAME:
; hires_blznrm
;     Version 1.1
;
; PURPOSE:
;  Extracts a spectrum from the trace flat, fits a function to the
;  shape and divides this into the data.  This is generally not
;  recommended as one often observes the flat through a filter that is
;  not used for science (e.g. ug5).
;
; CALLING SEQUENCE:
; hires_blznrm, hires, setup, obj_id, chip, /PLOT
;
; INPUTS:
;   hires    - HIRES structure
;   setup
;   obj_id  -  Object ID  (e.g. 1L)  (or STD index if /STD set)
;   chip    -  Blue (1), Green (2) OR Red (3) chip
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  /BSPLINE  -- Use a bspline to fit the blaze
;  /NONRM    -- Do not normalize the data.  Useful for simply
;               extracting the blaze
;  ZERO_FIL=  -- Array of file names that store the zero level
;                information.  This zero level is subtracted prior to
;                normalization
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
;   05-May-2006 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro hires_blznrm, hires, setup, obj_id, chip, PLOT=plot, BSPLINE=bspline, $
                  NONRM=nonrm, SV_BLZ=sv_blz, ZERO_FIL=zero_fil

  if  N_params() LT 3  then begin 
      print,'Syntax - ' + $
        'hires_blznrm, hires, setup, obj, [chip], /PLOT, ZERO_FIL=, /BSPLINE [v1.1]'
      return
  endif 

  if not keyword_set(NORD) then nord = 3L
  if not keyword_set(OBJ_NM) then obj_nm = 'a'
  if not keyword_set(CHIP) then chip = [1,2,3]
  nchip = n_elements(chip)

  for qq=0L,nchip-1 do begin
      ;; Dat fil
      allexp = where((hires.type EQ 'OBJ' OR hires.type EQ 'STD') AND $
                     hires.flg_anly NE 0 AND $
                     hires.setup EQ setup AND hires.obj_id EQ obj_id $
                     AND hires.chip EQ chip[qq])
      case chip[qq] of
          1: clrc = '_B' 
          2: clrc = '_G' 
          3: clrc = '_R' 
          else: stop
      endcase

      ;; Data
      if not keyword_set(NONRM) then begin
          subfil = strcompress(strtrim(hires[allexp[0]].Obj,2), $
                               /remove_all)+obj_nm
          datfil = 'FSpec/'+subfil+clrc+'.fits'
          spec = xmrdfits(datfil,1,/silent)
          nrmspec = spec
       endif

      if keyword_set(ZERO_FIL) then begin
         if strlen(zero_fil[qq]) GT 0 then begin
            zero_lvl = xmrdfits(zero_fil[qq]) 
            ;; Subtract
            print, 'hires_blznrm: Subtracting zero level: ', zero_fil[qq]
            spec.fx = spec.fx - zero_lvl
         endif
      endif

      ;; Flat
      flat = hires_getfil('qtz_fil', setup, CHIP=chip[qq])
      sz = size(flat, /dimensions)

      ;; Order structure
      ordr_str = hires_getfil('ordr_str', setup, chip=chip[qq])
      nordr = n_elements(ordr_str)

      ;; Save the blaze
      sv_blz = fltarr(nordr, sz[1])

      ;; Loop on orders
      for ii=0L,nordr-1 do begin
          oi = ordr_str[ii].order

          objcen = (ordr_str[ii].lhedg+ordr_str[ii].rhedg)/2.
          blz = extract_boxcar(flat, objcen, radius=2.)
          npix = n_elements(blz)

          ;; Save
          sv_blz[ii,0:npix-1] = blz
          if keyword_set(NONRM) then continue

          ;; Fit (smooth)
          if keyword_set(BSPLINE) then begin
              blz_set = bspline_iterfit(findgen(npix),blz,everyn=100, $
                                        /silent, yfit=blz_fit)
          endif else begin
              fitstr = x_setfitstrct(nord=nord)  ;; 3rd order POLY
              blz_fit = x_fitrej(findgen(npix),blz, FITSTR=fitstr)
          endelse

;;;;;;;Regina's (failed) attempt to get rid of the end wiggles:
;          blah = blz[where(blz NE 0.)]
;          npixblah = n_elements(blah)
;          blz_set = bspline_iterfit(findgen(npixblah),blah, nord=3,everyn=100, $
 ;                                  /silent, yfit=blz_fit)
          if keyword_set(PLOT) then x_splot, blz, ytwo=blz_fit, /bloc

          ;; Normalize
          pix = where(spec.wave[*,oi] GT 0.,npix)
          pix = pix[0] + lindgen(pix[npix-1]-pix[0]+1)

          ;; Scale the blaze
          nrm_blz = congrid(blz_fit, npix) 
          nrm_blz = nrm_blz/max(nrm_blz)
;          mwrfits, nrm_blz, 'Flats/Blz_B02_106.fits', /create

          ;; Write the order
          nrmspec.fx[pix,oi] = spec.fx[pix,oi]/nrm_blz
          nrmspec.var[pix,oi] = spec.var[pix,oi]/(nrm_blz^2)
      endfor

      ;; Write
      if not keyword_set(NONRM) then begin
          newout = 'FSpec/'+subfil+clrc+'bz.fits'
          mwrfits, nrmspec, newout, /create
          print, 'hires_nrm2dspec: Wrote ', newout
      endif
  endfor

  return
end
