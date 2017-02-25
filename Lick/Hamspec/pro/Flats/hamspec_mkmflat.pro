;+ 
; NAME:
; hamspec_mkmflat   
;     Version 1.1
;
; PURPOSE:
;    This set of routines takes a series of flats observed through the
;    diffuser or out of focus and creates a normalized Flat used to correct
;    pixel-to-pixel response variations.  
;
;    hamspec_mkmflat :: The main routine simply does some basic
;    accounting, organizes bias subtraction and performs I/O.
;    It now also stacks the images and does the normalization.  The
;    main driver is a bspline routine.
;
; CALLING SEQUENCE:
;   
;  hamspec_mkmflat, hamspec, setup, [side]
;
; INPUTS:
;   hamspec     -  MIKE structure
;   setup    -  Setup identifier 
;
; RETURNS:
;
; OUTPUTS:
;  One normalized flat per setup per side with names like
;  'Flats/Flat_1.fits.gz' 
;
; OPTIONAL KEYWORDS:
;   /CLOBBER - Overwrite Output MilkyFlat
;   /USEBIAS - Use bias frame in OV subtraction
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   hamspec_mkmflat, hamspec, 1
;
; PROCEDURES/FUNCTIONS CALLED:
;   hamspec_getfil
;   hamspec_subbias
;   hamspec_mkflat_work
;   x_statarray
;
; REVISION HISTORY:
;   16-May-2003 Adapted by JXP from existing programs by SB
;   24-Feb-2004 Switched to a series of median/linear interpolations (SB)
;   22-Jun-2005 Returned to simple median smooth as the default
;-
;------------------------------------------------------------------------------
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro hamspec_mkmflat, hamspec, setup, CLOBBER=clobber, USEBIAS=usebias, $
                  SMOOTH=smooth, _EXTRA=extra

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'hamspec_mkmflat, hamspec, setup, ' + $
        ' /USEBIAS, /CLOBBER, /SMOOTH [v2.0]'
      return
  endif 
  
  ;; QA
  if setup LT 10 then fqa = 'QA/Flats0'+strtrim(setup,2) $
  else fqa = 'QA/Flats'+strtrim(setup,2)
  a = findfile(fqa, count=count)
  if count EQ 0 then file_mkdir, fqa

;  Optional Keywords
  ;; Outfil
  outfil = 'Flats/Flat_'+strtrim(setup,2)+'_M.fits'
  a= findfile(outfil+'*', count=CHKFIL)
  if CHKFIL NE 0  AND not keyword_set( CLOBBER ) then begin
     print, 'hamspec_mkmflat: Flat exists, moving on..'
     return
  endif
     
  
  ;; Will need to loop over filters
  amflt = where(hamspec.flg_anly NE 0 AND $
                strtrim(hamspec.type,2) EQ 'MFLT' AND $
                hamspec.setup EQ setup, nflt)
  if nflt EQ 0 then begin
     print, 'hamspec_mkmflat: No Flats of type MFLT found!' 
     return
  endif
  unib = hamspec[amflt[uniq(hamspec[amflt].block,sort(hamspec[amflt].block))]].block
  nuni = n_elements(unib)
  ;nuni = -1L
  for ii=0L,nuni-1 do begin
     print, 'hamspec_mkmflat: Creating Milky flat for the filter', unib[ii]
     gdflt = where(hamspec.flg_anly NE 0 AND $
                   strtrim(hamspec.type,2) EQ 'MFLT' AND $
                   strmatch(strtrim(hamspec.block,2), unib[ii]) and $
                   hamspec.setup EQ setup, nflt)
     
     ;; Bias Subtract
     hamspec_subbias, hamspec, gdflt, CLOBBER=ovclob, USEBIAS=usebias, $
                      _EXTRA=extra
     
     ;; Read the images
     i0 = xmrdfits(hamspec[gdflt[0]].img_ov, /silen)
     sz = size(i0, /dime)
     
     flat = fltarr(sz[0], sz[1], nflt)
     flat[*,*,0] = i0
     for jj=1L,nflt-1 do flat[*,*,jj] = xmrdfits(hamspec[gdflt[jj]].img_ov, /silen)
     
     ;; Stack the images
     print, 'hamspec_mkmflat: Stacking...'
     x_statarray, flat, 3, mean=uni_final, stddev=flatstddev, sigrej=2., $
                  /OVRIDE
     if ii EQ 0 then begin
        all_flat = fltarr(sz[0], sz[1], nuni)
        all_stddv = fltarr(sz[0], sz[1], nuni)
     endif
     all_flat[*,*,ii] = uni_final
     all_stddv[*,*,ii] = flatstddev
  endfor

  ;; Combine multiple filters
  if nuni GT 1 then begin
     milky = total(all_flat,3) ;; Straight sum
     sig_milky = sqrt(total(all_stddv^2,3))
  endif else milky = all_flat[*,*,0]

  norm = fltarr(sz[0], sz[1])
  norm[*] = 1.
  for kk=0L,sz[0]-1 do begin
     if (kk MOD 50) EQ 0 then print, 'kk', kk
     dumh = (milky[kk,*])[*]
     medd = median(dumh, 62L)
     nrm = dumh / medd
     ;; Bad points
     bad = where(abs(nrm-1) GT 0.1, nbad)
     ivar = 1./dumh
     if nbad GT 0 then ivar[bad] = -1.
     
     ;; Bspline
     bfit = bspline_iterfit(findgen(sz[1]), dumh, invvar=ivar, $
                            everyn=32L, yfit=yval, upper=5.0, lower=5.0, $
                            outmask=mask, /silen)
     norm[kk,*] = yval
  endfor
  ;;norm = transpose( x_medianrow(image,width) )
  flat = milky/norm
  sig_flat = sig_milky/norm

  var = sig_flat^2 
  ivar = 1./(var + (var EQ 0)) * (sig_flat GT 0)

  ;; Output
  mwrfits, flat, outfil, /create, /silent
  mwrfits, ivar, outfil, /silent
  spawn, 'gzip -f '+outfil
  print, 'hamspec_mkmflat: Flat created ', outfil+'.gz'
  
  if not keyword_set(SVOV) then hamspec_delov, hamspec, gdflt


  print, 'hamspec_mkmflat: All done!'
  return
end
