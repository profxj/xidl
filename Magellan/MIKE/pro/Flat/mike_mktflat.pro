;+ 
; NAME:
; mike_mktflat   
;     Version 2.0
;
; PURPOSE:
;    Combines all trace flats for a given setup into one final image.
;    The images are first bias subtracted using mike_subbias.
;    The images are median combined after scaling by the median after
;    iteratively rejecting bad pixels.  The image and its inverse
;    variance array are written to one fits file. Because of the
;    thermal gradients, it is highly recommended that one use 
;    only a coeval set of trace flats.  
;
; CALLING SEQUENCE:
;   
;  mike_mktflat, mike, setup, [side], /SVOV, /REDOOV, /CLOBBER
;
; INPUTS:
;   mike     -  MIKE structure
;   setup    -  Setup identifier 
;   [side]   -  Blue (1), Red (2), or both [1,2L]    (Default: [1,2L])
;
; RETURNS:
;
; OUTPUTS:
;  One normalized flat per setup and side
;  (e.g. 'Flats/Flat/_B_01_T.fits.gz')
;
; OPTIONAL KEYWORDS:
;   /CLOBBER - Overwrite the final fits file
;   /USEBIAS - Use the bias frame in bias subtraction
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   mike_mktflat, mike, 1.0
;
; PROCEDURES/FUNCTIONS CALLED:
;  mike_subbias
;  xcombine
;  mike_delov
;
; REVISION HISTORY:
;   17-Apr-2003 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro mike_mktflat, mike, setup, side, CLOBBER=clobber, USEBIAS=usebias, NO_MFLAT=no_mflat, $
                  _EXTRA=extra

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'mike_mktflat, mike, setup, [side], /CLOBBER, /NO_MFLAT, ' + $
        '/USEBIAS [v2.0]'
      return
  endif 
  

;  Optional Keywords
  if not keyword_set( SIDE ) then side = [1L,2L]
  
  for ii=0L,n_elements(side)-1 do begin
      qq = side[ii]
      if qq EQ 1 then print, 'mike_mktflat: Creating BLUE trace flat' $
      else print, 'mike_mktflat: Creating RED trace flat'

      ;; Check outfil
      outfil = mike_getfil('tflat_fil', setup, SIDE=qq, /name, CHKFIL=chkf)
      if CHKF NE 0  AND not keyword_set( CLOBBER ) then begin
          print, 'mike_mktflat: Trace flat exists, moving on..'
          continue
      endif

      ;; Grab the flats
      flt = where(mike.side EQ qq AND mike.flg_anly NE 0 AND $
                  strtrim(mike.type,2) EQ 'TFLT' $
                  AND mike.setup EQ setup, nflt)
      if nflt EQ 0 then begin
          print, 'mike_mktflat: No Flats of type TFLT with setup', $
            setup, ' and side ', qq
          print, 'Trying IFLT'
          flt = where(mike.side EQ qq AND mike.flg_anly NE 0 AND $
                  strtrim(mike.type,2) EQ 'IFLT' $
                  AND mike.setup EQ setup, nflt)
          if nflt EQ 0 then begin
             print, 'Nope no IFLTS either'
             print, 'mike_mktflat: Returning..'
             return
          endif
      endif


      gdflt = flt

      ;; Bias Subtract
      mike_subbias, mike, gdflt, /CLOBBER, USEBIAS=usebias, _EXTRA=extra

      ;; Median Combine
      if n_elements(gdflt) GT 1 then begin
          xcombine, 'OV/ov_'+mike[gdflt].img_root, img_flat, head, $
            FCOMB=2, SCALE='MED', GAIN=mike[gdflt[0]].gain, $
            RN=mike[gdflt[0]].readno
      endif else $
        img_flat = xmrdfits('OV/ov_'+mike[gdflt].img_root,0,head,/silent)

      ;; Flatten
      if not keyword_set(NO_MFLAT) then begin
         print, 'mike_mktflat: Imposing Milky Flat'
         mflat = mike_getfil('mflat_fil', setup, $
                             SIDE=qq, FIL_NM=flatfil, HEAD=fhead) 
         imflat = mike_getfil('mflat_fil', setup, INDX=1L, $
                              SIDE=qq, FIL_NM=flatfil, HEAD=fhead) 
         gmflt = where(mflat GT 0.5 AND finite(mflat) EQ 1 AND imflat GT 0., $
                       complement=badpix, ncomplement=nbad)
      endif else begin ;; NOT recommended
         szif = size(img_flat, /dimen)
         mflat = fltarr(szif[0], szif[1])
         mflat[*] = 1.
         gmflt = where(mflat GT 0.5)
         nbad = 0
      endelse


      ;; Inverse variance
      ivar = 1./((img_flat>1.)*mike[gdflt[0]].gain + mike[gdflt[0]].readno^2) $
                * (img_flat GT 1 AND mflat GT 0.5) 

      img_flat[gmflt] = img_flat[gmflt] / mflat[gmflt] 
      ivar[gmflt] = ivar[gmflt] * mflat[gmflt]^2  
      
      ;; Bad pix
      if nbad NE 0 then begin
          img_flat[badpix] = 0.
          ivar[badpix] = 0.
      endif

      ;; Output
      mike_taghead, head
      mwrfits, img_flat, outfil, head, /create, /silent
      mwrfits, ivar, outfil, /silent
      spawn, 'gzip -f '+outfil
      print, 'mike_mktflat: Flat created ', outfil+'.gz.'
      
      ;; DEL OV
      if not keyword_set(SVOV) then mike_delov, mike, gdflt
  endfor

  print, 'mike_mktflat: All done!'
  return
end
