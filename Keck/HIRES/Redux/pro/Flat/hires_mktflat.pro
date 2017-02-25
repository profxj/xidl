;+ 
; NAME:
; hires_mktflat   
;     Version 1.1
;
; PURPOSE:
;    Combines all flats for a given setup into one final image.
;    The images are first bias subtracted using hires_subbias.
;    The images are median combined after scaling by the median after
;    iteratively rejecting bad pixels.  The image and its inverse
;    variance array are written to one fits file. Because of the
;    thermal gradients, it is highly recommended that one use 
;    only a coeval set of trace flats.  
;
; CALLING SEQUENCE:
;   
;  hires_mktflat, hires, setup, [chip], /REDOOV, /CLOBBER
;
; INPUTS:
;   hires     -  HIRES structure
;   setup    -  Setup identifier 
;   [chip]   -  Blue (1), Green (2), Red (3), or multiple (Default:
;              [1,2,3L])
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
;   hires_mktflat, hires, 1.0
;
; PROCEDURES/FUNCTIONS CALLED:
;  hires_subbias
;  xcombine
;  hires_delov
;
; REVISION HISTORY:
;   17-Apr-2003 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro hires_mktflat, hires, setup, chip, CLOBBER=clobber, USEBIAS=usebias, $
                   AVG=avg

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'hires_mktflat, hires, setup, [chip], /CLOBBER, ' + $
        '/USEBIAS [v1.1]'
      return
  endif 
  
  ;; QA
  if setup LT 10 then fqa = 'QA/Flats0'+strtrim(setup,2) $
  else fqa = 'QA/Flats'+strtrim(setup,2)
  a = findfile(fqa, count=count)
  if count EQ 0 then file_mkdir, fqa

;  Optional Keywords
  if not keyword_set( CHIP ) then chip = [1L,2L,3L]
  
  for ii=0L,n_elements(chip)-1 do begin
      qq = chip[ii]
      case qq of 
          -1: print, 'hires_mktflat: Creating Single flat'
          1: print, 'hires_mktflat: Creating BLUE flat'
          2: print, 'hires_mktflat: Creating GREEN flat'
          3: print, 'hires_mktflat: Creating RED flat'
          else: stop
      endcase

      ;; Check outfil
      outfil = hires_getfil('qtz_fil', setup, CHIP=qq, /name, CHKFIL=chkf)
      if CHKF NE 0  AND not keyword_set( CLOBBER ) then begin
          print, 'hires_mktflat: Flat exists, moving on..'
          continue
      endif

      ;; Grab the flats
      flt = where(hires.chip EQ qq AND hires.flg_anly NE 0 AND $
                  strtrim(hires.type,2) EQ 'TFLT' $
                  AND hires.setup EQ setup, nflt)
      if nflt EQ 0 then begin
          print, 'hires_mktflat: No Flats of type TFLT with setup', $
            setup, ' and chip ', qq
          print, 'Trying TFLT'
          return
;          flt = where(hires.chip EQ qq AND hires.flg_anly NE 0 AND $
;                  strtrim(hires.type,2) EQ 'TFLT' $
;                  AND hires.setup EQ setup, nflt)
;          if nflt EQ 0 then begin
;             print, 'Nope no TFLTs either'
;             print, 'hires_mktflat: Returning..'
;;             return
;          endif
      endif

      ;;
      gdflt = flt

      ;; Bias Subtract
      hires_subbias, hires, gdflt, USEBIAS=usebias, CLOBBER=clobber

      ;; Median Combine
      if n_elements(gdflt) GT 1 then begin
          ;; Filenames
          ovfil = ['']
          for q=0L,n_elements(gdflt)-1 do begin
              ovfil = [ovfil, $
                       hires_getfil('ov_fil', $
                                    CHIP=hires[gdflt[q]].chip, $
                                    FRAME=hires[gdflt[q]].frame, /name)]
          endfor
          ovfil = ovfil[1:*]
          ;; Median or average?
          if keyword_set(AVG) then fcomb = 3 else fcomb = 2
          ;; Combine
          xcombine, ovfil, img_flat, head, $
            FCOMB=fcomb, SCALE='MED', GAIN=hires[gdflt[0]].gain, $
            RN=hires[gdflt[0]].readno
      endif else $
        img_flat = xmrdfits(hires_getfil('ov_fil', $
                                    CHIP=hires[gdflt[0]].chip, $
                                    FRAME=hires[gdflt[0]].frame,/name), 0, head)

      ;; Inverse variance
      ivar = 1./((img_flat>1.)*hires[gdflt[0]].gain + hires[gdflt[0]].readno^2) $
                * (img_flat GT 1 )


      if qq EQ -1 then begin
         hires_badpix_single, badpix
         ;; Rebin
         sz_img = size(ivar,/dimen)
         msk = rebin(badpix, sz_img[0], sz_img[1])
         bad = where(msk GT 0.)
         ivar[bad] = 0.
      endif

      ;; Zero out edge of blue chip
      if not keyword_set( ALLBLUE ) and qq EQ 1 then begin
          lrow = round(1970./hires[gdflt[0]].colbin) 
          ivar[lrow:*,*] = -1
      endif

      ;; Output
      mwrfits, img_flat, outfil, head, /create, /silent
      mwrfits, ivar, outfil, /silent
      spawn, 'gzip -f '+outfil
      print, 'hires_mktflat: Flat created ', outfil+'.gz.'
      
      ;; DEL OV
      if not keyword_set(SVOV) then hires_delov, hires, gdflt, /silent
  endfor

  print, 'hires_mktflat: All done!'
  return
end
