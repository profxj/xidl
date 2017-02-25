;+ 
; NAME:
; apf_mktflat   
;     Version 1.1
;
; PURPOSE:
;    Combines all flats for a given setup into one final image.
;    The images are first bias subtracted using apf_subbias.
;    The images are median combined after scaling by the median after
;    iteratively rejecting bad pixels.  The image and its inverse
;    variance array are written to one fits file. Because of the
;    thermal gradients, it is highly recommended that one use 
;    only a coeval set of trace flats.  
;
; CALLING SEQUENCE:
;   
;  apf_mktflat, apf, setup, [chip], /REDOOV, /CLOBBER
;
; INPUTS:
;   apf     -  HIRES structure
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
;   apf_mktflat, apf, 1.0
;
; PROCEDURES/FUNCTIONS CALLED:
;  apf_subbias
;  xcombine
;  apf_delov
;
; REVISION HISTORY:
;   17-Apr-2003 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro apf_mktflat, apf, setup, CLOBBER=clobber, USEBIAS=usebias, AVG=avg

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'apf_mktflat, apf, setup, [chip], /CLOBBER, ' + $
        '/USEBIAS [v1.1]'
      return
  endif 
  
  ;; QA
  if setup LT 10 then fqa = 'QA/Flats0'+strtrim(setup,2) $
  else fqa = 'QA/Flats'+strtrim(setup,2)
  a = findfile(fqa, count=count)
  if count EQ 0 then file_mkdir, fqa

;  Optional Keywords
  print, 'apf_mktflat: Creating Single flat'

  ;; Check outfil
  outfil = apf_getfil('qtz_fil', setup, /name, CHKFIL=chkf)
  if CHKF NE 0  AND not keyword_set( CLOBBER ) then begin
     print, 'apf_mktflat: Flat exists, moving on..'
     return
  endif

  ;; Grab the flats
  flt = where(apf.flg_anly NE 0 AND $
              strtrim(apf.type,2) EQ 'TFLT' $
              AND apf.setup EQ setup, nflt)
  if nflt EQ 0 then begin
     print, 'apf_mktflat: No Flats of type TFLT with setup', setup
     print, 'Trying TFLT'
     return
  endif

  ;;
  gdflt = flt

  ;; Bias Subtract
  apf_subbias, apf, gdflt, USEBIAS=usebias, CLOBBER=clobber

  ;; Median Combine
  if n_elements(gdflt) GT 1 then begin
     ;; Filenames
     ovfil = ['']
     for q=0L,n_elements(gdflt)-1 do begin
        ovfil = [ovfil, $
                 apf_getfil('ov_fil', FRAME=apf[gdflt[q]].frame, /name)]
     endfor
     ovfil = ovfil[1:*]
     ;; Median or average?
     if keyword_set(AVG) then fcomb = 3 else fcomb = 2
     ;; Combine
     xcombine, ovfil, img_flat, head, $
               FCOMB=fcomb, SCALE='MED', GAIN=apf[gdflt[0]].gain, $
               RN=apf[gdflt[0]].readno
  endif else begin
     nm = apf_getfil('ov_fil', FRAME=apf[gdflt[0]].frame,/name)
     img_flat = xmrdfits(apf_getfil('ov_fil', FRAME=apf[gdflt[0]].frame,/name), 0, head)
  endelse
  
  ;; Inverse variance
  ivar = 1./((img_flat>1.)*apf[gdflt[0]].gain + apf[gdflt[0]].readno^2) $
         * (img_flat GT 1 )
  
  ;; Output
  mwrfits, img_flat, outfil, head, /create, /silent
  mwrfits, ivar, outfil, /silent
  spawn, 'gzip -f '+outfil
  print, 'apf_mktflat: Flat created ', outfil+'.gz.'
  
  ;; DEL OV
  if not keyword_set(SVOV) then apf_delov, apf, gdflt, /silent

  print, 'apf_mktflat: All done!'
  return
end
