;+ 
; NAME:
; uves_mktflat   
;     Version 1.1
;
; PURPOSE:
;    Combines all flats for a given setup into one final image.
;    The images are first bias subtracted using uves_subbias.
;    The images are median combined after scaling by the median after
;    iteratively rejecting bad pixels.  The image and its inverse
;    variance array are written to one fits file. Because of the
;    thermal gradients, it is highly recommended that one use 
;    only a coeval set of trace flats.  
;
; CALLING SEQUENCE:
;   
;  uves_mktflat, uves, setup, [side], /SVOV, /REDOOV, /CLOBBER
;
; INPUTS:
;   uves     -  HIRES structure
;   setup    -  Setup identifier 
;   [side]   -  Blue (1), Green (2), Red (3), or multiple (Default:
;              [1,2,3L])
;
; RETURNS:
;
; OUTPUTS:
;  One normalized flat per setup and side
;  (e.g. 'Flats/Flat/_B_01_T.fits.gz')
;
; OPTIONAL KEYWORDS:
;   /REDOOV  - Overwrite OV files if they exist for the flats
;   /SVOV    - Save the OV files created during this step
;   /CLOBBER - Overwrite the final fits file
;   /USEBIAS - Use the bias frame in bias subtraction
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   uves_mktflat, uves, 1.0
;
; PROCEDURES/FUNCTIONS CALLED:
;  uves_subbias
;  xcombine
;  uves_delov
;
; REVISION HISTORY:
;   17-Apr-2003 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro uves_mktflat, uves, setup, side, REDOOV=redoov, SVOV=svov, $
                  CLOBBER=clobber, USEBIAS=usebias

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'uves_mktflat, uves, setup, [side], /SVOV, /REDOOV, /CLOBBER, ' + $
        '/USEBIAS [v1.1]'
      return
  endif 

  if not keyword_set(LMPNM) then lmpnm = 'FF'
  

;  Optional Keywords
  if not keyword_set( SIDE ) then side = [1L]
  
  for ii=0L,n_elements(side)-1 do begin
      qq = side[ii]
      case qq of 
          1: print, 'uves_mktflat: Creating BLUE flat'
          2: print, 'uves_mktflat: Creating RED flat'
          else: stop
      endcase


      ;; Grab the flats
      flt = where(uves.side EQ qq AND uves.flg_anly NE 0 AND $
                  strtrim(uves.type,2) EQ 'TFLT' and $
                  strmid(strtrim(uves.lamp,2),0,2) EQ LMPNM $
                  AND uves.setup EQ setup, nflt)
;                  strmid(strtrim(uves.lamp,2),0,3) EQ 'Deu' $

      ;; QA
      wcen = round(uves[flt[0]].xdangl)
      cwcen = strtrim(wcen,2)
      if setup LT 10 then c_s = '0'+strtrim(setup,2) $
      else c_s = strtrim(setup,2)
      fqa = 'QA/Flats'+cwcen+'_'+c_s
      a = findfile(fqa, count=count)
      if count EQ 0 then file_mkdir, fqa

      ;; Check outfil
      outfil = uves_getfil('qtz_fil', setup, WCEN=uves[flt[0]].xdangl, $
                           /name, CHKFIL=chkf)
      if CHKF NE 0  AND not keyword_set( CLOBBER ) then begin
          print, 'uves_mktflat: Flat exists, moving on..'
          continue
      endif

      if nflt EQ 0 then begin
          print, 'uves_mktflat: No Flats of type TFLT with setup', $
            setup, ' and side ', qq
          print, 'Trying TFLT'
          return
;          flt = where(uves.side EQ qq AND uves.flg_anly NE 0 AND $
;                  strtrim(uves.type,2) EQ 'TFLT' $
;                  AND uves.setup EQ setup, nflt)
;          if nflt EQ 0 then begin
;             print, 'Nope no TFLTs either'
;             print, 'uves_mktflat: Returning..'
;;             return
;          endif
      endif

      ;;
      gdflt = flt

      ;; Bias Subtract
      uves_subbias, uves, gdflt, USEBIAS=usebias, CLOBBER=clobber

      ;; Median Combine
      if n_elements(gdflt) GT 1 then begin
          ;; Filenames
          ovfil = ['']
          for q=0L,n_elements(gdflt)-1 do begin
              ovfil = [ovfil, $
                       uves_getfil('ov_fil', $
                                    OBJN=uves[gdflt[q]].img_root, /name)]
          endfor
          ovfil = ovfil[1:*]
          ;; Combine
          xcombine, ovfil, img_flat, head, $
            FCOMB=2, SCALE='MED', GAIN=uves[gdflt[0]].gain, $
            RN=uves[gdflt[0]].readno
      endif else $
        img_flat = xmrdfits(uves_getfil('ov_fil', OBJN=uves[gdflt[0]].img_root ))

      ;; Inverse variance
      ivar = 1./((img_flat>1.)*uves[gdflt[0]].gain + uves[gdflt[0]].readno^2) $
                * (img_flat GT 1 )


      ;; Zero out edge of blue side
;      if not keyword_set( ALLBLUE ) and qq EQ 1 then begin
;          lrow = round(1970./uves[gdflt[0]].colbin) 
;          ivar[lrow:*,*] = -1
;      endif

      ;; Output
      mwrfits, img_flat, outfil, head, /create, /silent
      mwrfits, ivar, outfil, /silent
      spawn, 'gzip -f '+outfil
      print, 'uves_mktflat: Flat created ', outfil+'.gz.'
      
      ;; DEL OV
      if not keyword_set(SVOV) then uves_delov, uves, gdflt, /silent
  endfor

  print, 'uves_mktflat: All done!'
  return
end
