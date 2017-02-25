;+ 
; NAME:
; hamspec_mkflats   
;     Version 1.1
;
; PURPOSE:
;    Combines all flats for a given setup into one final image.
;    The images are first bias subtracted using hamspec_subbias.
;    The images are median combined after scaling by the median after
;    iteratively rejecting bad pixels.  The image and its inverse
;    variance array are written to one fits file. Because of the
;    thermal gradients, it is highly recommended that one use 
;    only a coeval set of trace flats.  
;
; CALLING SEQUENCE:
;   
;  hamspec_mktflat, hamspec, setup, /REDOOV, /CLOBBER
;
; INPUTS:
;   hamspec     -  HIRES structure
;   setup    -  Setup identifier 
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
;   hamspec_mktflat, hamspec, 1.0
;
; PROCEDURES/FUNCTIONS CALLED:
;  hamspec_subbias
;  xcombine
;  hamspec_delov
;
; REVISION HISTORY:
;   17-Apr-2003 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro hamspec_mkflats, hamspec, setup, type, CLOBBER=clobber, USEBIAS=usebias, $
                   AVG=avg

;
  if  N_params() LT 3  then begin 
      print,'Syntax - ' + $
        'hamspec_mkflats, hamspec, setup, type, /CLOBBER, ' + $
        '/USEBIAS [v1.1]'
      return
  endif 

  if type NE 'TFLT' and type NE 'PFLT' then begin
     print, 'hamspec_mkflats: Wrong flat type!', type
     return
  endif

  ;; QA
  if setup LT 10 then fqa = 'QA/Flats0'+strtrim(setup,2) $
  else fqa = 'QA/Flats'+strtrim(setup,2)
  a = findfile(fqa, count=count)
  if count EQ 0 then file_mkdir, fqa

  ;; Check outfil
  case type of 
     'TFLT': outfil = hamspec_getfil('qtz_fil', setup, /name, CHKFIL=chkf)
     'PFLT': outfil = hamspec_getfil('pixflt_fil', setup, /name, CHKFIL=chkf)
     else: stop
  endcase
  if CHKF NE 0  AND not keyword_set( CLOBBER ) then begin
     print, 'hamspec_mkflats: Flat exists, moving on..'
     return
  endif

  ;; Grab the flats
  flt = where(hamspec.flg_anly NE 0 AND $
              strtrim(hamspec.type,2) EQ type $
              AND hamspec.setup EQ setup, nflt)
  if nflt EQ 0 then begin
     print, 'hamspec_mkflats: No Flats of type ', type, ' with setup',  setup
     return
  endif
  
  ;;
  gdflt = flt
  
  ;; Bias Subtract  -- Get and set gain ratio if needed
  if hamspec[0].amp EQ 2 and type EQ 'PFLT' then begin
     print,'hamspec_mkflats:  Will calculate gain ratio'
  endif else gratio = hamspec[flt[0]].ratio_gain

  hamspec_subbias, hamspec, gdflt, USEBIAS=usebias, CLOBBER=clobber, $
                   GAIN_RATIO=gratio
  if hamspec[0].amp EQ 2 and type EQ 'PFLT' then begin
     print,'hamspec_mkflats:  Setting gain ratio to ', gratio
     allset = where(hamspec.setup EQ setup)
     hamspec[allset].ratio_gain = gratio
  endif
                   

  ;; Loop on Filters
  filts = hamspec[gdflt[uniq(hamspec[gdflt].block, sort(hamspec[gdflt].block))]].block
  nfilts = n_elements(filts)
;  ovfil0 = hamspec_getfil('ov_fil', $
 ;                FRAME=hamspec[gdflt[0]].frame, /name)]
 ; head0 = xheadfits(ovfil0)
 ; sv_stack = fltarr(sxpar(head0,'NAXIS1'), sxpar(head0,'NAXIS2'), nfilts)
;  sv_stack_iv = fltarr(sxpar(head0,'NAXIS1'), sxpar(head0,'NAXIS2'), nfilts)

  for ss=0L,nfilts-1 do begin
  
     mtflt = where(strmatch(hamspec[gdflt].block, filts[ss]))

     print, 'hamspec_mkflats: Stacking filter ', filts[ss]

     idx = gdflt[mtflt]

     ;; Median Combine
     if n_elements(mtflt) GT 1 then begin
        ;; Filenames
        ovfil = ['']
        for q=0L,n_elements(mtflt)-1 do begin
           ovfil = [ovfil, $
                    hamspec_getfil('ov_fil', $
                                   FRAME=hamspec[idx[q]].frame, /name)]
        endfor
        ovfil = ovfil[1:*]
        ;; Median or average?
        if keyword_set(AVG) then fcomb = 3 else fcomb = 2
        ;; Combine
        xcombine, ovfil, img_flat, head, $
                  FCOMB=fcomb, SCALE='MED', GAIN=hamspec[idx[0]].gain, $
                  RN=hamspec[idx[0]].readno
     endif else $
        img_flat = xmrdfits(hamspec_getfil('ov_fil', $
                                           FRAME=hamspec[idx[0]].frame,/name), 0, head)
  
     ;; Inverse variance
     ivar = 1./((img_flat>1.)*hamspec[idx[0]].gain + $
                hamspec[idx[0]].readno^2) $
            * (img_flat GT 1 )
     ;; Save
     if ss EQ 0 then begin
        fin_flat = img_flat 
        fin_ivar = ivar
     endif else begin
        fin_flat = img_flat + fin_flat
        fin_ivar = 1./(1./fin_ivar + 1./ivar)
     endelse
  endfor

  
  
  ;; Output
  mwrfits, fin_flat, outfil, head, /create, /silent
  mwrfits, fin_ivar, outfil, /silent
  spawn, 'gzip -f '+outfil
  print, 'hamspec_mktflat: Flat created ', outfil+'.gz.'
  
  ;; DEL OV
  if not keyword_set(SVOV) then hamspec_delov, hamspec, gdflt, /silent

  print, 'hamspec_mktflat: All done!'
  return
end
