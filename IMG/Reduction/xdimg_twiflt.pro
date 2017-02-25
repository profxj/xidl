;+ 
; NAME:
; xdimg_twiflt   
;  Version 1.1
;
; PURPOSE:
;    Creates twilight flats given the image list structure (assuming
;     twilight flats were taken!)
;
; CALLING SEQUENCE:
;   
;   xdimg_twiflt, struct, /SVOV, OUTROOT=
;
; INPUTS:
;   struct -- 'distruct' defining the images of interest
;
; RETURNS:
;
; OUTPUTS:
;   twiflats - fits files; 1 per filter
;
; OPTIONAL KEYWORDS:
;   SVOV - Save ov files
;   OUTROOT - Root name of Twilight flats (default is 'Flats/TwiFlat')
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   xdimg_twiflt, nght1_strct
;
;
; PROCEDURES/FUNCTIONS CALLED:
;  XDIMG_OVER
;  XCOMBINE
;  MWRFITS
;  XDIMG_DELOV
;  X_FILTERS
;
; REVISION HISTORY:
;   26-July-2001 Written by JXP
;   24-Apr-2002 Added nonlinearity correction
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro xdimg_twiflt, struct, SVOV=svov, OUTROOT=outroot

;
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'xdimg_twiflt, struct, /SVOV, OUTROOT= (v1.1)'
      return
  endif 
  
;  Optional Keywords

  if not keyword_set( OUTROOT ) then outroot = 'Flats/TwiFlat'
  outroot = strtrim(outroot,2)
  
;  Find the Twilight Flats

  tflts = where(struct.type EQ 'TWI' AND struct.flg_anly NE 0, ntflt)

;  Overscan

  if ntflt EQ 0 then begin
      print, 'xdimg_twiflt:  No TWI images!'
      return
  endif

  ovflt = where(struct[tflts].flg_ov EQ 0, nov)
  if nov NE 0 then xdimg_over, struct, tflts[ovflt], ORDR=4

;  Find all the filters involved
  
  x_filters, struct[tflts].filter, filt, nfilt
  
;  Loop on separate filters
  
  for q=0,nfilt-1 do begin
      wfilt = where(strtrim(struct[tflts].filter,2) EQ filt[q], dumi)

      outfil = strjoin([outroot, strtrim(filt[q],2), '.fits'])
      outfilN = strjoin([outroot, 'N', strtrim(filt[q],2), '.fits'])

      ; Status
      print, 'Combining images: '
      for i=0,dumi-1 do print, struct[tflts[wfilt[i]]].img_ov
      print, '             into ', outfil

      ; Combine Images
      xcombine, struct[tflts[wfilt]].img_ov, comb, head, FCOMB=2, $
        GAIN=struct[tflts[wfilt[0]]].gain, RN=struct[tflts[wfilt[0]]].readno, $
        SCALE=struct[tflts[wfilt]].med_ov

      ; Non-linearity
      if not keyword_set( NONONLIN) then $
        comb = xdimg_nonlinear(comb, struct[tflts[0]].ccd)

      ; Normalize
      ssec = struct[tflts[wfilt[0]]].statsec  ; Stat sec
      med = median(comb[ssec[0]:ssec[1],ssec[2]:ssec[3]])
      combN = comb / med

      ; Header
      headN = head
      fxaddpar, headN, 'MEDFLT', med

      ; Output
      if keyword_set( ALLOUT ) then mwrfits, comb, outfil, head, /create
      mwrfits, combN, outfilN, headN, /create

  endfor

; Delete images  
  if not keyword_set( SVOV ) then xdimg_delov, struct, tflts

; Resave the updated structure so that you can pick up where you left off easily.
  mwrfits, struct, 'struct.fits', /create

  print, 'xdimg_twiflt: All done with Twilight Flats!'

end
