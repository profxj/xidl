;+ 
; NAME:
; xdimg_dflat   
;  Version 1.1
;
; PURPOSE:
;    Creates dome flats given the image list structure
;
; CALLING SEQUENCE:
;   
;   xdimg_dflat, struct, /SVOV, OUTROOT=
;
; INPUTS:
;   struct -- dimg_strct defining the images of interest
;
; RETURNS:
;
; OUTPUTS:
;   flats - fits files in the dir Flats; 1 per filter
;
; OPTIONAL KEYWORDS:
;   SVOV - save ov files
;   OUTROOT - Root name of Dome flats (default is 'Flats/DFlat')
;   ALLOUT - Output unnormalized frame too
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   xdimg_dflat, nght1_strct
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
;   18-July-2001 Written by JXP
;   24-Apr-2002 Allow for non-linearity
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro xdimg_dflat, struct, SVOV=svov, OUTROOT=outroot, ALLOUT=allout

;
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'xdimg_dflat, struct, /SVOV, OUTROOT=, ALLOUT= (v1.1)'
      return
  endif 
  
;  Optional Keywords

  if not keyword_set( OUTROOT ) then outroot = 'Flats/DFlat'
  outroot = strtrim(outroot,2)
  
;  Find the Dome Flats

  dflts = where(struct.type EQ 'DFT' AND struct.flg_anly NE 0, ndflt)

;  Overscan

  ovflt = where(struct[dflts].flg_ov EQ 0, nov)
  if nov NE 0 then xdimg_over, struct, dflts[ovflt], ORDR=4

;  Find all the filters involved
  
  x_filters, struct[dflts].filter, filt, nfilt
  
;  Loop on separate filters
  
  for q=0,nfilt-1 do begin
      wfilt = where(struct[dflts].filter EQ filt[q], dumi)

      outfil = strjoin([outroot, strtrim(filt[q],2), '.fits'])
      outfilN = strjoin([outroot, 'N', strtrim(filt[q],2), '.fits'])

      ; Status
      print, 'Combining images: '
      for i=0,dumi-1 do print, struct[dflts[wfilt[i]]].img_ov
      print, '             into ', outfil

      ; Combine Images
      xcombine, struct[dflts[wfilt]].img_ov, comb, head, $
        SCALE=struct[dflts[wfilt]].med_ov

      ; Non-linearity
      if not keyword_set( NONONLIN) then $
        comb = xdimg_nonlinear(comb, struct[dflts[0]].ccd)

      ; Normalize
      ssec = struct[dflts[wfilt[0]]].statsec  ; Stat sec
      med = median(comb[ssec[0]:ssec[1],ssec[2]:ssec[3]])
      combN = comb / med

      ; Header
      headN = head
      fxaddpar, headN, 'MEDFLT', med

      ; Output
      if keyword_set( ALLOUT ) then mwrfits, comb, outfil, head, /create
      mwrfits, combN, outfilN, headN, /create

  endfor
  
;  Delete images

  if not keyword_set( SVOV ) then xdimg_delov, struct, dflts

  print, 'All done with Dome Flats!'

end
