;+ 
; NAME:
; xdimg_skyflt   
;   Version 1.1
;
; PURPOSE:
;    Creates Super-sky flat given the image list structure.  
;    The routine first OV subtracts each image as necessary.
;    It then calls xdimg_skymask create masks for each of the images.
;    It then calls xcombine which: (1) scales each image 
;    by the median calculated for the OV image (2) does a median
;    combine with 3sigma/2.5sigma (low/high) clipping.
;
; CALLING SEQUENCE:
;   
;   xdimg_skyflt, struct, MMEM=, /NOMSK, /INTER
;
; INPUTS:
;   struct -- dimg_strct defining the images of interest
;
; RETURNS:
;
; OUTPUTS:
;   skyflt - fits files; 1 per filter
;
; OPTIONAL KEYWORDS:
;  MMEM - Max memory to use with this routine (default = 200M)
;  NOMSK - Don't create new masks
;  OUTROOT - Root name of Sky flats (default is 'Flats/SkyFlt')
;  INTER - Interactive OV subtraction
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   xdimg_skyflt, nght1_strct
;
;
; PROCEDURES/FUNCTIONS CALLED:
;  XDIMG_OVER
;  XCOMBINE
;  MWRFITS
;  XDIMG_DELOV
;  X_FILTERS
;  XDIMG_SKYMASK
;
; REVISION HISTORY:
;   26-July-2001 Written by JXP
;   24-Apr-2002 Added nonlinearity correction
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro xdimg_skyflt, struct, MMEM=mmem, NOMSK=nomsk, INTER=inter,  $
                  OUTROOT=outroot, ORDR=ordr, XSIZE=xsize, YSIZE=ysize

;
  if  N_params() LT 1  then begin 
      print, 'Syntax - ' +$
        'xdimg_skyflt, struct, MMEM=mmem, /NOMSK, /INTER, OUTROOT=, ORDR= (v1.1)' 
      return
  endif 

  
;  Optional Keywords

  if not keyword_set( ORDR ) then ordr = 5
  if not keyword_set( OUTROOT ) then outroot = 'Flats/SkyFlt'
  outroot = strtrim(outroot,2)
  if not keyword_set(MMEM) then mmem = 200.
  
;  Find the Object Images

  skyflt = where(struct.type EQ 'OBJ' AND struct.flg_anly NE 0, nskyflt)


;  Overscan

  ovflt = where(struct[skyflt].flg_ov EQ 0, nov)
  if nov NE 0 then xdimg_over, struct, skyflt[ovflt], ORDR=ordr, INTER=inter

;  Create Mask images

  if not keyword_set( NOMSK ) then $
    xdimg_skymsk, struct, skyflt, XSIZE=xsize, YSIZE=ysize

;  Find all the filters involved
  
  x_filters, struct[skyflt].filter, filt, nfilt
  
;  Loop on separate filters
  
  for q=0L,nfilt-1 do begin
      wfilt = where(struct[skyflt].filter EQ filt[q], dumi)

      outfil = strjoin([outroot, strtrim(filt[q],2), '.fits'])
      outfilN = strjoin([outroot, 'N', strtrim(filt[q],2), '.fits'])

      ; Status
      print, 'Combining images: '
      for i=0L,dumi-1 do print, struct[skyflt[wfilt[i]]].img_ov
      print, '             into ', outfil

      ; Combine Images
      xcombine, struct[skyflt[wfilt]].img_ov, comb, head, FCOMB=4, MMEM=mmem, $
        SCALE=struct[skyflt[wfilt]].med_ov, GAIN=struct[skyflt[wfilt[0]]].gain, $
        RN=struct[skyflt[wfilt[0]]].readno, $
        MASKS=struct[skyflt[wfilt]].img_skymsk, SIGLO=3.0, SIGHI=2.5
        
      ; Deal with non-zero
      a = where(comb EQ 0., na)
      if na NE 0 then comb[a] = 1.

      ; Non-linearity
      if not keyword_set( NONONLIN) then $
        comb = xdimg_nonlinear(comb, struct[skyflt[0]].ccd)

      ; Normalize
      ssec = struct[skyflt[wfilt[0]]].statsec  ; Stat sec
      med = median(comb[ssec[0]:ssec[1],ssec[2]:ssec[3]])
      combN = comb / med


      ; Header
      headN = head
      fxaddpar, headN, 'MEDFLT', med

      ; Output
      if keyword_set( ALLOUT ) then mwrfits, comb, outfil, head, /create
      mwrfits, combN, outfilN, head, /create

  endfor
  
end
