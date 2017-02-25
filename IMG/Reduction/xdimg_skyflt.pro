;+ 
; NAME:
; xdimg_skyflt   
;   Version 1.1
;
; PURPOSE:
;    Creates Super-sky flat given the image list structure.  
;    The routine first OV subtracts (or dark subtracts, for IR images) each image as necessary.
;    It then calls xdimg_skymask to create masks for each of the images (unless the /NOMSK keyword is set).
;    It then calls xcombine which: (1) scales each image 
;    by the median calculated for the OV (or dark subtracted) image (2) does a median
;    combine with 3sigma/2.5sigma (low/high) clipping.
;
; CALLING SEQUENCE:
;   
;   xdimg_skyflt, struct, MMEM=, /NOMSK, /INTER, /IR, CCD='NIRC2W'
;
; INPUTS:
;   struct -- dimg_strct defining the images of interest
;
; RETURNS:
;
; OUTPUTS:
;   skyflt - fits files; 1 per filter (and per CCD, if CCD option is given)
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
;  XDIMG_DELDRKSUB
;  X_FILTERS
;  XDIMG_SKYMASK
;
; REVISION HISTORY:
;   26-July-2001 Written by JXP
;   24-Apr-2002 Added nonlinearity correction
;   25-Jun-2007  LKP added ability to use with IR, NIRC2 images
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro xdimg_skyflt, struct, IR=ir, CCD=ccd, MMEM=mmem, NOMSK=nomsk, INTER=inter,  $
                  OUTROOT=outroot, ORDR=ordr, XSIZE=xsize, YSIZE=ysize, SVOV=svov, SVDRKSUB=svdrksub

  if  N_params() LT 1  then begin 
      print, 'Syntax - ' +$
        'xdimg_skyflt, struct, /IR, CCD=ccd, MMEM=mmem, /NOMSK, /INTER, OUTROOT=, ORDR=, /SVOV, SVDRKSUB(v1.1)' 
      return
  endif 

  
;  Optional Keywords
  if not keyword_set( ORDR ) then ordr = 5
  if not keyword_set( OUTROOT ) then begin
     if keyword_set( CCD ) then begin
        CASE ccd of
           'NIRC2W': begin
              outroot = 'Flats/SkyFlt_w_'
           end
           'NIRC2N': begin
              outroot = 'Flats/SkyFlt_n_'
           end
        end
     endif else begin
        outroot = 'Flats/SkyFlt'
     endelse
     outroot = strtrim(outroot,2)
  endif
  if not keyword_set(MMEM) then mmem = 200.


; Dealing with IR data
  if keyword_set(IR) then irdata=1 else irdata=0
  if irdata eq 1 then NOMSK=1

;  Find the Object Images
  if keyword_set( CCD ) then skyflt = where(struct.type EQ 'OBJ' AND struct.ccd EQ ccd AND struct.flg_anly NE 0, nskyflt) $
                        else skyflt = where(struct.type EQ 'OBJ' AND struct.flg_anly NE 0, nskyflt)

  if irdata EQ 0 and strmid(struct[skyflt[0]].ccd,0,5) eq 'NIRC2' then begin
     print, 'When dealing with IR data you must use the /IR flag.  Aborting!'
     return
  endif

; Preprocessing
  ; Dark subtract if we're dealing with infrared CCDs.
  if irdata EQ 1 then begin
     preprocess = where(struct[skyflt].flg_drksub eq 0, npreprocess)
     if npreprocess gt 0 then xdimg_darksub, struct, skyflt[preprocess]
  endif else begin
     ;  Overscan subtract if we're dealing with optical CCDs.
     ovflt = where(struct[skyflt].flg_ov EQ 0, nov)
     if nov NE 0 then xdimg_over, struct, skyflt[ovflt], ORDR=ordr, INTER=inter
  endelse

;  Create Mask images
  if not keyword_set( NOMSK ) then $
    xdimg_skymsk, struct, skyflt, XSIZE=xsize, YSIZE=ysize

;  Find all the filters involved
    x_filters, struct[skyflt].filter, filt, nfilt

  
;  Loop on separate filters  
  for q=0L,nfilt-1 do begin
      wfilt = where(strtrim(struct[skyflt].filter,2) EQ filt[q], dumi)

      outfil = strjoin([outroot, strtrim(filt[q],2), '.fits'])
      outfilN = strjoin([outroot, 'N', strtrim(filt[q],2), '.fits'])

      ; Status
      print, 'Combining images: '
      if irdata EQ 1 then begin
         for i=0, dumi-1 do print, struct[skyflt[wfilt[i]]].img_drksub
      endif else begin
         for i=0L,dumi-1 do print, struct[skyflt[wfilt[i]]].img_ov
      endelse
      print, '             into ', outfil

      ; Combine Images
      ; Use dark subtracted images if dealing with IR data
      ; use overscan subtracted images if dealing with optical data.
      for hh=0, dumi-1 do begin
         if irdata EQ 1 then $
            struct[skyflt[wfilt[hh]]].img_drksub=strtrim(struct[skyflt[wfilt[hh]]].img_drksub,2) else $
               struct[skyflt[wfilt[hh]]].img_ov = strtrim(struct[skyflt[wfilt[hh]]].img_ov,2)
      endfor

      if irdata EQ 1 then $
         xcombine, struct[skyflt[wfilt]].img_drksub, comb, head, FCOMB=4, MMEM=mmem, $
                   SCALE=struct[skyflt[wfilt]].med_drksub, GAIN=struct[skyflt[wfilt[0]]].gain, $
                   RN=struct[skyflt[wfilt[0]]].readno, $
                   SIGLO=3.0, SIGHI=2.5 $
                   else xcombine, struct[skyflt[wfilt]].img_ov, comb, head, FCOMB=4, MMEM=mmem, $
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

; Delete images
  if irdata EQ 1 then begin
     if not keyword_set( SVDRKSUB ) then xdimg_deldrksub, struct, skyflt
  endif else begin
     if not keyword_set( SVOV ) then xdimg_delov, struct, skyflt
  endelse
  
  
; Resave the updated structure so that you can pick up where you left off easily.
  mwrfits, struct, 'struct.fits', /create
  print, 'All done with sky flats!'
  
end
