;+ 
; NAME:
; esi_imgskyflat   
;   Version 1.1
;
; PURPOSE:
;    Creates Super-sky flat given the image list esiure.  
;    The routine first OV subtracts each image as necessary.
;    It then calls esi_imgskymask create masks for each of the images.
;    It then calls xcombine which: (1) scales each image 
;    by the median calculated for the OV image (2) does a median
;    combine with 3sigma/2.5sigma (low/high) clipping.
;
; CALLING SEQUENCE:
;   
;   esi_imgskyflat, esi, MMEM=, /NOMSK, /INTER
;
; INPUTS:
;   esi -- dimg_strct defining the images of interest
;
; RETURNS:
;
; OUTPUTS:
;   skyflat - fits files; 1 per filter
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
;   esi_imgskyflat, nght1_strct
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

pro esi_imgskyflat, esi, MMEM=mmem, NOMSK=nomsk, INTER=inter,  $
                  OUTROOT=outroot, ORDR=ordr, XSIZE=xsize, YSIZE=ysize

;
  if  N_params() LT 1  then begin 
      print, 'Syntax - ' +$
        'esi_imgskyflat, esi, MMEM=mmem, /NOMSK, /INTER, OUTROOT=, ORDR= (v1.1)' 
      return
  endif 

  
;  Optional Keywords

  if not keyword_set( ORDR ) then ordr = 5
  if not keyword_set( OUTROOT ) then outroot = 'Flats/FlatIMG_S'
  outroot = strtrim(outroot,2)
  if not keyword_set(MMEM) then mmem = 200.
  
;  Find the Object Images

  skyflat = where(esi.type EQ 'OBJ' AND esi.flg_anly NE 0 $
                 AND esi.mode EQ 0, nskyflat)


;  Overscan

  if not keyword_set( REDOOV ) then bias = where(esi[skyflat].flg_ov EQ 0, nbias) $
  else begin
      nbias = n_elements(skyflat)
      bias = lindgen(nbias)
  endelse
  if nbias NE 0 then esi_subbias, esi, skyflat[bias]

;  Create Mask images

;  if not keyword_set( NOMSK ) then $
;    esi_imgskymsk, esi, skyflat, XSIZE=xsize, YSIZE=ysize

;  Find all the filters involved
  
  x_filters, esi[skyflat].imfilt, filt, nfilt
  
;  Loop on separate filters
  
  for q=0L,nfilt-1 do begin
      wfilt = where(esi[skyflat].imfilt EQ filt[q], dumi)

      outfil = strjoin([outroot, strtrim(filt[q],2), '.fits'])

      ;; Status
      print, 'esi_imgskyflat: Combining images: '
      for i=0L,dumi-1 do print, esi[skyflat[wfilt[i]]].img_ov
      print, '             into ', outfil

      ;; Combine Images
      xcombine, esi[skyflat[wfilt]].img_ov, comb, head, FCOMB=4, MMEM=mmem, $
        SCALE='MED', GAIN=esi[skyflat[wfilt[0]]].gain, $
        RN=esi[skyflat[wfilt[0]]].readno, SIGLO=3.0, SIGHI=2.5
;        MASKS=esi[skyflat[wfilt]].img_skymsk, SIGLO=3.0, SIGHI=2.5
        
      ;; Deal with non-zero
      a = where(comb EQ 0., na)
      if na NE 0 then comb[a] = 1.

      ; Non-linearity
;      if not keyword_set( NONONLIN) then $
;        comb = esi_imgnonlinear(comb, esi[skyflat[0]].ccd)

      ;; Normalize
      ssec = [140,580L, 90,1060L] 
      med = median(comb[ssec[0]:ssec[1],ssec[2]:ssec[3]])
      comb = comb / med


      ;; Header
      fxaddpar, head, 'MEDFLT', med

      ;; Output
      print, 'esi_imgdflat: Creating -- ', outfil
      mwrfits, comb, outfil, head, /create
  endfor
  
  print, 'esi_imgskyflat: All done with SuperSky Flat!'

end
