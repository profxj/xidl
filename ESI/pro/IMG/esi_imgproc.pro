;+ 
; NAME:
; esi_imgproc   
;    Version 1.0
;
; PURPOSE:
;    Processes a set of images  (OV, TRIM, FLAT)
;
; CALLING SEQUENCE:
;   
;   esi_imgproc, esi, img, flat, OUTPTH=, /INTER, /DELOV, /NOGAIN
;
; INPUTS:
;   esi -- dimg_strct defining the images of interest
;   img    -- intarr of images to process
;   flatnm -- Flat root to use (e.g. 'Flats/SkyFltN')
;
; RETURNS:
;
; OUTPUTS:
;   Processed image (e.g. Final/f_ccd001.fits)
;
; OPTIONAL KEYWORDS:
;  OUTPTH =  Output directory (default = 'Final/')
;  INTER  =  Interactive OV fitting
;  DELOV  =  Delete ov files when through
;  MSK    =  Default name for mask file
;  NOGAIN =  Do not apply the gain
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   esi_imgproc, dimg, STDS, 'Flats/SkyFltN'
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   30-Jul-2002 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro esi_imgproc, esi, obj_id, ALLOBJ=allobj, FLATNM=flatnm, $
                 OUTPTH=outpth, INTER=inter, DELOV=delov, $
                MSK=msk, NOGAIN=nogain, SILENT=silent, FLATFLG=flatflg

;
  if  N_params() LT 1  then begin 
      print, 'Syntax - ' +$
        'esi_imgproc, esi, obj_id, FLATNM=, OUTPTH=, MSK=, /DELOV, /NOGAIN'
      print, '       /SILENT, /ALLOBJ [v1.1]'
      return
  endif 

  
;  Optional Keywords

  if not keyword_set(OUTPTH) then outpth = 'Final/'
  if not keyword_set(NOGAIN) AND not keyword_set( SILENT ) then $
    print, 'esi_imgproc: Multiplying each image by the gain!'
  if not keyword_set( FLATNM ) then begin
      if not keyword_set( FLATFLG ) then flatflg = 0L
      case flatflg of
          0L: flatnm = 'Flats/FlatIMG_T'
          1L: flatnm = 'Flats/FlatIMG_S'
          2L: flatnm = 'Flats/FlatIMG_D'
          else: stop
      endcase
  endif
  
; IMG

  if keyword_set( ALLOBJ ) then begin
      img = where(esi.flg_anly NE 0 AND esi.mode EQ 0 AND $
                  (esi.type EQ 'STD' OR esi.type EQ 'OBJ'), nimg)
  endif else begin
      if not keyword_set( OBJ_ID ) then obj_id = 0L
      img = where(esi.flg_anly NE 0 AND esi.mode EQ 0 AND $
                  esi.obj_id EQ obj_id, nimg)
  endelse

;  OV

  bias = where(esi[img].flg_ov EQ 0, nbias) 
  if nbias NE 0 then esi_subbias, esi, img[bias]
  

;  Find all the filters involved
  
  x_filters, esi[img].imfilt, filt, nfilt
  
;;;;;;;;;;;;;;;;;;;;
;  Loop on separate filters
  for q=0L,nfilt-1 do begin

      ; Read Flat file
      flatfil = strjoin([strtrim(flatnm,2), filt[q], '.fits'])
      flat = xmrdfits(flatfil, /silent, /fscale)
      szflat = size(flat, /dimensions)

      wfilt = where(strtrim(esi[img].imfilt,2) EQ filt[q], nimg)

      if not keyword_set( SILENT ) then $
        print, 'esi_imgproc: Processing filter '+filt[q]

      ;;;;;;;;
      ; Loop on Images
      for j=0L,nimg-1 do begin

          ;; Check sizes
          data = xmrdfits(esi[img[wfilt[j]]].img_ov, 0, himg, /silent, /fscale)
          szimg = size(data, /dimensions)
          if szimg[0] NE szflat[0] OR szimg[1] NE szflat[1] then begin
              print, 'esi_imgproc: Inconsistent image sizes! Aborting...'
              return
          endif
          
          ;; Non-linearity
;          if not keyword_set( NONONLIN) then $
;            data = esi_imgnonlinear(data, esi[img[wfilt[j]]].ccd)

          ;; Flatten
          fin = data / flat
          print, 'esi_imgproc: Flattening ', esi[img[wfilt[j]]].img_ov, ' with ', $
            flatfil

          ;; Apply the Gain and reset Saturation

          if not keyword_set( NOGAIN ) then begin
              fin = fin * esi[img[wfilt[j]]].gain
;              esi[img[wfilt[j]]].satur = esi[img[wfilt[j]]].satur * $
;                esi[img[wfilt[j]]].gain
          endif

          ;; Modify header
          sxaddpar, himg, 'FLAT', flatfil
          if keyword_set( MSK) then begin
              ; check for existing BPM
              a = sxpar(himg, 'BPM', count=acnt)
              ; Add parameter
              if acnt EQ 0 then sxaddpar, himg, 'BPM', msk
          endif
          
          

          ;; Output
          outfil = strjoin([outpth, 'f_', esi[img[wfilt[j]]].img_root])
          mwrfits, fin, outfil, himg, /create
          spawn, 'gzip -f '+outfil

          ;; Update flg_final
          if keyword_set( NOGAIN ) then esi[img[wfilt[j]]].flg_final = 4 $
          else esi[img[wfilt[j]]].flg_final = 7

          esi[img[wfilt[j]]].img_final = outfil
      endfor
  endfor

  ; Delete OV files
  if keyword_set( DELOV ) then esi_imgdelov, esi, img
  ; All done
  if not keyword_set( SILENT ) then $
    print, 'esi_imgproc: All done'
end
