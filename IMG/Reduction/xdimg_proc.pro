;+ 
; NAME:
; xdimg_proc   
;    Version 1.0
;
; PURPOSE:
;    Processes a set of images  (OV, TRIM, FLAT)
;
; CALLING SEQUENCE:
;   
;   xdimg_proc, struct, img, flat, OUTPTH=, /INTER, /DELOV, /NOGAIN
;
; INPUTS:
;   struct -- dimg_strct defining the images of interest
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
;   xdimg_proc, dimg, STDS, 'Flats/SkyFltN'
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   03-Aug-2001 Written by JXP
;   03-Jan-2002 Dealt with Gain
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro xdimg_proc, struct, img, flatnm, OUTPTH=outpth, INTER=inter, DELOV=delov, $
                MSK=msk, NOGAIN=nogain, SILENT=silent

;
  if  N_params() LT 3  then begin 
      print, 'Syntax - ' +$
        'xdimg_proc, struct, img, flat, OUTPTH=, MSK=, /INTER, /DELOV, /NOGAIN'
      print, '       /SILENT [v1.1]'
      return
  endif 

  
;  Optional Keywords

  if not keyword_set(OUTPTH) then outpth = 'Final/'
  if not keyword_set(NOGAIN) AND not keyword_set( SILENT ) then $
    print, 'xdimg_proc: Multiplying each image by the gain!'
  

;  OV

  ovflt = where(struct[img].flg_ov EQ 0, nov)
  if nov NE 0 then xdimg_over, struct, img[ovflt], ORDR=5, INTER=inter
  

;  Find all the filters involved
  
  x_filters, struct[img].filter, filt, nfilt
  
;;;;;;;;;;;;;;;;;;;;
;  Loop on separate filters
  for q=0L,nfilt-1 do begin

      ; Read Flat file
      flatfil = strjoin([flatnm, filt[q], '.fits'])
      flat = mrdfits(flatfil, /silent, /fscale)
      szflat = size(flat, /dimensions)

      wfilt = where(strtrim(struct[img].filter,2) EQ filt[q], nimg)

      if not keyword_set( SILENT ) then $
        print, 'xdimg_proc: Processing filter '+filt[q]

      ;;;;;;;;
      ; Loop on Images
      for j=0L,nimg-1 do begin

          ; Check sizes
          data = mrdfits(struct[img[wfilt[j]]].img_ov, 0, himg, /silent, /fscale)
          szimg = size(data, /dimensions)
          if szimg[0] NE szflat[0] OR szimg[1] NE szflat[1] then begin
              print, 'xdimg_proc: Inconsistent image sizes! Aborting...'
              return
          endif
          
          ; Non-linearity
          if not keyword_set( NONONLIN) then $
            data = xdimg_nonlinear(data, struct[img[wfilt[j]]].ccd)

          ; Flatten
          fin = data / flat
          print, 'xdimg_proc: Flattening ', struct[img[wfilt[j]]].img_ov, ' with ', $
            flatfil

          ; Apply the Gain and reset Saturation

          if not keyword_set( NOGAIN ) then begin
              fin = fin * struct[img[wfilt[j]]].gain
              struct[img[wfilt[j]]].satur = struct[img[wfilt[j]]].satur * $
                struct[img[wfilt[j]]].gain
          endif

          ; Modify header
          sxaddpar, himg, 'FLAT', flatfil
          if keyword_set( MSK) then begin
              ; check for existing BPM
              a = sxpar(himg, 'BPM', count=acnt)
              ; Add parameter
              if acnt EQ 0 then sxaddpar, himg, 'BPM', msk
          endif
          
          

          ; Output
          outfil = strjoin([outpth, 'f_', struct[img[wfilt[j]]].img_root])
          mwrfits, fin, outfil, himg, /create

          ; Update flg_final

          if keyword_set( NOGAIN ) then struct[img[wfilt[j]]].flg_final = 4 $
          else struct[img[wfilt[j]]].flg_final = 7

          struct[img[wfilt[j]]].img_final = outfil
      endfor
  endfor

  ; Delete OV files
  if keyword_set( DELOV ) then xdimg_delov, struct, img
  ; All done
  if not keyword_set( SILENT ) then $
    print, 'xdimg_proc: All done'
end
