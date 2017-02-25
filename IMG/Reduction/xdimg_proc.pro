;+ 
; NAME:
; xdimg_proc   
;    Version 1.1
;
; PURPOSE:
;    Processes a set of images  (OV, TRIM, FLAT)
;
; CALLING SEQUENCE:
;   
;   xdimg_proc, struct, img, flat, OUTPTH=, /IR, /INTER, /SVOV, /SVCODIV, 
;               /SVDRKSUB, /GAIN, /SILENT, MSK=, BPM=, FILEOUT=
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
;  OUTPTH  =  Output directory (default = 'Final/')
;  /INTER  =  Interactive OV fitting
;  /SVOV  =  Delete ov files when through
;  /SVCODIV = Delete coadd divided files when through
;  /SVDRKSUB = Delete dark subtracted files when through
;  MSK     =  Default name for mask file
;  /GAIN =  Multiply by the gain
;  BPM     =  Name of bad pixel mask file to apply to image. (including directory)
;
; OPTIONAL OUTPUTS:
;  FILEOUT = Name of output file which has the names of the reduced images.
;            Useful for call to stacking procedure.
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
;   14-June-2007 Updated by LKP to allow for IR images, like NIRC2.
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro xdimg_proc, struct, img, flatnm, IR=ir, OUTPTH=outpth, INTER=inter, SVOV=delov, $
                SVCODIV=delcodiv, SVDRKSUB=deldrksub, MSK=msk, GAIN=gain, $
                SILENT=silent, BPM=bpm, FILEOUT=fileout

;
  if  N_params() LT 3  then begin 
      print, 'Syntax - ' +$
        'xdimg_proc, struct, img, flat, OUTPTH=, MSK=, /IR, /INTER, /SVOV, /GAIN'
      print, '       /SVCODIV, /SVDRKSUB, BPM=bpm, FILEOUT=fileout, /SILENT [v1.1]'
      return
  endif 

  
;  Optional Keywords
  if not keyword_set(OUTPTH) then outpth = 'Final/'
  if keyword_set(GAIN) AND not keyword_set( SILENT ) then $
    print, 'xdimg_proc: Multiplying each image by the gain!'

; Dealing with IR data
  if keyword_set(IR) then irdata=1 else irdata=0

  if irdata EQ 0 and strmid(struct[img[0]].ccd,0,5) eq 'NIRC2' then begin
     print, 'When dealing with IR data you must use the /IR flag.  Aborting!'
     return
  endif
    
  ;  Preprocessing
  ;  Dark subtract if we're dealing with infrared CCDs.
  if irdata EQ 1 then begin
     preprocess = where(struct[img].flg_drksub eq 0, npreprocess)
     if npreprocess gt 0 then xdimg_darksub, struct, img[preprocess]
  endif else begin
     ; Overscan subtract if we're dealing with optical CCDs.
     ovflt = where(struct[img].flg_ov EQ 0, nov)
     if nov NE 0 then xdimg_over, struct, img[ovflt], ORDR=5, INTER=inter
  endelse

;  Find all the filters involved
    x_filters, struct[img].filter, filt, nfilt
  
;  Make an array which stores the names of the files processed.
    writeout=strarr(n_elements(img))
    counter=0


; Loop on separate filters
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
         partial_chip=0

         ; Check sizes
         if irdata EQ 1 then data = mrdfits(strtrim(struct[img[wfilt[j]]].img_drksub, 2), 0, himg, /silent, /fscale) else $
            data =  mrdfits(strtrim(struct[img[wfilt[j]]].img_ov, 2), 0, himg, /silent, /fscale)

         szimg = size(data, /dimensions)
         if szimg[0] NE szflat[0] OR szimg[1] NE szflat[1] then begin
            window=sxpar(himg, 'WINDOW')
            ; Must use correct portion of flat if img is partial read-out.
            if window eq '0,0,550,2048,1000' then begin
               partial_chip=1
            endif else begin
               print, 'xdimg_proc: Inconsistent image sizes! Aborting...'
               print, 'NOTE: Image sizes cannot be reconciled by assuming LRISR science image was accidentally done in partial read-out mode.'
               return
            endelse
         endif
         
         ; Non-linearity
         if not keyword_set( NONONLIN) then $
            data = xdimg_nonlinear(data, struct[img[wfilt[j]]].ccd)
         
         ; Flatten
         if partial_chip eq 1 then begin
            fin = data / flat[*,498:1497]
         endif else begin
            fin = data / flat
         endelse
         
         if irdata eq 1 then $
            print, 'xdimg_proc: Flattening ', struct[img[wfilt[j]]].img_drksub, ' with ', flatfil else $
               print, 'xdimg_proc: Flattening ', struct[img[wfilt[j]]].img_ov, ' with ', flatfil
         
         
         ; Apply the Gain and reset Saturation
         if keyword_set( GAIN ) then begin
            fin = fin * struct[img[wfilt[j]]].gain
            struct[img[wfilt[j]]].satur = struct[img[wfilt[j]]].satur * $
                                          struct[img[wfilt[j]]].gain
         endif
         
         ; Modify header
         sxaddpar, himg, 'FLAT', flatfil
         
         if keyword_set( BPM ) then begin
            ; check for existing BPM
            a = findfile(strtrim(BPM, 2), count=count)
            if count EQ 0 then begin
               print, 'Requested bad pixel mask does not exit.  Proceeding.'
            endif else begin
               badpix = mrdfits(strtrim(BPM,2), /silent)
               
               ; Set bad pixels to zero.
               fin = fin * badpix
               
               ; Modify header
               sxaddpar, himg, 'BPM', strtrim(BPM,2)
               
            endelse
         endif
         
         ; Not sure what X was doing here.
         if keyword_set( MSK ) then begin
            ; check for existing BPM
            a = sxpar(himg, 'BPM', count=acnt)
            ; Add parameter
            if acnt EQ 0 then sxaddpar, himg, 'BPM', msk
         endif
         
         
         ; Output
         outfil = strjoin([outpth, 'f_', strtrim(struct[img[wfilt[j]]].img_root, 2)])
         mwrfits, fin, outfil, himg, /create

         ; Document which files have been written out.
         writeout[counter]='f_'+strtrim(struct[img[wfilt[j]]].img_root, 2)
         counter=counter+1

         ; Update flg_final
         ; NOTE: not sure what X's intention was here, with different values flags.
         if keyword_set( GAIN ) then struct[img[wfilt[j]]].flg_final = 4 $
         else struct[img[wfilt[j]]].flg_final = 7
         
         struct[img[wfilt[j]]].img_final = outfil
      endfor
   endfor
  
  ; Delete files
  if irdata EQ 1 then begin
     if not keyword_set( SVDRKSUB ) then xdimg_deldrksub, struct, img
  endif else begin
     if not keyword_set( SVOV ) then xdimg_delov, struct, img
  endelse

  
  ; When I tried to implement this in the loop the unit number kept getting messed up.
  ; So I've saved an array with the info and write the file now.
  if keyword_set(FILEOUT) then begin
     openw, unit, outpth+fileout, /get_lun
     print, 'Opening file for writing: ' + outpth+fileout
     wait, 1

     for l=0, n_elements(img)-1 do begin
        printf, unit, writeout[l]
     endfor

     close, unit
     free_lun, unit
  endif


   ; Resave the updated structure so that you can pick up where you left off easily.
   mwrfits, struct, 'struct.fits', /create

  ; All done
  if not keyword_set( SILENT ) then $
    print, 'xdimg_proc: All done'


end
