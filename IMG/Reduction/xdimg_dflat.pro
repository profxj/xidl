;+ 
; NAME:
; xdimg_dflat   
;  Version 1.1
;
; PURPOSE:
;    Creates dome flats given the direct image structure (assuming
;    there are dome flats in the structure)
;
; CALLING SEQUENCE:
;   
;   xdimg_dflat, struct, /IR, /SVOV, /SVDRKSUB, BPM=, OUTROOT=, /ALLOUT
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
;   /SVOV   - save ov files
;   OUTROOT - Root name of Dome flats (default is 'Flats/DFlat')
;   /ALLOUT - Output unnormalized frame too
;   BPM   - Name of the bad pixel mask to apply, including directory name.
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
;  XDIMG_DELCODIV
;  XDIMG_DELDRKSUB
;  X_FILTERS
;
; REVISION HISTORY:
;   18-July-2001 Written by JXP
;   24-Apr-2002 Allow for non-linearity
;   14-June-2007 Updated by LKP to allow for making flats for IR data
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro xdimg_dflat, struct, IR=ir, CCD=ccd, SVOV=svov, SVDRKSUB=svdrksub, BPM=bpm, OUTROOT=outroot, ALLOUT=allout

  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'xdimg_dflat, struct, /IR, /SVOV, /SVDRKSUB, BPM=bpm, OUTROOT=, /ALLOUT (v1.1)'
      return
  endif 
  
;  Optional Keywords
  if not keyword_set( OUTROOT ) then begin
     if keyword_set( CCD ) then begin
        CASE ccd of
           'NIRC2W': begin
              outroot = 'Flats/DFlat_w_'
           end
           'NIRC2N': begin
              outroot = 'Flats/DFlat_n_'
           end
        end
     endif else begin
        outroot = 'Flats/DFlat'
     endelse
     outroot = strtrim(outroot,2)
  endif

; Dealing with IR data
  if keyword_set(IR) then irdata=1 else irdata=0

;  Find the Dome Flats
  if keyword_set( CCD ) then dflts = where(struct.type EQ 'DFT' AND struct.ccd EQ ccd AND struct.flg_anly NE 0, ndflt) $
                        else dflts = where(struct.type EQ 'DFT' AND struct.flg_anly NE 0, ndflt) 
  if ndflt eq 0 then begin
     print, 'There are no dome flats to process.  Aborting!'
     return
  endif

  if irdata EQ 0 and strmid(struct[dflts[0]].ccd,0,5) eq 'NIRC2' then begin
     print, 'When dealing with IR data you must use the /IR flag.  Aborting!'
     return
  endif

; Preprocessing
  ;  Dark subtract if we're dealing with infrared CCDs.
  if irdata EQ 1 then begin
     preprocess = where(struct[dflts].flg_drksub eq 0, npreprocess)
     if npreprocess gt 0 then xdimg_darksub, struct, dflts[preprocess]
  endif else begin
     ; Overscan subtract if we're dealing with optical CCDs.
     ovflt = where(struct[dflts].flg_ov EQ 0, nov)
     if nov NE 0 then xdimg_over, struct, dflts[ovflt], ORDR=4
  endelse

;  Find all the filters involved
    x_filters, struct[dflts].filter, filt, nfilt

;  Loop on separate filters
    for q=0,nfilt-1 do begin
      wfilt = where(strtrim(struct[dflts].filter,2) EQ filt[q], dumi)

      outfil = strjoin([outroot, strtrim(filt[q],2), '.fits'])
      outfilN = strjoin([outroot, 'N', strtrim(filt[q],2), '.fits'])

      ; Status
      print, 'Combining images: '
      if irdata EQ 1 then begin
         for i=0, dumi-1 do print, struct[dflts[wfilt[i]]].img_drksub
      endif else begin
         for i=0,dumi-1 do print, struct[dflts[wfilt[i]]].img_ov
      endelse
      print, '             into ', outfil
      

      ; Combine Images
      ; Use dark subtracted images if dealing with IR data
      ; use overscan subtracted images if dealing with optical data.
      for hh=0, dumi-1 do begin
         if irdata EQ 1 then $
            struct[dflts[wfilt[hh]]].img_drksub=strtrim(struct[dflts[wfilt[hh]]].img_drksub,2) else $
               struct[dflts[wfilt[hh]]].img_ov=strtrim(struct[dflts[wfilt[hh]]].img_ov,2)
      endfor

      if irdata EQ 1 then $
         xcombine, struct[dflts[wfilt]].img_drksub, comb, head, SCALE=struct[dflts[wfilt]].med_drksub else $
            xcombine, struct[dflts[wfilt]].img_ov, comb, head, SCALE=struct[dflts[wfilt]].med_ov

      ; Non-linearity
      if not keyword_set( NONONLIN ) then $
        comb = xdimg_nonlinear(comb, struct[dflts[0]].ccd)

      ; Normalize
      ssec = struct[dflts[wfilt[0]]].statsec  ; Stat sec
      med = median(comb[ssec[0]:ssec[1],ssec[2]:ssec[3]])
      combN = comb / med
      
      ; Set bad pixels to a value of -1000. so we don't have problems when we flat field
      if keyword_set(BPM) then begin
         badpix = readfits(strtrim(BPM,2), /silent)
         bad = where(badpix eq 0.)
         combN[bad] = -10000.
      endif
      
      ; Header
      headN = head
      fxaddpar, headN, 'MEDFLT', med

      ; Output
      if keyword_set( ALLOUT ) then mwrfits, comb, outfil, head, /create
      mwrfits, combN, outfilN, headN, /create
      
   endfor


;   Delete images
    if irdata EQ 1 then begin
       if not keyword_set( SVDRKSUB ) then xdimg_deldrksub, struct, dflts
    endif else begin
       if not keyword_set( SVOV ) then xdimg_delov, struct, dflts
    endelse
    

; Resave the updated structure so that you can pick up where you left off easily.
    mwrfits, struct, 'struct.fits', /create
    
    print, 'All done with Dome Flats!'
    
 end
