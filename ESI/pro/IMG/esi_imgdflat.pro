;+ 
; NAME:
; esi_imgdflat   
;  Version 1.1
;
; PURPOSE:
;    Creates dome flats given the image list structure
;
; CALLING SEQUENCE:
;   
;   esi_imgdflat, struct, /SVOV, OUTROOT=
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
;   esi_imgdflat, nght1_strct
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
;   30-July-2002 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro esi_imgdflat, esi, SVOV=svov, OUTROOT=outroot, ALLOUT=allout

;
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'esi_imgdflat, esi, /SVOV, OUTROOT=, ALLOUT= (v1.1)'
      return
  endif 
  
;  Optional Keywords

  if not keyword_set( OUTROOT ) then outroot = 'Flats/FlatIMG_D'
  outroot = strtrim(outroot,2)
  
;  Find the Dome Flats

  flat = where(esi.type EQ 'DFLT' AND esi.flg_anly NE 0 $
                AND esi.mode EQ 0, ndflt)

;  Overscan

  if not keyword_set( REDOOV ) then bias = where(esi[flat].flg_ov EQ 0, nbias) $
  else begin
      nbias = n_elements(flat)
      bias = lindgen(nbias)
  endelse
  if nbias NE 0 then esi_subbias, esi, flat[bias]

;  Find all the filters involved
  
  x_filters, esi[flat].imfilt, filt, nfilt
  
;  Loop on separate filters
  
  for q=0,nfilt-1 do begin
      wfilt = where(esi[flat].imfilt EQ filt[q], dumi)

      outfil = strjoin([outroot, strtrim(filt[q],2), '.fits'])

      ;; Status
      print, 'Combining images: '
      for i=0,dumi-1 do print, esi[flat[wfilt[i]]].img_ov
      print, '             into ', outfil

      ;; Combine Images
      xcombine, esi[flat[wfilt]].img_ov, comb, head, scale='MED'

      ;; Non-linearity
;      if not keyword_set( NONONLIN) then $
;        comb = esi_imgnonlinear(comb, esi[flat[0]].ccd)

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
  
;  Delete images

  if not keyword_set( SVOV ) then esi_delov, esi, flat

  print, 'esi_imgdflat: All done with Dome Flats!'

end
