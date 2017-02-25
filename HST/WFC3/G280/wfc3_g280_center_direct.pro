;+ 
; NAME:
; wfc3_g280_center_direct
;
; PURPOSE:
;   Apply a simple algorithm to center the object in the direct
;   image.  This is necessary to generate the trace and wavelength
;   solutions. 
;
; CALLING SEQUENCE:
;  wfc3_g280_center_direct, wfc3_g280_strct, EXTENDBOX=, SRCH=, $
;                             BADPIX=
; INPUTS:
;   wfc3_g280_strct -- the wfc3_g280 structure
;
; RETURNS:
;
; OUTPUTS:
;   x0 -- Best centered x position
;   y0 -- Best centered y position
;
; OPTIONAL KEYWORDS:
;   EXTENDBOX= -- Number of pixels to extend analysis for object
;   SRCH=  -- Number of pixels to extend search for QSO 
;   BADPIX= -- Pixel threshold above which the pixel values are
;              ignored for fitting, useful for CR+hot pixel rejection
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;  wfc3_g280_center_direct, wfc3_g280_strct, EXTENDBOX=extendbox, SRCH=srch, $
;                             BADPIX=badpix
;
; PROCEDURES CALLED:
;  cntrd
;
; REVISION HISTORY:
;   23-Dec-2010 Written by JXP/JMO
;   10-Jun-2016 Update to deal with structure input and improved
;               search algorithm to deal with low S/N data by MN
;------------------------------------------------------------------------------

pro wfc3_g280_center_direct, wfc3_g280_strct, EXTENDBOX=extendbox, SRCH=srch, $
                             BADPIX=badpix

  if (N_params() LT 1) then begin 
    print,'Syntax - ' + $
          'wfc3_g280_center_direct, img_fil, spec_fil, fin_strct, QADIR=, NAME=, ' + $
          'XGUESS=, YGUESS=, SEARCH= [v1.0]'
    return
  endif 

  for ii=0L, n_elements(wfc3_g280_strct)-1 do begin
  ;;;;;;;;;;;
     ;; Direct Image
     img = xmrdfits(wfc3_g280_strct(ii).img_fil, $
                    (7-3*wfc3_g280_strct(ii).chip), head)
     sz_img = size(img, /dim)

     ;; trim the image
     img_trim=img[wfc3_g280_strct(ii).xguess-srch:wfc3_g280_strct(ii).xguess+srch, $
                  wfc3_g280_strct(ii).yguess-srch:wfc3_g280_strct(ii).yguess+srch]
     
     ;; lets tweak the image here for bad pixels that could dominate
     ;; if only a single exposure was taken (CR)
     ;; i.e. remove bad pixels with high (greater than badpix) value
     img_trim(where(img_trim gt badpix))=0d
     
     ;; remove small structures i.e. hot pixels, CRs, etc
     ;; this will only do something if the search FOV has
     ;; 2 sources in it, then it will choose the largest
     ;; 3-sigma source and mask out the rest.
     
     his=histogram(img,bin=1,min=-300,max=300,loc=loc)
     tloc=loc(where(his gt 0.5*max(his)))
     sig=(max(tloc)-min(tloc))/2.35482d
     b = label_region(img_trim GT 3*sig, /all) ; Get blob indices.
     h = histogram(b, REVERSE_INDICES=r)       ; Get population and members of each blob.
     if n_elements(h) ge 3 then begin
        source = (reverse(sort(h)))[1]         ; Find the largest source after the sky
        ind = r[r[source]:r[source+1]-1]       ; Find subscripts of the source pixels
        mask=img_trim                          ; Mask out the rest of the image
        mask(*)=0d
        mask(ind)=1d
        img_trim=mask*img_trim
        img(*)=0
        img[wfc3_g280_strct(ii).xguess-srch:wfc3_g280_strct(ii).xguess+srch, $
            wfc3_g280_strct(ii).yguess-srch:wfc3_g280_strct(ii).yguess+srch] = img_trim
     endif
     
     ;; Find maximum near the guess
     mx = max(img_trim,imx)
     xi = wfc3_g280_strct(ii).xguess - srch + (imx mod (2*srch+1))
     yi = wfc3_g280_strct(ii).yguess - srch + (imx)/(2*srch + 1)
     
     ;; Centroid
     cntrd, img, xi, yi, x0, y0, 5., EXTENDBOX=extendbox
     
     ;;assign
     wfc3_g280_strct(ii).x0=x0
     wfc3_g280_strct(ii).y0=y0 
     
  endfor
     
  return
end
