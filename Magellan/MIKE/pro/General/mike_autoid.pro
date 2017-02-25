;+ 
; NAME:
; mike_autoid   
;     Version 2.0
;
; PURPOSE:
;    Guess the type of MIKE image based on exposure time, counts,
;    and other characteristics of the image.  While this algorithm is
;    quite clever, there are many times when it makes the wrong
;    guess.  As such, you must be careful to check the output types.
;
; CALLING SEQUENCE:
;   
;  guess = mike_autoid(data, xbinguess, ybinguess, $
;               filename=img[q], hdr=head, /silent)
;
; INPUTS:
;    data      - 2D image (generally set to 0 and read from filename)
;
; RETURNS:
;  Image type::  'UNK', 'BAD', 'ZRO', 'ARC', 'IFLT', 'TFLT', 'STD',
;  'OBJ', 'MFLT', IFLT', 'TWI'
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   FILENAME   - Name of the image to be considered
;   XBIN       - x binning
;   YBIN       - y binning
;   SILENT     - Suppress print statements
;   SATURATED  - Saturation level (default = 50000.)
;   
;
; OPTIONAL OUTPUTS:
;   HDR     -  Image header if filename specified
;   GUESS_EXPTIME -- Guess of exposure time based on CR hits
;
; COMMENTS:
;   Need to do some consistency checks between red and blue side,
;         i.e. they should be the same type
;
; EXAMPLES:
;   guess = mike_autoid(data, xbin, ybin, filename='file')
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   ??-??-2003 Written by SB
;-
;------------------------------------------------------------------------------
function mike_autoid, image, xbin, ybin, filename=filename, silent=silent, $
       saturated=saturated, ilun=ilun, guess_exptime=guess_exptime, hdr=hdr

  ;;
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'guess = mike_autoid(image, [xbin, ybin], FILENAME=, /SILENT,'+ $
        'SATURATED=, ILUN=, GUESS_EXPTIME=, HDR=, (v2.0)'
      return, -1
  endif 

    if keyword_set(filename) then mike_mproc, filename, image, hdr=hdr, /silent
    if NOT keyword_set(saturated) then saturated = 50000.0
    if NOT keyword_set(ilun) then ilun =-1
    if NOT keyword_set(silent) then s=1 else s=0
    
    saturated_mask = (image GT saturated)
    ncol = (size(image))[1]
    nrow = (size(image))[2]
 
    exptype = 'UNK'
;
;       central ~10% here:
;
    img = (image * (saturated_mask EQ 0))[ncol/3:2*ncol/3,nrow/3:2*nrow/3]
    nsat = total(saturated_mask)

        
;
;     guessing here, are these pre- or postpost overscan removal?
;     let's say post!
;
    if NOT keyword_set(xbin) then xbin = long(2200.0 / ncol)
    if NOT keyword_set(ybin) then ybin = long(4500.0 / nrow)

;
;   Let's compare three images, the original img
;   A median-filtered image along rows:      imgrow
;   A median-filtered image along columns:   imgcol
;
    mn = mean(img)
    if s then printf,ilun, filename, '= ', format='(a,a,$)'

    if mn LT -10.0 then begin
        if s then printf, ilun, 'Negative counts...bad exposure'
        exptype = 'BAD'
        return, exptype
    endif

    ;; Check for Bias  (low counts)
    ;; JXP -- Error in bad images 
    p75 = where(img GE mn, npx)
    if npx LT 10 then begin
        print, 'mike_autoid: Image has all the same value. Bad readout..'
        exptype = 'BAD'
        return, exptype
    endif

    m75 = mean(img[p75])
    if mn GT -5 AND mn LT 5 AND m75 LT 5 then begin
         if s then printf,ilun, '... This must be a bias'
         return, 'ZRO'
    endif

    his_mn = histogram(img[p75], min=mn, bins=10.)
    his_bn = findgen(n_elements(his_mn)) + mn
    his_cumul = total(his_mn, /cumul)
    p90 = min(where(his_cumul GT npx*0.8)) > 0
    if his_bn[p90] GT 16000. then $
       printf, ilun, 'high counts at 90% level: ', his_bn[p90]



    kernel = fltarr(2*(16L/xbin)+1) + 1
    kernel = kernel/total(kernel)

    imgrow = median_row(img, 2*(16L/xbin)+1, /reflect)
    imgt = median_row(transpose(img), 2*(16L/ybin) + 1, /reflect)
    imgcol = transpose(imgt)

    mns = [mn,mean(imgrow),mean(imgcol)]
    s25 = [mean(img[where(img LT mns[0])]), $
           mean(imgrow[where(imgrow LT mns[1])]), $
           mean(imgcol[where(imgcol LT mns[2])])]
    s75 = [mean(img[where(img GT mns[0])]), $
           mean(imgrow[where(imgrow GT mns[1])]), $
           mean(imgcol[where(imgcol GT mns[2])])]

    if s then printf,ilun, mns, s25, s75 

    coldiff = img - imgcol
    rowdiff = img - imgrow
    rowcol  = imgrow - imgcol

    sx = (size(coldiff))[1]    
    sy = (size(coldiff))[2]    
    sxmin = 16L/xbin
    sxmax = sx-1-16L/xbin
    symin = 16L/ybin
    symax = sy-1-16L/ybin
 
    meandiff = [mean(coldiff[sxmin:sxmax, symin:symax]), $
                mean(rowdiff[sxmin:sxmax, symin:symax]), $
                mean(rowcol[sxmin:sxmax, symin:symax])]
    if s then printf,ilun,meandiff

    crmask = where(coldiff[sxmin:sxmax, symin:symax] GT $
                0.2*imgcol[sxmin:sxmax, symin:symax]+100 AND $
                   rowdiff[sxmin:sxmax, symin:symax] GT $
                    imgrow[sxmin:sxmax, symin:symax]*0.15, ncr)

    guess_exptime = 4*(ncr * xbin * ybin - 400)/10.  > 0 ; * 3 if red
    smooth_operator = meandiff[2]*abs(meandiff[2])/(abs(mns[0])+20.0)

    if s then printf,ilun,ncr, guess_exptime, smooth_operator


    ; print, smooth_operator, mns[0], meandiff
    if abs(smooth_operator) LT 100.0 AND $
       abs(smooth_operator) GT 1.0  AND $
       abs(mns[0]) GT 500. then begin
       if abs(meandiff[0]) LT 3.0 AND $
             abs(meandiff[1]-meandiff[2]) GT 80.0*abs(meandiff[0]) then begin
         if s then printf,ilun, '... this is an internal Quartz Lamp'
         return, 'IFLT'
       endif
    endif
    if smooth_operator GT 0.5 AND abs(mns[0]) GT 10 then begin
       if s then printf,ilun, 'Calibration exposure? ', format='(a,$)'
       if total(meandiff[0:1]) GT 8.0 then begin
         if s then printf,ilun, '... this is an Arc'
         return, 'ARC'
       endif
    endif
    if smooth_operator GT 2.0 AND abs(mns[0]) GT 15 then begin
       if total(meandiff[0:1]) LT -5.0 then begin
         if s then printf,ilun, '.. this is daylight blue sky spectrum'
         return, 'TFLT'
       endif
    endif 


    if smooth_operator -1.0 GT -0.001 * mns[0] then begin
         if s then printf,ilun, '... Diffuse spectrum (sky, internal, or O-star?)'
         return, 'MFLT'
    endif 


    if smooth_operator LT -2.0 then begin
       if s then printf,ilun, 'Object is brighter than sky'

       if meandiff[0] LT -10.0 then begin
         if s then printf,ilun, '.. this is daylight red sky spectrum'
         return, 'TFLT'
       endif

       if mns[0] GT 100.0 AND guess_exptime LT 1000.0 then begin
          if s then printf,ilun, '...Standard star, or a very high S/N science exposure'
          return, 'STD'
        endif else begin
          if s then printf,ilun, '...Science target, good counts'  
          return, 'OBJ'
        endelse

    endif else begin
       if s then printf,ilun, 'Smooth spectrum ', format='(a,$)'
       if max(s25) GT 2000.0 then begin
         if s then printf,ilun, '... Diffuse spectrum (sky, internal, or O-star?)'
         return, 'MFLT'
       endif 

       if meandiff[0] LT -1.0 AND meandiff[1] GT 10.0 then begin
         if s then printf,ilun, '... Low signal internal quartz?'
         return, 'IFLT'
       endif 


       if meandiff[1] LT 0.0 AND meandiff[1] LT 0.0 then begin
         if s then printf,ilun, '... Low signal night sky?'
         return, 'TWI'
       endif else begin
         if s then printf,ilun, '... Science target, low counts?'
         return, 'OBJ'
       endelse


    endelse

    return, 'UNK'

end
    
