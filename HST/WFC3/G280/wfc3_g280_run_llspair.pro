;+ 
; NAME:
;  wfc3_g280_run_llspair
;
; PURPOSE:
;   wrapper program for the LLS pair program
;
; CALLING SEQUENCE:
;   
;   wfc3_g280_run_llspair, list, IMG_ROOT=img_root, BADPIX=badpix, BOXCAR=boxcar, $
;                           USEUPPERSKY=useuppersky, USELOWERSKY=uselowersky, NOFIT=nofit
;
; INPUTS:
;
;   list -- list with objects to reduce
; 
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  BADPIX= -- Threshold to mark bad pixels in direct image
;  BOXCAR= -- Boxcar used for extraction
;  USEUPPERSKY= -- Use upper sky for sky subtraction
;  USELOWERSKY= -- Use lower sky for sky subtraction
;  NOFIT= -- Do not run the centroid routine on the direct image
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
; wfc3_g280_run_llspair, list, IMG_ROOT=img_root, BADPIX=badpix, BOXCAR=boxcar, $
;                           USEUPPERSKY=useuppersky, USELOWERSKY=uselowersky, NOFIT=nofit
;
; PROCEDURES CALLED:
; adxy
; wfc3_g280_reduce_qso
;
; REVISION HISTORY:
;   10-Jun-2016 Written by MN
;------------------------------------------------------------------------------

pro wfc3_g280_run_llspair, list, IMG_ROOT=img_root, BADPIX=badpix, BOXCAR=boxcar, $
                           USEUPPERSKY=useuppersky, USELOWERSKY=uselowersky, NOFIT=nofit, $
                           _extra=_extra
  
  if not keyword_set(IMG_ROOT) then img_root='clean/final/'
  
  ;; read the list
  readcol, list, name, img_fil, spec_fil, img_ext, spec_ext, name1, ra1, dec1, name2, ra2, dec2, $
           comment='#', format='A,A,A,A,A,A,A,A,A,A,A'
  
  for ii=0L, n_elements(ra1)-1 do begin
     ;; read the headers
     h1=headfits(img_root+name(ii)+'_'+img_ext(ii)+'.fits.gz',ext=4)
     h2=headfits(img_root+name(ii)+'_'+img_ext(ii)+'.fits.gz',ext=1)
     dim=[sxpar(h1,'NAXIS1'),sxpar(h1,'NAXIS2')]
     
     ;; convert ra and dec to decimal
     ra=[ra1(ii),ra2(ii)]
     dec=[dec1(ii),dec2(ii)]
     ra=hms2dec(ra)*15.
     dec=hms2dec(dec)
     
     guess=dblarr(3,n_elements(ra))
     for jj=0L, n_elements(ra)-1 do begin
        ;; find their position on the chips
        adxy, h1, ra(jj), dec(jj), x1, y1
        adxy, h2, ra(jj), dec(jj), x2, y2

        if (x1 ge 0 and x1 lt dim(0)) and $
           (y1 ge 0 and y1 lt dim(1)) $
        then $
           guess(*,jj)=round([x1,y1,1]) $
        else $
           if (x2 ge 0 and x2 lt dim(0)) and $
              (y2 ge 0 and y2 lt dim(1)) $
           then $
              guess(*,jj)=round([x2,y2,2]) $
           else $
              stop
     endfor
     
     wfc3_g280_reduce_qso, img_root+name(ii)+'_'+img_ext(ii)+'.fits.gz', $
                           img_root+name(ii)+'_'+spec_ext(ii)+'.fits.gz', $
                           guess=guess, srch=15, badpix=badpix, boxcar=boxcar, $
                           useuppersky=useuppersky, uselowersky=uselowersky, $
                           nofit=nofit, _extra=_extra
  endfor

end

        
