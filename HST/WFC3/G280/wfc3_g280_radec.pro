;+ 
; NAME:
;  wfc3_g280_radec
;
; PURPOSE:
;   Assigns the ra and dec of the object and creates the name from
;   this information.
;
; CALLING SEQUENCE:
;   
;   wfc3_g280_radec, wfc3_g280_strct
;
; INPUTS:
;
;   wfc3_g280_strct -- the wfc3_g280 structure
; 
; RETURNS:
;
; OUTPUTS:
;   Updates structure with ra, dec, and name
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
; wfc3_g280_radec, wfc3_g280_strct
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   10-Jun-2016 Written by MN
;------------------------------------------------------------------------------

;; from x and y values of the centroid image, it calculates the
;; ra and dec and assigns a name based on this

pro wfc3_g280_radec, wfc3_g280_strct

  for ii=0L,n_elements(wfc3_g280_strct)-1 do begin
     head=headfits(wfc3_g280_strct(ii).img_fil,exten=(7-3*wfc3_g280_strct(ii).chip))

     ;; calculate ra and dec
     xyad, head, wfc3_g280_strct(ii).x0, wfc3_g280_strct(ii).y0, ra, dec

     ;; make a name
     if dec lt 0 then sign='-' else sign='+'
     name='J'+strcompress(dec2hms(ra/15.),/r)+sign+$
          strcompress(dec2hms(abs(dec)),/r)
     
     ;; assign to the structure
     wfc3_g280_strct(ii).ra=ra
     wfc3_g280_strct(ii).dec=dec
     wfc3_g280_strct(ii).name=name

  endfor  
end
