;+ 
; NAME:
; wfc3_g280_find_sources
;
; PURPOSE:
;  Test program that finds sources in the direct image to reduce
;  in the spectral image. This program has NOT been tested very
;  well, and should not be used blindly.
;
; CALLING SEQUENCE:
;   wfc3_g280_find_sources, wfc3_g280_strct, fraction=, $
;                           npix=
;
; INPUTS:
;   wfc3_g280_strct -- the wfc3_g280 structure
;
; RETURNS:
;
; OUTPUTS:
;   Update structure with guess positions of possible sources
;
; OPTIONAL KEYWORDS:
;   FRACTION= -- the fraction that sets the threshold level to be
;                considered a detection
;   NPIX= -- the minimum number of pixels in an object to be
;            considered a real feature (not a hot pixel, etc.)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;  wfc3_g280_find_sources, wfc3_g280_strct, fraction=fraction, $
;                          npix=npix
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   10-Jun-2016 Written by MN
;------------------------------------------------------------------------------

pro wfc3_g280_find_sources, wfc3_g280_strct, FRACTION=fraction, $
                            NPIX=npix
  
  if size(fraction,/type) EQ 0 then fraction = 0.001
  if size(npix,/type) EQ 0 then npix = 200
  
  cnt=0L
  for ii=0,1 do begin
     img = xmrdfits(wfc3_g280_strct(0).img_fil, 3*ii+1L, head)
     
     s = size(img)
     a = median(img,5)
     j = sort(a)
     level = a[j[s[4]*(1.0-fraction)]]
     

     b = label_region(a GT level, /all)  ; Get blob indices above level value set by fraction.
     h = histogram(b, REVERSE_INDICES=r) ; Get population and members of each blob.
     tmp=where(h gt npix, nsources)       ; npix above the threshold level
     for jj=1L,nsources-1 do begin
        source = (reverse(sort(h)))[jj]   ; Find the nth largest source (after the sky)
        ind = r[r[source]:r[source+1]-1] ; Find subscripts of the source pixels
        tmp=max(img(ind),maxpix)         ; locate the maximum pixel
        wfc3_g280_strct[cnt].xguess = double(ind(maxpix) MOD s[1]) ; locate x_guess
        wfc3_g280_strct[cnt].yguess = double(ind(maxpix)/s[1])     ; locate y_guess
        wfc3_g280_strct[cnt].chip = 2-ii                           ; chip number
        cnt=cnt+1
     endfor
  endfor

  wfc3_g280_strct=wfc3_g280_strct(0:cnt-1)
  
end
