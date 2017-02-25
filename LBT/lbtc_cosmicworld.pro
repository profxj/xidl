
;pro that deals with cosmic rays


pro lbtc_cosmicworld, infits, mask, creject=creject, currside=currside, $
                      currchip=currchip, medianimg=medianimg, scale=scale,$
                      shift=shift, objectpix=objectpix
  



 ;;;;;;;;;;;;BE CAREFUL WHEN ADDING BAD MASK .. problem with the edges
  ;;;;;;   FOR COSMIC, COPY ISOAREA STUFF FROM LRIS
 

  ;----------------------------------------
  ;;Cosmics plus masks
  ;----------------------------------------

  ;now deal with cosmics
  splog, 'Mask cosmic rays...'
 
  ;set tolerance 
  if ~keyword_set(creject) then creject=6.

  ;compare median and image, to flag outliers 
  diff_image=medianimg-infits*scale
  djs_iterstat, diff_image, maxiter=5, median=mean_level, sigma=sig_level
  flag_cosmic=where(abs(diff_image) GT mean_level+creject*sig_level,cfnd)
     
  ;mask all the rejected
  if(cfnd gt 0) then begin
     splog, 'Found ', cfnd, ' cosmic rays'
     mask[flag_cosmic]=mask[flag_cosmic]*0.
  endif


   ;also open bad column masks
   ;if on red side, get the achived mask mask
   if(currside eq 'B') then $
      mas_arch=mrdfits(getenv('LBT_IMG')+'/archive/Bbadpix2455231.3_.fits')
   if(currside eq 'R') then $
      mas_arch=mrdfits(getenv('LBT_IMG')+'/archive/Rbadpix2455231.3_.fits')
 
   ;get right chip
   archmask=reform(mas_arch[currchip,*,*])
   undefine, mas_arch
  
 
  ;;LBT filed of view is huge... To clean object by object is not ideal.
  ;;Use a object mask obtained with the median, instead.
 
   mask[objectpix]=mask[objectpix]*0.+1.

  ;;---
  ;;;restore pixel within aperture  around bright stars
  ;;;use the previous extracted sources
  ; splog, 'Restore star centers...(it takes a bit)'
  ; 
  ; ;read once again star list
  ; readcol, 'tmp.dat', x_ref, y_ref, a1, b1, c1, d1, e1, f1, g1,$
  ;         h1, i1, isoarea, k1, l1, class_star, xg1, yg1,$
  ;         format='F,F,I,F,F,F,F,F,F,F,F,I,I,F,F,F,F', /silent
  ;
  ; 
  ;;take the 60 brightes stars (file is already sorted)
  ; nsextra=n_elements(x_ref)
  ; if(nsextra gt 60) then nsextra=60
  ;
  ; axis_x=n_elements(mask[*,0])
  ; axis_y=n_elements(mask[0,*])
  ; 
  ; for star=0, nsextra-1 do begin  
  ;    rad_sta=1.5*sqrt(isoarea[star]/!PI)
  ;    pval=mpix_circ(x_ref[star],y_ref[star],rad_sta,XSIZE=axis_x,YSIZE=axis_y, count=nunf,/silent)
  ;    if(nunf gt 0) then   mask[pval]=mask[pval]*0.+1
  ; endfor
   
         
   ;add bad columns
   ;shift bad column mask accordingly
   shifmask=shiftf(archmask,shift[0],shift[1])
   undefine, archmask
   
   ;correct mask (cut the interpolation at 0.5)
   bignum=where(shifmask gt 0.5, numbig)
   if(numbig gt 0) then shifmask[bignum]=shifmask[bignum]*0+1 else shifmask=shifmask*0
   
   ;now add bad columns 
   mask=temporary(mask)*shifmask
   undefine, shifmask


end
