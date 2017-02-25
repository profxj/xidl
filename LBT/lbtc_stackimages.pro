;look for pair of images in the science frame and stack images with
;same filters

;NOTE nocleancr is disabled 

;gzip  --> if set, gzip images
;nowcs   --> if set, no wcs are added to the field.
;[blue/red]only --> reduces only the red or blue side            
;chipnum --> the number of the chip to be reduced [0,1,2,3]
;noskysub --> if set, sky background is not subtracted
;creject  -- > set the N of sigma to reject cosmic rays in median comparison
;manshift --> do the manual shift selecting stars

PRO lbtc_stackimages, path=path, gzip=gzip, nowcs=nowcs, blueonly=blueonly,$
                      redonly=redonly, chipnum=chipnum,$
                      noskysub=noskysub, creject=creject,$
                      manshift=manshift

  splog, 'Start stacking on '+systime(), filename='stack_image.log', /append
  
  if ~keyword_set(path) then path='./Raw/'

  ;make sci folder
  spawn, 'mkdir Sci'



  ;loop on red and blue sides separately

  ;;----
  currside='B'
  for sd=0, 1 do begin
      
     if keyword_set(blueonly) then sd=1
     
     if keyword_set(redonly) then begin
        currside='R'
        sd=1
     endif

  ;loop on  chips
     
     if keyword_set(chipnum) then begin
        if(chipnum gt 3) then begin
           splog, 'chipnum must be less than 3!'
           return
        endif
        currchip=chipnum 
     endif

     currchip=0
     for ch=0, 3 do begin
      

        if keyword_set(chipnum) then begin
           currchip=chipnum 
           ch=3
        endif
        
   
  ;get the list of reduced frames
        if(currside eq 'B') then $
           spawn, 'ls '+path+'lbcb.*_redux.fits', img_name else $
              spawn, 'ls '+path+'lbcr.*_redux.fits', img_name  
            
        numimg=n_elements(img_name)

        splog, 'Found ',numimg,' images on side ',currside,' chip ',currchip 

        ;allocate info arrays
        img_obje=strarr(numimg)
        img_filt=strarr(numimg)

        for i=0, numimg-1 do begin
        ;acquire info
           header=headfits(img_name[i]) 
           img_obje[i]=sxpar(header,"OBJNAME",/silent)
           img_filt[i]=sxpar(header,"FILTER",/silent)
        endfor
        
        ;now process
        lbtc_stc_process, img_name, img_obje, img_filt, numimg, path=path, $
                          currside=currside, currchip=currchip, noskysub=noskysub, $
                          creject=creject, manshift=manshift

       
        ;move to next chip chip
        currchip=currchip+1
     endfor
  ;;----

 ;switch to red side
     currside='R'
  endfor
  ;;----

  
   ;fit new wcs
   ;if ~keyword_set(nowcs) then lris_fitwcs, save_sci_name, $
   ;                                      currchip=currchip

   ;gzip
    
           ;   if keyword_set(gzip) then begin
           ;   spawn, 'gzip '+save_sci_name
           ;   spawn, 'gzip '+save_mask_name
           ;   spawn, 'gzip '+save_median
           ;endif
  
        
  
  ;clean 
    spawn, 'rm -fr defaul.sex'
    spawn, 'rm -fr defaul.param'
    spawn, 'rm -fr defaul.nnw'
    spawn, 'rm -fr defaul.conv'
    spawn, 'rm -fr tmp.dat'
    spawn, 'rm -fr check.fits'
    

    splog, /close
    
 END









