
;;procedure that gets a lbtc mosaic in inputs and makes a final fits
;;for science data. For test flats and focus frame, just disply the
;;image how it is and do not create any structure/mosaic.



;imagename: fits file as it comes from the telescope
;imgstruc:  if set to a variable, in output a structure with each
;           chip+header 
;mosaic:    if set to a variable, in output a single array full mosaic
;           (not included fractional rotation)
;nodisplay:   if keyword set, no xatv will be launched
;save        if keyword set create fits file of mosaic and structure
;redux       if set, display the image structure that as been
;             processed by lbtc_ccdproc



;written by MF Sept 2009



PRO lbtc_readfits, imagename, IMGSTRUC=imgstruc, MOSAIC=mosaic, $
                   NODISPLAY=nodisplay, SAVE=save, REDUX=redux



  if ~keyword_set(redux) then begin
     ;treat this as a raw image

;size chips 2304*4608
;prescan 0:50, overscan 2099:2304
;data [51:2098,1:4608]

;check if science frame
     science=1
     hed0=headfits(imagename,exten=0)
     numext=sxpar(hed0,'NEXTEND')
     if (numext eq 1) then science=0 else science=1
     

;no science
     if (science eq 0) then mosaic=mrdfits(imagename,1,/fscale)

;science
     if (science eq 1) then begin
 ;set structure
        fit=make_array(2048,4608,/double)
        pres=make_array(50,4608,/double)
        overs=make_array(206,4608,/double)
        hea=make_array(1D4,/string)
        imgstruc={header:hea,data:fit,prescan:pres,overscan:overs}  
        imgstruc=replicate(imgstruc,5)
   
   
;iterate over levels
        for i=0, 4 do begin
;open file and load in structure
           fits=mrdfits(imagename,i,hea,/silent,/fscale)
           imgstruc[i].header=hea
           if (i ne 0) then begin
              imgstruc[i].data=fits[50:2097,*]
              imgstruc[i].prescan=fits[0:49,*]
              imgstruc[i].overscan=fits[2098:2303,*]
           endif
        endfor
   
     endif
 
  endif




  if keyword_set(redux) then begin
     ;treat this as a reduced image
     
     ;open the header
     fits=mrdfits(imagename,0,hea,/silent)
     fit=make_array(2048,4608,/float)
     imgstruc={header:hea,data:fit}  
     imgstruc=replicate(imgstruc,5)
     
     ;get the chips
     for i=1, 4 do begin
        imgstruc[i].data=fits[i-1,*,*]
     endfor
  endif

;create fits mosaic (6178*6673)
  mosaic=make_array(6178,6673)
   
   
;fill in third chip
  mosaic[0:2047,0:4607]=imgstruc[3].data
;fill in second chip, leaving 18 pixx for chip gap
  mosaic[2065:4112,0:4607]=imgstruc[2].data
;fill in first chip, leaving 18 pixx for chip gap
  mosaic[4130:6177,0:4607]=imgstruc[1].data
;rotate and fill in 4th chip,  leaving 18 pixy for chip gap
  chip4=rotate(imgstruc[4].data,1)
  mosaic[770:5377,4625:6672]=chip4

;save
   if keyword_set(save) then begin
      mosname="mos_"+imagename
      strname="str_"+imagename
      mwrfits, mosaic, mosname, imgstruc[0].header, /create
      mwrfits, imgstruc, strname, /create
   endif
   



;display
if ~keyword_set(nodisplay) then atv, mosaic


end
