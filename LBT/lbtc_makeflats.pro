
;procedure that make flats from flat frames 


;name---> list of file names
;filter --> list of filters
;index ---> array of indexes which contains only flats frame
;date-->  observation date
;bias --> the median bias from makebias if it exists
;MINMAX ---> array of min and max value to esclude faint or saturated flats  
;STATUS ---> the status structre from ccdproc 

PRO lbtc_makeflats, name, filter, index,  bias=bias, status=status, date=date, minmax=minmax, side=side, path=path
 


  ;loop over flats to find group of common filters 
  
  for i=0, n_elements(index)-1 do begin
     if(i eq 0) then flat=strcompress(filter[index[i]],/remove_all)
     if(i gt 0) then begin
        ;check if new filter
        res=where(flat eq strcompress(filter[index[i]],/remove_all),nummatch)
        ;if new, append
        if(nummatch eq 0) then flat=[flat,strcompress(filter[index[i]],/remove_all)]
     endif
  endfor
 
  splog, "For side ", side, " found filters: ", flat 
  
 ;make flats
  for flfra=0, n_elements(flat)-1 do begin
  ;check if flat done
     sss=where(strcompress(status.filflat,/remove_all) eq flat[flfra],numbfl)
     
     if(numbfl gt 0) then splog, "Found flat ", flat[flfra] else begin
        splog, "Working on filter ", flat[flfra]

     ;find images for this filter
        thisfilter=where(strcompress(filter[index],/remove_all) eq flat[flfra],nfilt)
     
    ;check if everything is ok (nfilt MUST be > 0)
        if(nfilt eq 0) then begin
           splog, "I cannot find images with this filter: ", flat[flfra], " Returning..."
           return
        endif
     
     ;call loop flats to make final image

        lbtc_loopflats, name, nfilt, thisfilter, index, finalflat, $
                        bias=bias, minmax=minmax, path=path
        
     ;save final output
        mwrfits, finalflat, path+side+date+"_"+flat[flfra]+"flat.fits", /create
     ;update status structure
        status.filflat[flfra]=flat[flfra]
        mwrfits, status, side+date+"status.str", /create
     endelse
  endfor
  
  splog, "All done with flats for side ",side 
  
   
  
END
