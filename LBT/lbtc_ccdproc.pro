
;prcedure that reduce images from LBT/LBTC. A structure status keep track
;of what has been done and what not.	


;PLAN ---> plan file created by lbtc_makeplan
;SKYFIT ---> if set remove a second order polinomial surface
;gzip  --> if set gzip the final images
;path --> path with the data
;blueonly,redonly --> reduces only blue/red side
;archflat --> set to an array of name with archival flats if no flats are taken
;wcs  --> if set, run astrometry.net on reduced frame (this takes a lot)
;fitbias --> perform a linear interpolation across prescan/postscan
;            (not quite good results)


PRO lbtc_ccdproc, plan, skyfit=skyfit, gzip=gzip, path=path, blueonly=blueonly,$
                  redonly=redonly, archflat=archflat, wcs=wcs, fitbias=fitbias

  ;read the logfile
  readcol, plan, name, filt, time, obj, type, side, $
           format='A,A,F,A,A,A,A', /silent
  
  
  side=strtrim(side,2)

  ;get date
  date=strtrim(STRMID(name[0],7,6),2)

  ;set path 
  if ~keyword_set(path) then path='./Raw/'
   
  ;set log file
  splog, 'Start reduction at '+systime(), filename=date+'redux.log', /append


  ;make structures  on both sides
  
  ;;----
  currside='B'
  for sd=0, 1 do begin
      
     if keyword_set(blueonly) then sd=1
     
     if keyword_set(redonly) then begin
        currside='R'
        sd=1
     endif
        
     
  ;find science, bias and flat
     bias=where(type eq 'zer' and side eq currside, nbias)
     flat=where(type eq 'fla' and side eq currside, nflats)
     science=where(type eq 'obj' and side eq currside, nsci)
    
  ;-----   

  ;check if status str exists
     aa=file_search(currside+date+"status.str",count=cc)
     if (cc eq 0 ) then begin
   ;create status structure 
        splog, "Creating new STATUS structure for side ", currside
        status={BIASMED:0,FLAT:0,FILFLAT:STRARR(10),IMG:make_array(n_elements(name),/integer,value=0)}
        mwrfits, status, currside+date+"status.str", /create
     endif else begin
        splog, "STATUS structure found. Continuing reduction ", currside+date+"status.str"
        status=mrdfits(currside+date+"status.str",1,/silent)
     endelse

  ;;----
  ;switch to red side
     currside='R'
     
  endfor
  ;;----

;-----------------------------------------------------------------------
;make bias on both sides

  ;;------  
  ;loop sides
  currside='B'
  for sd=0, 1 do begin
     
     if keyword_set(blueonly) then sd=1
     
     if keyword_set(redonly) then begin
        currside='R'
        sd=1
     endif
     
      status=mrdfits(currside+date+"status.str",1,/silent)
     
     ;find bias
     bias=where(type eq 'zer' and side eq currside, nbias)
     
  ;;----

 

  ;check if bias already exists 
  
     if(status.biasmed eq 1) then  begin
        splog, "Using exisiting median bias."
        biasmedian=mrdfits(path+currside+date+"medbias.fits",/silent)
     endif   
     
     if(status.biasmed eq 2) then  splog, "No bias frame. Use overscan only!"
     
     if(nbias gt 0 and status.biasmed eq 0) then begin
        lbtc_makebias, name, biasmedian, bias, date=date, side=currside, $
                       path=path, fitbias=fitbias
        splog, "Done with median bias for side ", currside
        status.BIASMED=1
     endif 
  
     if(nbias eq 0 and status.biasmed eq 0) then begin
        splog, "No bias frame. Use overscan only"
        status.BIASMED=2
     endif   
     
     ;update structure
     mwrfits, status, currside+date+"status.str", /create

      ;switch to red side
     currside='R'
  endfor
  
  splog, 'All done with Bias!'

  
;---------------------------------------------------
;make flat for both sides

 ;;------  
  ;loop sides
  currside='B'
  for sd=0, 1 do begin
     
     if keyword_set(blueonly) then sd=1
     
     if keyword_set(redonly) then begin
        currside='R'
        sd=1
     endif
     
     status=mrdfits(currside+date+"status.str",1,/silent)
     
     ;find flat
     flat=where(type eq 'fla' and side eq currside, nflats)
     
  ;;----
  
  
   ;check if flats already exists 
  
     if(status.flat eq 1) then  splog, "Using exisiting flats."
     if(status.flat eq 2) and ~keyword_set(archflat) then begin
        splog, "No flats found..."
        splog, "Please, provide archival flats"
        return
     endif

     if(status.flat eq 0 and nflats GT 0) then begin
  
    ;set min and max values to accept a flat as good 
    ;(bias~1000, saturation 65535 ADU )
        minmax_flat=[1100.,64000.]
     
        ;lbtc_makeflats, name, filt, flat, status=status, bias=biasmedian, $
        ;                date=date, minmax=minmax_flat, side=currside, path=path
        ;Bias does not show big patterns. Do not subtract them
        lbtc_makeflats, name, filt, flat, status=status, $
                        date=date, minmax=minmax_flat, side=currside, path=path
        splog, "All done with flats"
        status.flat=1
     endif  

     if(status.flat eq 0 and nflats EQ 0) and ~keyword_set(archflat) then begin
        splog, "No flats found... I cannot reduce data without flats!"
        splog, "Please, provide archival flats"
        return
     endif
 
     if keyword_set(archflat) then begin
        splog, 'Using archival flats '+archflat 
        status.flat=2
     endif

  ;update structure
  mwrfits, status, currside+date+"status.str", /create
   
       ;switch to red side
     currside='R'
  endfor
  
  splog, 'All done with Flats!'



;------------------------------------------------------  
;make science frame
  
  ;loop sides
  currside='B'
  for sd=0, 1 do begin
     
     if keyword_set(blueonly) then sd=1
     
     if keyword_set(redonly) then begin
        currside='R'
        sd=1
     endif
     
     status=mrdfits(currside+date+"status.str",1,/silent)
     
     ;find science
     science=where(type eq 'obj' and side eq currside, nsci)
  ;;----
  

     if(nsci gt 0) then begin
     ;loop on each science frame
        for img=0, nsci-1 do begin
        
           if(status.img[img] eq 1) then  splog, "Found previously reduced frame ",$
                                                 name[science[img]] else begin
              
              ;lbtc_makescience, name, filt, science[img], status=status,bias=biasmedian,$
              ;                  fitbias=fitbias,skyfit=skyfit, gzip=gzip,date=date, $
              ;                  side=currside, wcs=wcs, path=path, object=obj
              ;Bias shows no large scale pattern. Do not subtract them
              lbtc_makescience, name, filt, science[img], status=status,$
                                fitbias=fitbias,skyfit=skyfit, gzip=gzip,date=date, $
                                side=currside, wcs=wcs, path=path, object=obj
      
                                
     ;update structure
              status.img[img]=1
              mwrfits, status, currside+date+"status.str", /create
           endelse
        endfor
     endif else begin
        splog, "No science frame for ", currside
        return
     endelse
     
     
     ;switch to red side
     currside='R'
  endfor
  

  splog, 'All done with data reduction at '+systime(), /close
     



END
