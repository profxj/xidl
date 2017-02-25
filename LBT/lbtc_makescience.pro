
;procedure to make science frame (one each time)
;bias,zero,flat,gain,exptime


;str---> structure from lbtc make_plan 
;index ---> index of the image current under exam
;bias --> the median bias from makebias if it exists
;FITBIAS --> perform a linear interpolation across prescan/postsca. If
;            not set only median of 0.5*(postcan+prescan) is
;            subtracted from each row
;STATUS ---> the status structre from ccdproc 
;PLAN---> what comes out of makeplan

PRO lbtc_makescience, name, filter, index, status=status, bias=biasmedian,$
                      fitbias=fitbias, skyfit=skyfit, gzip=gzip,date=date, $
                      side=side, wcs=wcs, path=path, object=object

  splog, "Processing sci frame ", name[index]


  ;make storage 
  chips=make_array(4,2304,4608,/float)


  ;iterate over levels
  for kk=0, 3 do begin
  ;open file 
     chips[kk,*,*]=mrdfits(path+name[index],kk+1,he,/silent,/fscale)
  endfor
  header=headfits(path+name[index],exten=0,/silent)  
  

 ;go for oscan subtraction
  lbtc_oscan, chips, imageout, /silent, fitbias=fitbias
          
 ;remove zero
  if keyword_set(bias) then imageout=temporary(imageout-bias)

                
 ;find corresponding flat field
  flat=where(strcompress(status.filflat,/remove_all) eq strcompress(filter[index],/remove_all),nfl)

  if(nfl gt 0) then begin
     nameflat=side+date+"_"+strcompress(status.filflat[flat],/remove_all)+"flat.fits"
     flatfits=mrdfits(path+nameflat[0],/silent)     
  endif else begin
     ;use archival.. (not yet implemented...)
     ;or keep going with no flats
     splog, 'No flats found for ', filter[index]
     nameflat='No flats!'
     flatfits=imageout-imageout+1.
  endelse
  

  ;get exptime
  exptime=fxpar(header,'EXPTIME')


  ;loop over chips and do stuff (flat,gain,texp)
  for chnum=0, 3 do begin
  ;flat
     tmpimg=reform(imageout[chnum,*,*])
     tmpfla=reform(flatfits[chnum,*,*])
     
     divideflat, tmpimg, tmpfla, minval=0.001, /quiet
  ;set time exp to 1
     tmpimg=tmpimg/exptime
  ;backl
     imageout[chnum,*,*]=tmpimg[*,*]
  ;clean
     undefine, tmpimg, tmpfla
 endfor


  ;apply gain
  gain=replicate(fxpar(header,'GAIN'),4)
  lbtc_applygain, imageout, eleimg, outgain=outgain,  gain_value=gain
  imageout=eleimg


 if keyword_set(skyfit) then begin
  ;remove luminosity pattern
     splog, "Skyfit does not work! :( "
     ;skyfit,imageout,skyimage, XORDER=2, YORDER=2
     ;atv, skyimage, /block
  endif


  ;update header 
 sxaddpar, header, "OBJECT", object[index] 
 sxaddpar, header, "COMMENT", "IMG reduced using LBTC_CCDPROC", after='COMMENT'
 sxaddpar, header, "COMMENT", string("Flat field ",nameflat), after='COMMENT'
 sxaddpar, header, "COMMENT", string("Image in counts/second. Diveded by T=",exptime), after='COMMENT'
 sxaddpar, header, "COMMENT",strjoin(strcompress(string("Gain applied: ",outgain))), after='COMMENT'
 sxaddpar, header, "COMMENT", string("Redux ended ", SYSTIME()),  after='COMMENT'

  ;save science image
  posit=strpos(name[index],".fits")
  substring=strmid(name[index],0,posit) 
  
  ;write chips
  mwrfits, imageout, path+substring+"_redux.fits", header, /create, /silent
    
  ;gzip image
  if keyword_set(gzip) then begin 
     command="gzip "+substring+"_redux.fits"
     spawn, command
  endif
     
  ;wcs
  if keyword_set(wcs) then begin 
     splog, 'wcs not coded yet :('
  endif
 
  ;clean stuff
  undefine, flatfits, imageout, eleimg


END
