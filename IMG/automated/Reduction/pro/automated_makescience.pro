;+
; NAME:
;   automated_makescience
;
; PURPOSE:
; Procedure that prepares the final science frames
;
;    
; CALLING SEQUENCE:
;
; INPUTS:
;   
; 
; OPTIONAL INPUTS:
;  
;
;
; OUTPUTS: 
;
; Creates a bunch of reduced science frames  
;
; OPTIONAL OUTPUTS:
;   
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;    05-Jun-2012  Written and revised by MF
;
;-
;------------------------------------------------------------------------------
;; 


pro automated_makescience, str, index, status=status, plan=plan, bias=biasmedian, skyfit=skyfit, gzip=gzip

  splog, "Processing sci frame ", str.name[index]

  ;;open file 
  chips=mrdfits(strtrim(str.name[index],2),0,header,/silent,/fscale)
  
  ;;go for oscan subtraction
  automated_oscan, chips, imageout
      
  ;;remove zero if set
  if keyword_set(BIAS) then imageout=imageout-bias
   
  ;;remove dark current
  if(status.dark[0] NE -9999) then begin
     darkval=status.dark[0]*str.extime[index]+status.dark[1]
     splog, "Subtract dark current ", darkval
     imageout=imageout-darkval
  endif

  
  
  ;;find corresponding flat field
  flat=where(strcompress(status.filflat,/remove_all) eq strcompress(str.filter[index],/remove_all))
  nameflat=strcompress(status.filflat[flat],/remove_all)+plan+"flat.fits"
  flatfits=mrdfits(nameflat[0],/silent)     
  

  ;;do stuff over the chip
  ;;flat
  divideflat, imageout, flatfits, minval=0.001
  ;;set time exp to 1
  imageout=imageout/str.EXTIME[index]
  ;;apply gain
  imageout=imageout*str.gain[index]


  if keyword_set(skyfit) then begin
      ;;remove luminosity pattern
      splog, "Skyfit does not work! :( "
  endif

  ;;update WCS information accounting for datasec (needed?)
  crpix1=fxpar(header,"CRPIX1")
  crpix1=crpix1-51.
  SXADDPAR, header, "CRPIX1", crpix1


  ;;update header 
  SXADDPAR, header, "HISTORY", "IMG reduced using SUPL_CCDPROC", AFTER='COMMENT'
  SXADDPAR, header, "HISTORY", string("Flat field ",nameflat), AFTER='COMMENT'
  SXADDPAR, header, "HISTORY", string("Image in counts/second"), AFTER='COMMENT'
  SXADDPAR, header, "HISTORY", string("Gain applied: ",str.GAIN[index]), AFTER='COMMENT'
  SXADDPAR, header, "HISTORY", string("Redux ended ", SYSTIME()), AFTER='COMMENT'

  ;;update geometry of images
  sxdelpar, header, 'BIASSEC'
  SXADDPAR, header, 'TRIMSEC', '[1:2048,1:2048]'
  SXADDPAR, header, 'DATASEC', '[1:2048,1:2048]'


  ;;save science image
  posit=STRPOS(str.NAME[index],".fit")
  substring=STRMID(str.NAME[index],0,posit) 
  
  ;;write chips
  mwrfits, imageout, substring+"_redux.fits", header, /CREATE, /SILENT
 

  if keyword_set(gzip) then begin
      ;;gzip image
     command="gzip "+substring+"_redux.fits"
     spawn, command
     splog, "Done with "+substring+"_redux.fits.gz"
  endif else splog, "Done with "+substring+"_redux.fits"


end
