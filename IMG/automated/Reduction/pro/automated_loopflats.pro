;+
;
; NAME:
;   automated_loopflats
;
; PURPOSE:
;   Creates flats for every filter
;    
; CALLING SEQUENCE:
;
; INPUTS:
;   
; str        - the structure with relevant information
; index      - the index where the bias is
; nfilt      - number of images to stack together
; flat_chip  - the final flat
;
; OPTIONAL INPUTS:
;
;  bias   - the median bias from makebias if it exists
;  MINMAX - array of min and max value to esclude faint or saturated flat
;  STATUS - the status structre from ccdproc 
;  PLAN   - plan structure
;
; OUTPUTS: 
;
; Make the actual flats
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



pro automated_loopflats, str, nfilt, thisfilter, index, flat_chip, $
                         bias=bias, minmax=minmax, status=status
  
  ;;make storage
  flat1=make_array(nfilt,2048,2048,/float)
  dump=0
  
  ;;now loop over all the images
  for img=0, nfilt-1 do begin
      
      
      ;;open
      splog, "Working on image ", str.NAME[index[thisfilter[img]]]  
      chips=mrdfits(strtrim(str.NAME[index[thisfilter[img]]],2),0,hea,/silent,/fscale)
      
      ;;check if bright enough or saturated 
      if(djs_median(chips) LT minmax[1] AND djs_median(chips) GT minmax[0]) then begin
          
          ;;go for oscan subtraction
          automated_oscan, chips, imageout
          
          ;;remove zero
          if keyword_set(BIAS) then imageout=imageout-bias
          
          ;;remove dark current
          if(status.dark[0] NE -9999) then begin
              darkval=status.dark[0]*str.EXTIME[index[thisfilter[img]]]+status.dark[1]
              splog, "Subtract dark current ", darkval
              imageout=imageout-darkval
          endif
          
          ;;store chip after normalisation for flats 
          flat1[img-dump,*,*]=imageout[*,*]/djs_median(imageout) ;normal action
          
      endif else begin
          ;;take out spot for saturated image
          splog, "Dump flat ",  str.NAME[index[thisfilter[img]]]
          ;;if bad flat, then keep same index "img-dump"
          dump=dump+1
          flat1=reform(flat1[0:N_ELEMENTS(flat1[*,0,0])-2,*,*])
      endelse       
  endfor
  
  ;;If there are no good flats then find last good ones and try this over again
  if dump eq nfilt then begin
      splog, "No working flats. Extracting archive."
      ;;Read in archive of superflat
      flat_chip = mrdfits(getenv("AUTOM_DIR")+'/archive/archivedflat'+str.FILTER[index[thisfilter[0]]]+'.fits') 
  endif else begin
      ;;here stack them with a simple median
      flat_chip=djs_median(flat1,1)
      ;;write to the archive for good flats 
      mwrfits, flat_chip, getenv("AUTOM_DIR")+"/archive/archivedflat"+str.FILTER[index[thisfilter[0]]]+".fits", /create
  endelse
     
     

END
