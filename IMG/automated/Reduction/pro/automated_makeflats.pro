;+
;
; NAME:
;   automated_makeflats
;
; PURPOSE:
;   Creates flats for every filter
;    
; CALLING SEQUENCE:
;
; INPUTS:
;   
; str     - the structure with relevant information
; index   - the index where the bias is
; 
; OPTIONAL INPUTS:
;
;  bias   - the median bias from makebias if it exists
;  MINMAX - array of min and max value to esclude faint or saturated flat
;  STATUS - the status structre from ccdproc 
;  PLAN   - plan structure
;  
;
; OUTPUTS: 
;
; Find the dark current value as function of time
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


pro automated_makeflats, str, index,  bias=bias, status=status, $
                         plan=plan, minmax=minmax
 

  ;;if no flats, archive
  if(index[0] eq -1) then begin
      
      ;;find group of common filters 
      for i=0, n_elements(str.filter)-1 do begin
          if(i eq 0) then flat=strcompress(str.filter[i],/remove_all)
          if(i gt 0) then begin
              ;;check if new filter
              res=where(flat EQ strcompress(str.FILTER[i],/remove_all),nummatch)
              ;;if new, append
              if(nummatch EQ 0) then flat=[flat,strcompress(str.FILTER[i],/remove_all)]
          endif
      endfor
      
      for flfra=0, n_elements(flat)-1 do begin
          finalflat = mrdfits(getenv("AUTOM_DIR")+'/archive/archivedflat'+flat[flfra]+'.fits') 
          ;;save final output
          mwrfits, finalflat, flat[flfra]+plan+"flat.fits", /create
          ;;update status structure
          status.filflat[flfra]=flat[flfra]
          mwrfits, status, plan+"status.str", /create
      endfor

  endif else begin
      
      ;;else, loop over flats to find group of common filters 
      for i=0, n_elements(index)-1 do begin
          if(i eq 0) then flat=strcompress(str.filter[index[i]],/remove_all)
          if(i gt 0) then begin
              ;;check if new filter
              res=where(flat EQ strcompress(str.FILTER[index[i]],/remove_all),nummatch)
              ;;if new, append
              if(nummatch EQ 0) then flat=[flat,strcompress(str.FILTER[index[i]],/remove_all)]
          endif
      endfor
      
      splog, "Filters found: ", flat 
      
      
      ;;make flats
      for flfra=0, n_elements(flat)-1 do begin
          ;;check if flat done
          sss=where(strcompress(status.filflat,/remove_all) eq flat[flfra],numbfl)
          
          if(numbfl gt 0) then splog, "Found flat ", flat[flfra] else begin
              splog, "Working on filter ", flat[flfra]
              
              ;;find images for this filter
              thisfilter=where(strcompress(str.FILTER[index],/remove_all) EQ flat[flfra],nfilt)
              
              ;;check if everything is ok (nfilt MUST be > 0)
              if(nfilt EQ 0) then begin
                  splog, "I cannot find images with this filter: ", flat[flfra], " Stop..."
                  stop
              endif
              
              ;;call loop flats to make final image 
              automated_loopflats, str, nfilt, thisfilter, index, finalflat,$
                      BIAS=bias, MINMAX=minmax, status=status
              
              ;;save final output
              mwrfits, finalflat, flat[flfra]+plan+"flat.fits", /create
              
              ;;update status structure
              status.filflat[flfra]=flat[flfra]
              mwrfits, status, plan+"status.str", /create
          endelse
      endfor
      
  endelse
  
  splog, "All done with flats!!"
  
end
