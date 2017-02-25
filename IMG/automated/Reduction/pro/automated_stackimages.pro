;+
; NAME:
;   automated_stackimages
;
; PURPOSE:
; Procedure that alignes different frames and produced a final stacked
; image with wcs
;
;    
; CALLING SEQUENCE:
;
; INPUTS:
;   
; 
; 
; OPTIONAL INPUTS:
;  
; 
; OUTPUTS: 
;
; Creates a bunch of stacked science frames with wcs  
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


pro automated_stackimages, plan, gzip=gzip, wcs=wcs, date=date, badseeing=badseeing
  

  ;;open the plan structure 
  str=mrdfits(plan,1,/silent)
    
  ;;find science, bias and flat
  ;;remember you changed ima to IMA
  science=where(str.type EQ 'IMA', nsci)

  ;;open status 
  status=mrdfits(plan+"status.str",1,/silent)


  ;;make associations (for each object find frames, then stack according
  ;;to filters)

  ;;find targets
  for obj=0, nsci-1 do begin
     if (obj eq 0) then objlist=strtrim(str.object[science[obj]],2) 
     if (obj gt 0) then begin
        find=where(objlist EQ strtrim(str.object[science[obj]],2),nfound)
        ;;if new, append
        if(nfound EQ 0) then objlist=[objlist,strtrim(str.object[science[obj]],2)]
     endif
  endfor
  
  ntarg=n_elements(objlist)

  splog, "Found ", ntarg ," targets: ", objlist
  
  ;;if no targets found, then done
  if(ntarg eq 0) then return
  
  ;;make final directory
  spawn, 'mkdir -p Sci'

  ;;loop over targets; make associations with filters
  for targ=0, ntarg-1 do begin
      
      splog, "Working with ", objlist[targ]
      
      ;;find single fits
      ind=where(strtrim(str.object,2) eq strtrim(objlist[targ],2),nnn)
      if(nnn gt 0) then begin
          images=str.name[ind]
          filter=str.filter[ind]
      endif
      splog, "Found fits: ",  images
      
      ;;indetify filter 
      for fil=0, n_elements(images)-1 do begin
          if (fil eq 0) then currfilt=strtrim(filter[fil],2) 
          if (fil gt 0) then begin
              find=where(currfilt EQ filter[fil] ,nfound)
              ;;if new, append
              if(nfound EQ 0) then  currfilt=[currfilt,strtrim(filter[fil],2)]
          endif
      endfor
      
      splog, "For this target, found filters: ",  currfilt
      
      nfilt=n_elements(currfilt)
      
      ;;now loop over filters and stack
      for fil=0, nfilt-1 do begin
          splog, "Considering filter ", currfilt[fil]
          
          ;;make save name
          savename='sci_'+strtrim(objlist[targ],2)+'_'+strtrim(currfilt[fil],2)+'.fits'
          

          ;;open ps check file
          x_psopen, "checkstack_"+strtrim(objlist[targ],2)+'_'+strtrim(currfilt[fil],2)+'.ps'

          ;;check if already done
          file=file_info('./Sci/'+savename)
          if(file.exists eq 1) then splog, savename+' already processed!' else begin
              
              ind=where(strtrim(filter,2) eq currfilt[fil])
              stackname=images[ind]
              
              splog, "Stacking ", stackname
              nstack=n_elements(stackname)
              
              ;;make storage
              chips=make_array(nstack,2048,2048)
              nameredux=stackname
              

              ;;open fits
              for st=0, nstack-1 do begin
                  ;;make name
                  posit=STRPOS(stackname[st],".fit")
                  substring=STRMID(stackname[st],0,posit) 
                  openfile=substring+"_redux.fits"
                  nameredux[st]=openfile
                  if(st eq 0) then fits=mrdfits(openfile,0,hea,/silent)$
                  else fits=mrdfits(openfile,0,/silent)
                  chips[st,*,*]=fits[*,*]
              endfor
              
              ;;align frames
              splog, 'Align frames...'
              if keyword_set(badseeing) then fwhm=replicate(3.,nstack) else fwhm=replicate(1.5,nstack)
     
              automated_offset, nameredux, outx, outy
              
              ;;shift
              for st=0, nstack-1 do begin
                  fits=reform(chips[st,*,*])
                  shiftimg=shiftf(fits,outx[st],outy[st])
                  chips[st,*,*]=shiftimg
              endfor
              
              ;;stacking with a median to reject cosmic rays (simple, no scailing)
              splog, 'Combine... '
              combimg=djs_median(chips,1)
              
              ;;save science frame
              sxaddpar, hea, "HISTORY", string("Combined multiple frames "), AFTER='COMMENT'
              sxaddpar, hea, "NCOMBINE", nstack
              sxaddpar, hea, "HISTORY", strjoin(strcompress(string("Shift applied x: ",outx))),$
                      AFTER='COMMENT'
              sxaddpar, hea, "HISTORY", strjoin(strcompress(string("Shift applied y: ",outy))),$
                      AFTER='COMMENT'
           
              mwrfits, combimg, savename, hea, /create
              
              ;;fit new wcs
              if keyword_set(wcs) then automated_fitwcs, savename
              
              ;;gzip
              if keyword_set(gzip)then  spawn, "gzip  "+savename
              
              ;;mv final file
              spawn, 'mv '+savename+'* ./Sci/. '
              
              ;;close output
              x_psclose

          endelse
      endfor
  endfor
  
  
;;clean
spawn, "rm -f default.sex"
spawn, "rm -f default.param"
spawn, "rm -f default.nnw"
spawn, "rm -f default.conv"

 
end









