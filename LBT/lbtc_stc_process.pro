
;procedure that does the job of stacking and calls child processes



pro lbtc_stc_process, img_name, img_obje, img_filt, numimg, $
                      path=path, currside=currside, $
                      currchip=currchip, noskysub=noskysub, creject=creject,$
                      manshift=manshift 
  

;find combination of name
  for obj=0, numimg-1 do begin
     if (obj eq 0) then objlist=strtrim(img_obje[obj],2) 
     if (obj gt 0) then begin
        find=where(objlist eq strtrim(img_obje[obj],2),nfound)
        ;if new, append
        if(nfound EQ 0) then objlist=[objlist,strtrim(img_obje[obj],2)]
     endif
  endfor
  
  ntarg=n_elements(objlist)

  splog, "Found ", ntarg ," targets: ", objlist
 
  if(ntarg eq 0) then begin
     splog, 'No targets found!'
     return
  endif



    
;loop over targets; make associations with filters

  for targ=0, ntarg-1 do begin
 
     
     ;check if already done
     file=file_info('alldone_'+objlist[targ]+currside+string(currchip)+'.utility')
     if(file.exists eq 1) then splog, objlist[targ]+' already done!' else begin
        
        splog, "Working with ", objlist[targ]

     ;find single fits
        ind=where(strtrim(img_obje,2) eq strtrim(objlist[targ],2),nnn)
        if(nnn gt 0) then begin
           images=img_name[ind]
           filter=img_filt[ind]
        endif
        ;splog, "Found fits: ",  images
        
     ;indetify filter 
        for fil=0, n_elements(images)-1 do begin
           if (fil eq 0) then currfilt=strtrim(filter[fil],2) 
           if (fil gt 0) then begin
              find=where(currfilt EQ strtrim(filter[fil],2),nfound)
              ;if new, append
              if(nfound EQ 0) then  currfilt=[currfilt,strtrim(filter[fil],2)]
           endif
        endfor
     
        splog, "For this target, found filters: ",  currfilt
   
        nfilt=n_elements(currfilt)


     ;now loop over filters and stack
        for fil=0, nfilt-1 do begin
           splog, "Considering filter ", currfilt[fil]
           ind=where(strtrim(filter,2) eq currfilt[fil])
           stackname=images[ind]
           
           splog, "Stacking ", stackname
           nstack=n_elements(stackname)
           
           ;get the header 
           header=headfits(stackname[0])
                      
           sxaddpar, header, "NCOMBINE", nstack, before='DATASUM'
              
           ;go to the actual stack procedure
           lbtc_getcombined, stackname, outimg, outmask, header,$
                             medianout=medianout, path=path, noskysub=noskysub, $
                             creject=creject,$
                             currside=currside, currchip=currchip,$
                             manshift=manshift  

                      
           ;write the final file and the total mask
           save_sci_name='./Sci/'+strcompress(objlist[targ]+'_'+currfilt[fil]+'_sci'+string(currchip)+'.fits',/remove_all)
           save_mask_name='./Sci/'+strcompress(objlist[targ]+'_'+currfilt[fil]+'_mask'+string(currchip)+'.fits',/remove_all)
           save_median='./Sci/'+strcompress(objlist[targ]+'_'+currfilt[fil]+'_median'+string(currchip)+'.fits',/remove_all)
           
           sxaddpar, header, "COMMENT", 'Img combined with LBTC_STACKIMAGES', $
                     before='DATASUM'
           sxaddpar, header, "COMMENT", 'Img generated on '+systime(),$
                     before='DATASUM'

           ;delete NAXIS3
           sxdelpar, header, 'NAXIS3'

           ;fix nan pixels in outimg (already did in getcombined)
           nan_pix=where(finite(outimg,/nan),ngd)
           if(ngd gt 0) then outimg[nan_pix]=0.

             
           mwrfits, outimg, save_sci_name, header, /create
           mwrfits, outmask*100., save_mask_name, /create
           mwrfits, medianout, save_median, header, /create
           
           ;clean 
           undefine, outimg, outmask, medianout
           
        endfor                  ;filters
     
           
     endelse                    ;if new process

     
     ;flag object done
     namedone=strjoin('alldone_'+objlist[targ]+currside+string(currchip)+'.utility')
     spawn, 'touch  '+namedone
       

  endfor                        ;loop on targets
  
   ;clean
  file_delete, 'default.conv'
  file_delete, 'default.sex'
  file_delete, 'default.param'
  file_delete, 'default.nnw'


end
