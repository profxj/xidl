;procedure called by lbtc_makeflats that does the dirty job of
;creating the flats


;str---> structure from lbtc make_plan 
;index ---> array of indexes which contains only flats frame
;bias --> the median bias from makebias if it exists
;FITBIAS --> perform a linear interpolation across prescan/postsca. If
;            not set only median of 0.5*(postcan+prescan) is
;            subtracted from each row
;MINMAX ---> array of min and max value to esclude faint or saturated flats  
;nfilt  ---> from make_flats, the number of images to stack together
;thisfilter  ---> from make_flats, the subindexes of index[...] which
;                 contain fltas for a same filter
;flat_chip----> the final flat [4,2048,4608]


pro lbtc_loopflats,  name, nfilt, thisfilter, index, flat_chip,fitbias=fitbias, bias=bias, minmax=minmax, path=path


  ;cut down to a resonable number
  if(nfilt gt 15) then begin
     splog, 'Too many flats.. Considering only 15!'
     nfilt=15.
  endif


     ;make storage
     flat1=make_array(nfilt,2048,4608,/float)
     flat2=make_array(nfilt,2048,4608,/float)
     flat3=make_array(nfilt,2048,4608,/float)
     flat4=make_array(nfilt,2048,4608,/float)
     
     dump=0

     ;now loop over all the images
     for img=0, nfilt-1 do begin
        
   
     ;open
        splog, "Working on image ", name[index[thisfilter[img]]]  
        chips=make_array(4,2304,4608,/float)
     ;iterate over levels
        for kk=0, 3 do begin
     ;open file 
           chips[kk,*,*]=mrdfits(path+name[index[thisfilter[img]]],kk+1,hea,/silent,/fscale)
        endfor
        

        
    ;check if bright enough or saturated 
        if(djs_median(chips) LT minmax[1] AND djs_median(chips) GT minmax[0]) then begin
          
        ;go for oscan subtraction
           lbtc_oscan, chips, imageout, /silent, fitbias=fitbias
           
        
     ;remove zero
           if keyword_set(BIAS) then imageout=temporary(imageout)-bias
               
   
     ;store each chip after normalisation 
           flat1[img-dump,*,*]=imageout[0,*,*]/djs_median(imageout[0,*,*])
           flat2[img-dump,*,*]=imageout[1,*,*]/djs_median(imageout[1,*,*])
           flat3[img-dump,*,*]=imageout[2,*,*]/djs_median(imageout[2,*,*])
           flat4[img-dump,*,*]=imageout[3,*,*]/djs_median(imageout[3,*,*])
           
        endif else begin
         ;take out spot for saturated image
           splog, "Dump flat ",  name[index[thisfilter[img]]]
           dump=dump+1
           flat1=temporary(reform(flat1[0:n_elements(flat1[*,0,0])-2,*,*]))
           flat2=temporary(reform(flat2[0:n_elements(flat2[*,0,0])-2,*,*]))
           flat3=temporary(reform(flat3[0:n_elements(flat3[*,0,0])-2,*,*]))
           flat4=temporary(reform(flat4[0:n_elements(flat4[*,0,0])-2,*,*]))
        endelse
           
        ;clean stuff
        undefine, chips, imageout

     endfor

    ;here stack them with a simple median
     flat_chip=make_array(4,2048,4608,/double)
     tmp=djs_median(flat1,1)
     flat_chip[0,*,*]=tmp[*,*]
     undefine, tmp
     tmp=djs_median(flat2,1)
     flat_chip[1,*,*]=tmp[*,*]
     undefine, tmp
     tmp=djs_median(flat3,1)
     flat_chip[2,*,*]=tmp[*,*]
     undefine, tmp
     tmp=djs_median(flat4,1)
     flat_chip[3,*,*]=tmp[*,*]

     
     ;clean stuff
     undefine, flat1, flat2, flat3, flat4, tmp


END
