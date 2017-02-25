
;procedure that make a median bias from different bias frames


;name---> list of names lbtc_ccdproc
;index ---> array of indexes which contains only bias frame
;fitbias --> perform a linear interpolation across prescan/postsca. If
;            not set only median of 0.5*(postcan+prescan) is
;            subtracted from each row
;biasmedian --> in output the median bias
   


PRO lbtc_makebias, name, biasmedian, index, date=date, side=side, path=path,$
                   fitbias=fitbias
 

  
  if(n_elements(index) gt 10) then begin
     index=index[0:10]
     splog, "Considering only 11 bias frames!"
  endif


  
;make storage for 4 chips
  bias1=make_array(n_elements(index),2048,4608,/float)
  bias2=make_array(n_elements(index),2048,4608,/float)
  bias3=make_array(n_elements(index),2048,4608,/float)
  bias4=make_array(n_elements(index),2048,4608,/float)
  

  for pos=0, n_elements(index)-1 do begin
;open files 
     chips=make_array(4,2304,4608,/float)
     splog, "working on  bias ", name[index[pos]]  

;iterate over levels
     for kk=0, 3 do begin
;open file 
        chips[kk,0:2303,0:4607]=mrdfits(path+name[index[pos]],kk+1,hea,/silent,/fscale)
     endfor
     
    

;go for oscan subtraction
     lbtc_oscan, chips, biasimage, /silent, fitbias=fitbias
    
     ;clean memory
     undefine, chips
      
;store each chip
     bias1[pos,0:2047,0:4607]=biasimage[0,0:2047,0:4607]
     bias2[pos,0:2047,0:4607]=biasimage[1,0:2047,0:4607]
     bias3[pos,0:2047,0:4607]=biasimage[2,0:2047,0:4607]
     bias4[pos,0:2047,0:4607]=biasimage[3,0:2047,0:4607]
     
  endfor
  
;combine the final bias
  splog, "creating median bias... "  

  biasmedian=make_array(4,2048,4608,/float)
  
  biasmedian[0,0:2047,0:4607]=djs_median(bias1,1) 
  biasmedian[1,0:2047,0:4607]=djs_median(bias2,1) 
  biasmedian[2,0:2047,0:4607]=djs_median(bias3,1) 
  biasmedian[3,0:2047,0:4607]=djs_median(bias4,1) 

  ;save median bias
  mwrfits, biasmedian, path+side+date+"medbias.fits",  /create

  splog, "all done with median bias!"

    
  ;clean stuff
  undefine, bias1, bias2, bias3, bias4


end
