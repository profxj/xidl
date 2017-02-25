
;procedure that take the 4 chips as input and subtract the overscan bias
;return a [4,2048,4608] data image bias subtracted 

;FITBIAS --> perform a linear interpolation across prescan/postsca. If
;            not set only median of 0.5*(postcan+prescan) is
;            subtracted from each row
;    
;   




pro lbtc_oscan, chips, dataimage, fitbias=fitbias, silent=silent
 
  dataimage=make_array(4,2048,4608,/float)
  
;run over all the chips
  
  FOR ch=0, 3 DO BEGIN
     
     data=reform(chips[ch,50:2097,0:4607])
     
     ;get bias info
     prescan=reform(chips[ch,0:49,0:4607])
     ;cut first pixels in the postscan
     overscan=reform(chips[ch,2150:2303,0:4607])
     diff=median(djs_median(overscan,1) - djs_median(prescan,1))
     
     IF(diff GT 5) THEN $
        splog, "Prescan postscan differ by", diff ,"counts. Try FITBIAS"

     
     if keyword_set(fitbias)  then begin
        splog, "WARNING: It looks like there is a positive slope in the bias"
    ;median array
        over=djs_median(overscan,1)
        pres=djs_median(prescan,1)
    ;interoplate
        m=(over-pres)/2048.
        x=make_array(2048,/float,/index)
        x=rebin(x,2048,4608)
        m=rebin(transpose(m),2048,4608)
        pres=rebin(transpose(pres),2048,4608)
        bias=m*x+pres
     endif else begin
  
        ;equal for both side
        ;get mean bias
        bia=0.5*(djs_median(overscan,1)+djs_median(prescan,1))
        ;smoot it (25 pix window, 3 deg poly)
        ;plot, bia, /ynozero
        smoo_pol=SAVGOL(25,25,0,3) 
        bia_sm=convol(bia,smoo_pol,/EDGE_TRUNCATE)
        ;oplot, bia_sm, color=150
        ;stop
        bias=rebin(transpose(bia_sm),2048,4608)
       
     endelse
        
     data=data-bias
     
  
     if ~keyword_set(silent) then splog, "Done chip ", ch+1
     dataimage[ch,0:2047,0:4607]=data[0:2047,0:4607]
     
  endfor
  	
  
  ;clean stuff
  undefine, data, prescan, overscan, diff, bias, bia
  undefine, pres, m, x, bia_sm

  
  
end
