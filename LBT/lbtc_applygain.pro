



pro lbtc_applygain, adu_image, electron_image, outgain=outgain,  gain_value=gain_value

  ;tweak the gain to level in each
  ;chip. Consider central chip (number 2,index1) the reference.
  ;--------------------------------------------------
 


 med_amp=fltarr(4)
  
  ;get the scale you need
  for amp=0, 3 do begin
     ;trim some [50:2000,50:4550] 
     med_amp[amp]=djs_median(adu_image[amp,50:2000,50:4550]*gain_value[amp])
  endfor


  ;correct
  splog, 'Tweak gain... Corr_12: ', med_amp[0]/med_amp[1]
  splog, 'Tweak gain... Corr_32: ', med_amp[2]/med_amp[1]
  splog, 'Tweak gain... Corr_42: ', med_amp[3]/med_amp[1]

  outgain=gain_value
  outgain[0]=gain_value[0]*med_amp[1]/med_amp[0]
  outgain[2]=gain_value[2]*med_amp[1]/med_amp[2]
  outgain[3]=gain_value[3]*med_amp[1]/med_amp[3]

  ;check if resonable 
  nogood=where(outgain/outgain[1] gt 1.2 or  outgain/outgain[1] lt 0.8, ncheck)
  if(ncheck gt 0) then begin
     splog, 'Variation for chips ', nogood+1, ' bigger than 20%!.'
     splog, 'Restore gain ', outgain[1]
     outgain[nogood]=outgain[1]
  endif


  ;apply gain
  electron_image=adu_image
  for amp=0, 3 do begin
     electron_image[amp,*,*]=adu_image[amp,*,*]*outgain[amp]
  endfor
  
end
