

;get a  median bias and find bad column

pro lbtc_mapbadcolum, list

  ;final mask
  chip_mask=fltarr(4,2048,4608)


  ;consider one chip at time
  for chip=0, 3 do begin
     splog, 'Look at chip ', chip+1
     
  ;open  fits
     chip_image=mrdfits(list)
     
  ;get the chip    
     median_all=reform(chip_image[chip,*,*])
  
  ;sigma_clip to make the mask
     djs_iterstat, median_all, mask=mask, sigrej=4.5

     ;unmask first 50 pixel on blue side
     ;mask[0:50,*]=mask[0:50,*]-mask[0:50,*]+1.

     atv, median_all, /block
     atv, mask, /block
     
  ;store
     chip_mask[chip,*,*]=mask[*,*]
     undefine, fits, median_all, column, image, mask
 
  endfor  ;loop over chip


  ;write entire mask
   time=string(SYSTIME(/julian))
   allname=strcompress("badpix"+time+"_.fits",/remove_all)
   mwrfits, chip_mask, allname,  /CREATE, /SILENT
 
  splog, 'Output written in ', allname


end
