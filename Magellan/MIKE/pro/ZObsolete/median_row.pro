
;
;  Returns a median smooth image with width provided along
;   each individual row in image.  
;
;  Uses buffering to avoid looping, but works efficiently
;   only when ncol >> width
;
;

function median_row, image, width, reflect=reflect

      ncol = (size(image))[1]
      nrow = 1 
      if (size(image))[0] EQ 2 then nrow = (size(image))[2] 

      hw = (long(width) /2) > 1
      if keyword_set(reflect) then begin
         prebuffer = reverse(image[1:hw,*])
         postbuffer = reverse(image[ncol-hw-1:ncol-2,*])
      endif else begin
         prebuffer = replicate(1,hw) # image[0,*]
         postbuffer = replicate(1,hw) # image[ncol-1,*]
      endelse

      work_array = ([prebuffer, image, postbuffer])[*]

      result = median(work_array,width)
      final = (reform(result,ncol+hw*2, nrow))[hw:hw+ncol-1,*]

      return, final
end
