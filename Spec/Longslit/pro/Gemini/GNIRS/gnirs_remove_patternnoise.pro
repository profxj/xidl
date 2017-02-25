pro gnirs_remove_patternnoise, image

   sz = size(image)
   nx = sz[1]
   ny = sz[2]

   result = 0 * image

   pattern = fltarr(nx/2)

; Upper left quadrant

   for j=0, nx/2-1 do begin
      pattern[j] = median(image[j, ny/2:ny-1])
   endfor
   model = gnirs_fit_patternnoise(pattern[0:200], nx/2)
;   plot, pattern
;   oplot, model, color=getcolor('green')
;   oplot, pattern-model, color=getcolor('red')
   for j=0, nx/2-1 do begin
       result[j,ny/2:ny-1] = image[j,ny/2:ny-1] - model[j]
   endfor
;   stop

; Lower left quadrant

   for j=0, nx/2-1 do begin
       pattern[j] = median(image[j,0:ny/2-1])
   endfor
   model = gnirs_fit_patternnoise(pattern[0:200], nx/2)
;   plot, pattern
;   oplot, model, color=getcolor('green')
;   oplot, pattern-model, color=getcolor('red')
   for j = 0, nx/2-1 do begin
      result[j, 0:ny/2-1] = image[j, 0:ny/2-1] - model[j]
   endfor
;   stop

; Upper right quadrant

   for j=0, nx/2-1 do begin
       pattern[j] = median(image[j+nx/2,ny/2:ny-1])
    endfor
   model = gnirs_fit_patternnoise(pattern[nx/4:nx/2-1], nx/2)
;   plot, pattern
;   oplot, model, color=getcolor('green')
;   oplot, pattern-model, color=getcolor('red')
   for j=0, nx/2-1 do begin
       result[j+nx/2,ny/2:ny-1] = image[j+nx/2,ny/2:ny-1] - model[j]
    endfor
   
;   stop

; Lower right quadrant
   for j=0, nx/2-1 do begin
       pattern[j] = median(image[j+nx/2,0:nx/2-1])
   endfor
   model = gnirs_fit_patternnoise(pattern[nx/4:nx/2-1], nx/2)
;   plot, pattern
;   oplot, model, color=getcolor('green')
;   oplot, pattern-model, color=getcolor('red')
   for j=0, nx/2-1 do begin
       result[j+nx/2,0:ny/2-1] = image[j+nx/2,0:ny/2-1] - model[j]
   endfor
   image = result
   return

end
