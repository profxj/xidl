
pro statarray, image, iaxis, mean=mean, stddev=stddev, $
       sigrej=sigrej, niter=niter

     if N_PARAMS() LT 2 then begin
        print, 'Please give an axis to calculate'
        print, 'Use statarray, statarray, image, 0, mean=mn.... for full array'
        return
     endif


     ndim = (size(image))[0] 
     if ndim LT iaxis then begin
        print, "iaxis is larger than ndim of image"
        return
     endif
     if NOT keyword_set(sigrej) then sigrej=3.0
     if NOT keyword_set(niter) then niter=1

     nbefore = 1L
     nafter  = 1L
     if (iaxis EQ 0) OR (ndim EQ 1) then nsum = n_elements(image)  $
     else begin
       nsum = (size(image))[iaxis]
       for i=1,iaxis-1 do    nbefore = nbefore * (size(image))[i]
       for i=iaxis+1,ndim do nafter  = nafter * (size(image))[i]
     endelse

     ; first iteration
     mn1 = djs_median(image,iaxis)
     
     mnfull = transpose(reform((mn1[*]) # $
                replicate(1,nsum),nbefore,nafter,nsum),[0,2,1])
     stddev = sqrt(djs_median((image - mnfull)^2,iaxis))

     ;; Deal with zero stddev points (rare)
     a = where(stddev EQ 0., na)
     if na NE 0 then begin
         for ii=0L,na-1 do $
           stddev[a[ii]] = mean( stddev[a[ii]-lindgen(5)]+ $
                                 stddev[a[ii]+lindgen(5)] )
     endif
     bad = where(stddev EQ 0., nbad)
     if nbad NE 0 then stop

     for iter=1,niter do begin
       sdfull = transpose(reform((stddev[*]) # $
                replicate(1,nsum),nbefore,nafter,nsum),[0,2,1])
       maskfull = (image GT mnfull - sigrej*sdfull) AND $
               (image LT mnfull + sigrej *sdfull)

       norm = total(maskfull, iaxis)

       mean = total(image*maskfull,iaxis) /(norm + (norm EQ 0)) * (norm GT 0)
       mnfull = transpose(reform((mean[*]) # $
                 replicate(1,nsum),nbefore,nafter,nsum),[0,2,1])
       stddev = sqrt(total((image - mnfull)^2*maskfull,iaxis) / $
                  (norm - (norm NE 1))) * (norm GT 1)
    endfor

return
end
