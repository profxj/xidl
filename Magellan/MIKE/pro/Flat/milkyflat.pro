
;
;  TODO:  flativar is not constructed from the medium flat
;      :  If fitsout is used, then read back in and construct flativar
;      :  Use default flat mask file depending on hdr/filename
;           This is need to mask the extra bad holes on the red/blue sides
;
;  INPUTS:
;      milkyflat   - filenames of raw MIKE files which have been observed with 
;                        diffuser in beam, producing a "milky" flat
;
;  COMMENTS:  This is a different version, which used 2d medians and two
;             passes of b-spline to attempt to get a smooth model for the
;             milky image.  It has only really been tested on one set of
;             red and blue milky daylight sky flats, so it could easily
;             break or fail for different milky images.

 
function milky_med, milkyfile, flatstddev=flatstddev, GAIN=gain, RN=rn

     if NOT keyword_set(rn) then rn = 4.1
     if NOT keyword_set(gain) then gain = 1.0
     if NOT keyword_set(fitsout) then fitsout=0

     nfiles = n_elements(milkyfile)

     if nfiles GT 1 then begin
         print, 'Doing 1st of ', nfiles
         first = milkyflat(milkyfile[0]) 
         nx = (size(first))[1]
         ny = (size(first))[2]
         final = fltarr(nx,ny, nfiles)
         final[*,*,0] = first
         for i=1,nfiles -1 do begin
           print, 'Doing ', i+1, (i EQ 1 OR i EQ 2) ? 'rd of' : 'th of', nfiles
           final[*,*,i] = milkyflat(milkyfile[i])
         endfor
         statarray, final, 3, mean=flat, stddev=flatstddev
         return, flat
     endif


     ; Open file
     milky = xmrdfits(milkyfile, /silent)
     milky = milky*gain
     milkyivar =  1./(abs(milky)*gain + rn^2)

     ; Transpose for quicker execution
 
     image  = milky
     invvar = milkyivar

;     image  = transpose(milky)
;     invvar = transpose(milkyivar)

     ncol = (size(image))[1]
     nrow = (size(image))[2]

     width = 2*(ncol/256L) + 1L
     print, ncol, width
     prebuffer = image[reverse(findgen(width/2))+1,*]
     postbuffer = image[ncol - 2 - reverse(findgen(width/2)),*]
     tempmed = median([prebuffer,image,postbuffer], width)
     medimage = tempmed[width/2:width/2+ncol-1,*]

     imagefit = image * 0.0
     workfit = image/(medimage + (medimage EQ 0))
     invvarfit = invvar * medimage^2

     print, 'First b-spline run'
     x = findgen(nrow) & $
     for icol= 0,ncol-1 do begin & $
       sset = bspline_iterfit(x, workfit[icol,*], invvar=invvarfit[icol,*], $
                   lower=3, upper=7, everyn=nrow/25, yfit=tempfit, /silent, $
                   /groupbadpix, maxrej=3) & $
       imagefit[icol,*] = tempfit  &  $
     endfor

     ratio = workfit/(imagefit + (imagefit EQ 0))
     bsplinemask = smooth(ratio,width) GT (1 - 0.045/width)
     ratioivar = invvarfit * imagefit^2 * (imagefit GT 0.5) * bsplinemask
     finalfit = ratio*0.0

     print, 'Second b-spline run'
     x = findgen(nrow) & $
     for icol= 0,ncol-1 do begin & $
       sset = bspline_iterfit(x, ratio[icol,*], invvar=ratioivar[icol,*], $
                   lower=3, upper=7, everyn=nrow/25, yfit=tempfit, /silent, $
                   /groupbadpix, maxrej=3) & $
       finalfit[icol,*] = tempfit  &  $
     endfor


     flat = ratio/(finalfit + (finalfit EQ 0)) * (finalfit GT 0)
 
     return, flat

end

