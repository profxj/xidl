function BRI_getID, bris, xpix, ypix

; parse_phot -- Reads in BRI data to a structure

  if (N_params() LT 3) then begin 
    print,'Syntax - ' + $
             'BRI_getID(structure, xpix, ypix)'
    return, -1
  endif 

;

  nobj = n_elements(bris)
  j = 0
  min = 99.
  while(j LT nobj-1 AND min GT 2.0) do begin
      min = sqrt( (bris[j].xpix-xpix)^2 + (bris[j].ypix-ypix)^2 )
      j=j+1
  endwhile

  if(j LT nobj-1 AND min LT 2.0) then begin
      print, 'Successful:', j-1
      return, j-1 
  endif else begin
      print, 'No object at ', xpix, ypix
      return, -1
  endelse
end
