;+
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
; spec1d_ascii.pro
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
; PURPOSE
; 	A routine to convert spec1d FITS files into ascii tables
; 	containing the flux, ivar, and lambda values in a 1-d
; 	spectrum. 
;
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
; SYNTAX
;	spec1d_ascii, files=files, /horne, /optimal, /boxcar, $
;                     /boxsprof, /nlsky
;	             
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
; INPUTS
;       files = a list of spec1d FITS file names (type = STRING).
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
; KEYWORDS
;       /horne, 
;       /optimal, 
;       /boxcar, 
;       /boxsprof = select the corresponding extraction from the
;                   spec1d FITS files and print to ascii files. If
;                   more than one of these keywords is set, then only
;                   one will be recognized. 
;       /nlsky = select the non-local sky-subtracted extraction if
;                available. 
;
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
; OUTPUTS
;       ascii files written to the current directory with the file
;       name spec1d.xxxx.nnn.ooooooo.dat where xxxx is the mask name,
;       nnn is the slit number and ooooooo is the object name. The
;       ascii data is put in columns with a single initial line to
;       label the columns: lambda, flux, ivar (in that order).
;
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
; PROCEDURES CALLED
;       
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
; EXAMPLES
;	None.
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
; HISTORY
;       Created January 21, 2003 by mcc.
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;-


pro spec1d_ascii, files=files, optimal=optimal, horne=horne, $
                  boxcar=boxcar, boxsprof=boxsprof, nlsky=nlsky

; check the various keyword settings.
  optimal = keyword_set(optimal) ? 1:0
  horne = keyword_set(horne) ? 1:0
  boxcar = keyword_set(boxcar) ? 1:0
  nlsky = keyword_set(nlsky) ? 1:0
  if not(optimal) and not(horne) and $
    not(boxcar) then boxsprof = 1 else boxsprof = 0

  if n_elements(files) eq 0 then files = findfile('spec1d.*.fits*')
  nfiles = n_elements(files)
  for i=0,nfiles-1 do begin
; fill the gap between the red and blue sides of spectrum.
      if optimal then spec1d = fill_gap(files[i], /optimal)
      if horne then spec1d = fill_gap(files[i], /horne)
      if boxcar then spec1d = fill_gap(files[i], /boxcar)
      if boxsprof then spec1d = fill_gap(files[i], /boxsprof)
      if size(spec1d, /tname) eq 'STRUCT' then begin
; construct outut file name.
          len = strlen(files[i])
          filename = strmid(files[i], 0, len-5) + '.dat'
; open file unit for writing.
          openw, unit, filename, /get_lun
; print a comment line at the top of file.
          printf, unit, '#   lambda        flux        ivar'
; print the flux (spec), lambda, and ivar values to the file.
          npts = n_elements(spec1d.spec)
          for j=0,npts-1 do $
            printf, unit, spec1d.lambda[j], '   ', spec1d.spec[j], '   ', $
            spec1d.ivar[j], format='(F10.4, A, F10.4, A, E12.5)'
; close and free the file unit.
          close, unit
          free_lun, unit
      endif else print, 'Unable to create ascii file for ' + files[i]
  endfor

end





