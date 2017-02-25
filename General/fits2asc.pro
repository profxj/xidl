;+ 
; NAME:
; fits2asc
; V1.1
;
; PURPOSE:
;    Converts fits spectrum to ASCII.   The wavelength array is 
;    read from the header.  Particularly useful for VPFIT
;
; CALLING SEQUENCE:
;   
;   fits2asc, file, [error], OUTFIL=
;
; INPUTS:
;   file    - Data Filename [spectrum; data in extension 0]
;   [error] - Error filename
;
; RETURNS:
;
; OUTPUTS:
;   outfil  - ASCII file with wavelength, file, [error] in columns
;
; OPTIONAL KEYWORDS:
;  WVMIN -- Minimum wavelength to print out [default: 0.]
;  WVMAX -- Maximum wavelength to print out [default: 1e6]
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   fits2asc, 'Blah.fits'
;
; PROCEDURES CALLED:
;  x_fitswave
;  writecol
;
; REVISION HISTORY:
;   27-Aug-2001 Written by JXP
;-
;------------------------------------------------------------------------------

pro fits2asc, file, error, OUTFIL=outfil, WVMIN=wvmin, WVMAX=wvmax

  if (N_params() LT 1) then begin 
    print,'Syntax - ' + $
             'fits2asc, file, [error], OUTFIL=, WVMIN=, WVMAX=, [v1.1]'
    return
  endif 

; Keywords

  if not keyword_set( WVMIN ) then wvmin = 0.d
  if not keyword_set( WVMAX ) then wvmax = 100000.d

;  Outfile

;  Check for fits extension  
  len = strlen(file)
  if strmid(file,len-4,4) NE 'fits' then begin
	print, 'File must be a fits file!'
	return
  endif	

  if not keyword_set( OUTFIL ) then outfil = strjoin([strmid(file,0,len-5), $
                                                     '.asc'])

;  Read data

  data = xmrdfits(file, 0, head)
  naxis = sxpar(head, 'NAXIS')
  if naxis NE 1 then begin
      print, 'NAXIS NE 1', naxis
      return
  endif

  if keyword_set( ERROR ) then sig = mrdfits(error)

; Make Wavelength array

  wave = x_fitswave(head)

; Cut on wavelengths as necessary

  goodwv = where( wave LE wvmax AND wave GE wvmin, ngd)
  if ngd EQ 0 then begin
	  print, 'fits2asc: No wavelength values satisfy wvmin, wvmax'
	  return
  endif
  

; Output

  if keyword_set( ERROR ) then writecol, FMT='(F12.6,2F10.4)', outfil, $
    wave[goodwv], data[goodwv], sig[goodwv] $
    else writecol, FMT='(F12.6, F10.4)', outfil, wave[goodwv], data[goodwv]

end

