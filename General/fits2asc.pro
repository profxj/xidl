;+ 
; NAME:
; fits2asc
;
; PURPOSE:
;    Converts fits to ASCII.  Useful for VPFIT
;
; CALLING SEQUENCE:
;   
;   fits2asc, file, [error], OUTFIL=
;
; INPUTS:
;   file    - Data Filename
;   [error] - Error filename
;
; RETURNS:
;
; OUTPUTS:
;   outfil  - ASCII file with wavelength, file, [error]
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   fits2asc, 'Blah.fits'
;
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   27-Aug-2001 Written by JXP
;-
;------------------------------------------------------------------------------

pro fits2asc, file, error, OUTFIL=outfil, WVMIN=wvmin, WVMAX=wvmax

  if (N_params() LT 1) then begin 
    print,'Syntax - ' + $
             'fits2asc, file, [error], OUTFIL=, WVMIN=, WVMAX='
    return
  endif 

; Keywords

  if not keyword_set( WVMIN ) then wvmin = 0.0
  if not keyword_set( WVMAX ) then wvmax = 100000.0

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

  data = mrdfits(file, 0, head)
  naxis = sxpar(head, 'NAXIS')
  if naxis NE 1 then begin
      print, 'NAXIS NE 1', naxis
      return
  endif

  if keyword_set( ERROR ) then sig = mrdfits(error)

; Make Wavelength array

  wave = x_fitswave(head)

; Cut on wavelengths as necessary

  goodwv = where( wave LE wvmax AND wave GE wvmin)
  

; Output

  if keyword_set( ERROR ) then writecol, FMT='(F12.6,2F10.4)', outfil, $
    wave[goodwv], data[goodwv], sig[goodwv] $
    else writecol, FMT='(F12.6, F10.4)', outfil, wave[goodwv], data[goodwv]

end

