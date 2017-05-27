;+
; NAME:
;   deimos_read_chip
;
; PURPOSE:
;   Read, bias-subtract, and rotate image from one chip
;
; CALLING SEQUENCE:
;   image = deimos_read_chip(fname, chipno, nobias=nobias, $
;                            norotate=norotate, header=header)
;
; INPUTS:
;   fname    - FITS file name
;   chipno   - chip number (corresponds to CCDLOC in header)
;
; KEYWORDS:
;   nobias    - set to suppress bias subtraction
;   norotate  - set to suppress image rotation
;   noconvert - set to suppress linearization and conversion to e-
;   fixup     - fix up electronics problem in chip 5
;   quick     - only apply gain correction, but not linearity
;
; OUTPUTS:
;   image    - <returned> 1-chip image, bias-subtracted and rotated
;
; OPTIONAL OUTPUTS:
;   header   - header of relevent extension (1st HDU in case of 2 HDU/chip)
;   satmask  - mask of saturated pixels (A/D saturation at 65535)
;
; RESTRICTIONS:
;
; EXAMPLES:
;   image = deimos_read_chip('GridOfHoles_8amp_Mosaic.fits', 3, /norotate)
;
; COMMENTS:
;
; REVISION HISTORY:
;
;   2002-Feb-27  Written by D. Finkbeiner, Princeton
;   2002-Mar-27  16 amp mode, read only 1st block of headers
;   2002-Apr-17  nonlinear correction, mapping to correct HDU
;   2002-Jul-14  Saturation mask
;   2002-Oct-06  Read dual-amp mode correctly - DPF
;   2002-Oct-17  Rotate the saturation mask also!
;   2002-Oct-18  Added checksum test - DPF
;
;----------------------------------------------------------------------
function deimos_read_1amp, fname, hno, nobias=nobias, norotate=norotate, $
               satmask=satmask, detsec=detsec, header=header

; -------- Read in the image from that header as unsigned int (if approp)
  fullimage = mrdfits(fname, hno, /unsigned, /silent)

; -------- Read in header (do NOT use the header altered by mrdfits!)
  header = headfits(fname, ext=hno)

; -------- Check for data for corruption
  datasum = ulong(sxpar(header, 'DATASUM', count=count))
  if count EQ 1 then begin 
     data = fullimage+32768US   ; flip high bit back for checksum
     host_to_ieee, data
     ulongdata = ulong(data,0,n_bytes(data)/4)
     host_to_ieee, ulongdata
     checksum32, ulongdata, sum32
     delvarx, data, ulongdata    ; free some memory
     if datasum NE sum32 then begin 
        print, '  File name: ', fname
        print, '  Header number:', hno
        print, '  DATASUM (header):', datasum, '    computed checksum:', sum32
        message, 'Data may be corrupted!!'
     endif else begin 
        print, 'CHECKSUM correct: ', sum32
     endelse 
  endif else begin 
     print, 'NOTE: no checksum field in header - cannot test data integrity'
  endelse 

; -------- Check header for corruption
  checksum = sxpar(header, 'CHECKSUM', count=countc)
  if countc EQ 1 then begin 
     nline = n_elements(header) 
     nmod = nline mod 36   ; units of 2880 bytes - FITS standard
     if nmod NE 0 then begin 
        npad = 36-nmod
        blankline = string(bytarr(80)+32B)
        pad = replicate(blankline, npad)
        hpad = [header, pad]
     endif else hpad = header
     byte_header = byte(hpad)
     ulong_header = ulong(byte_header, 0, n_bytes(byte_header)/4)
     host_to_ieee, ulong_header
     if count EQ 0 then message, 'Header has CHECKSUM but not DATASUM!'
     checksum32, [sum32, ulong_header], headsum
     if (headsum + 1UL) NE 0UL then begin 
        print, '  File name: ', fname
        print, '  Header number:', hno
        print, '  Computed header checksum (should be 4294967295):', headsum
        message, 'Corrupted file header!!'
     endif else begin 
        print, '  Header intact'
     endelse
  endif 

; -------- Extract datasec from header
  datasec = sxpar(header, 'DATASEC',count=datact)
  detsec  = sxpar(header, 'DETSEC')

; -------- Parse datasec using a (clumsy) regexp
  x = (stregex(datasec, '\[([0-9]+):([0-9]+),([0-9]+):([0-9]+)\]', /subexp, /extract))[1:4] - 1

; -------- Pull out science image and overscan
  image = fullimage[x[0]:x[1], x[2]:x[3]]

  israwdata = n_elements(image) ne n_elements(fullimage)
; -------- if there is no overscan, subtract no bias
  if israwdata then over  = fullimage[x[1]+1:*, *] $
	else over=0.

; -------- Check for A/D saturation (at 65535; no check for full well)
  if arg_present(satmask) then begin 
     satmask = image EQ 65535U
     satfrac = total(satmask)/n_elements(satmask)
     if satfrac GT 0.01 then message, $
       'WARNING:  Saturated image!  Fraction saturated:'+ $
       string(satfrac*100, format='(F6.2,"%")'),/info
  endif 

; -------- Check that there is some signal in the image(!)
  if stdev(image[0:99999]) eq 0 then message, $
    'WARNING:  There is no signal in this image!', /info

; -------- Bias subtract (unless requested not to) (before rotation!)
;               We should make this position dependent. 
  if (NOT keyword_set(nobias)) AND israwdata then begin
     djs_iterstat, over,  mean=bias, sigrej=5.
     image = image- bias
  endif

; -------- Rotate (unless requested not to)
  if (NOT keyword_set(norotate)) AND israwdata  then begin 
     detsec1 = detsec
     detsec_rotate, detsec, image
     if keyword_set(satmask) then detsec_rotate, detsec1, satmask
  endif 

  return, image
end



function deimos_read_chip, fname, chipno, satmask=satmask, nobias=nobias, $
         norotate=norotate, noconvert=noconvert, fixup=fixup, $
               header=header, quick=quick


  if n_elements(quick) eq 0 then quick = 0


  phdr = headfits(fname) ; primary FITS header
  numamps = sxpar(phdr, 'NVIDINP') ; numamps 
  ampmode = sxpar(phdr, 'AMPMODE') ; SINGLE:A / SINGLE:B / DUAL
; TBD -- need to decode next keyword AMPLIST = '8,1,0,0 ' to properly
;        decide which amplifiers to use.
  single = strmid(ampmode,0,6) EQ 'SINGLE' ; is it single-amp mode?

;check to see if mapping from amplifier to HDU has been done.
  if single then defsysv, '!deimos_hduamp1', exists=exists $
	else defsysv, '!deimos_hduamp2', exists=exists

  if NOT  exists then deimos_red_mosaic, ampmode
  
	if single then deimos_hduamp = !deimos_hduamp1 $
	  else deimos_hduamp = !deimos_hduamp2

;read subimage
  hdu = (single ? 1 : 2)*(chipno) ; 1,3,5... or 1,2,3...
  subimage = deimos_read_1amp(fname, hdu, nobias=nobias, norotate=norotate, $
       satmask=satmask, detsec=detsec, header=header)
;convert to electrons
;does this work in 8 amp as well as 16 amp mode??

; -------- fix up chip 5 (this is a hack)
  if keyword_set(fixup) AND chipno eq 5 then begin 
     chip6 = deimos_read_chip(fname, 6, /noconvert)
     subimage = subimage+((chip6-median(chip6)) < 50)*.04
  endif 

  if NOT keyword_set(noconvert) then $
    subimage = deimos_adu2e(subimage,deimos_hduamp[hdu-1], single=single, $
                            quick=quick)

  if single then return, subimage

; -------- For 2-amp readout:

  subimage1 = deimos_read_1amp(fname,deimos_hduamp[hdu-2], $
      nobias=nobias, norotate=norotate, satmask=satmask1, detsec=detsec)
  satmask = [satmask, satmask1]

  if NOT keyword_set(noconvert) then  $
    subimage1 = deimos_adu2e(subimage1,deimos_hduamp[hdu-2], single=single, $
                            quick=quick)

  return, [subimage,subimage1]
end

  




