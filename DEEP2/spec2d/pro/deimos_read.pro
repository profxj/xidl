;+
; NAME:
;   deimos_read
;
; PURPOSE:
;   Read one chip from a DEIMOS 8-chip file
;
; CALLING SEQUENCE:
;   image = deimos_read(fname, chipno,[header=header], [flat=flat] )
;
; INPUTS:
;   fname      - file name (with path, if required)
;   chipno     - chip number (1-8). 
;   
; OPTIONAL INPUTS: 
;   flat       - flat image: if supplied then flat fielding is performed
; 
; OUTPUTS:
;   image      - requested image array[2048, 4096]
;   header     - optional FITS header
;
; COMMENTS: 
;   arrangement of chips assumed to be:   5678
;                                         1234
;   All image locations refer to a display where pixel [0,0] is lower
;   left.  All arrays are 0-indexed in IDL. 
; 
;   Number of prescan and overscan pixels taken from header. 
;   This routine does NOT bias subtract. 
; 
; EXAMPLES:
;   arc = deimos_read('arcimage.fits', 4)
;
; REVISION HISTORY:
;   01-Dec-2000  Written by D. Finkbeiner, Berkeley
;  
;-
;----------------------------------------------------------------------------


function deimos_read, fname, chipno, header=header, flat=flat

; number of chips in array (4x2 for DEIMOS)
  nchip = [4, 2]
  header = headfits(fname)

; extract useful information from the header
  naxis1  = sxpar(header, 'NAXIS1') 
  naxis2  = sxpar(header, 'NAXIS2')
  prepix  = sxpar(header, 'PREPIX')   ; number of prescan cols per chip (on left)
  if prepix EQ 0 then prepix = sxpar(header, 'PRECOL')
  if prepix EQ 0 then stop
  postpix = sxpar(header, 'POSTPIX')  ; number of postscan cols per chip (on rt)
;
; 
; these numbers are wrong in DEIMOS fits headers at the moment. FUDGE IT!
  prepix = 184/4
  postpix = (naxis1-naxis2)/4 -prepix

; chip numbering - 1-4 across bottom (left-right) (5-8 on top)
  i = (chipno-1) mod nchip[0]
  j = (chipno-1)  /  nchip[0]

; size of data area for each chip
  dx = naxis1/nchip[0]-(prepix+postpix)
  dy = naxis2/nchip[1]

; determine location of desired chip within full image
  x0 = nchip[0]*prepix+i*dx
  y0 = j*dy

; read subimage from FITS file
  fxread, fname, im, h, x0, x0+dx-1, y0, y0+dy-1

  return, im
end



