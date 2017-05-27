function deimos_read_file, filename, nobias=nobias
;+
; NAME:
;    deimos_read_file
;
; PURPOSE:
;    reads a full DEIMOS red mosaic frame, putting it together into
;    8K*8K arrary
;
; CALLING SEQUENCE:
;    image = deimos_read_file(filename)
; 
; INPUTS:
;    filename -- string full path or current directory
;
; OPTIONAL INPUTS:
;	
; KEYWORDS:
;  nobias -- set to not bias subtract
; OUTPUTS:
;   
; COMMENTS:
;   assumes for the moment that data is in 8 amp mode
;
; REVISION HISTORY:
;   md 30apr02
;
;----------------------------------------------------------------------
  if n_elements(nobias) eq 0 then nobias = 0
  image = fltarr(8192, 8192)
  ii = 2048
  jj = 4096
  
  image[0:ii-1,      0:jj-1] = deimos_read_chip(filename, 1, /noconvert, nobias=nobias)
  image[ii:2*ii-1,   0:jj-1] = deimos_read_chip(filename, 2, /noconvert, nobias=nobias)
  image[2*ii:3*ii-1, 0:jj-1] = deimos_read_chip(filename, 3, /noconvert, nobias=nobias)
  image[3*ii:4*ii-1, 0:jj-1] = deimos_read_chip(filename, 4, /noconvert, nobias=nobias)
  image[0:ii-1,      jj:2*jj-1] = deimos_read_chip(filename, 5, /noconvert, nobias=nobias)
  image[ii:2*ii-1,   jj:2*jj-1] = deimos_read_chip(filename, 6, /noconvert, nobias=nobias)
  image[2*ii:3*ii-1, jj:2*jj-1] = deimos_read_chip(filename, 7, /noconvert, nobias=nobias)
  image[3*ii:4*ii-1, jj:2*jj-1] = deimos_read_chip(filename, 8, /noconvert, nobias=nobias)


  return,  image
end
