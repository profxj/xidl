;+ 
; NAME:
; hires_cut_35
;    Version 1.1
;
; PURPOSE:
;    Zero out the first 35 points of a HIRES spectrum.  I cant
;  even recall why I do it.
;
; CALLING SEQUENCE:
;   hires_cut_35, fil
;
; INPUTS:
;   fil -- Name of HIRES file to edit
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   ??  Written by JXP
;-
;------------------------------------------------------------------------------
pro hires_cut_35, fil1

  if  N_params() LT 5  then begin 
      print,'Syntax - ' + $
        'hires_cut_35, fil [v1.1]'
      return
  endif 
	img = readfits(fil1, head)
	img[0:34,*] = -1.
	writefits, fil1, img, head
return
end
