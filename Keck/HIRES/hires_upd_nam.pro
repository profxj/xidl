;+ 
; NAME:
; hires_upd_nam
;    Version 1.1
;
; PURPOSE:
;   Set the FILENAME, SIGFILE, and CONTFILE header key cards
;
; CALLING SEQUENCE:
;   hires_upd_nam, fil, filnm, signm, contnm
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
pro hires_upd_nam, fil, filnm, signm, contnm

  if  N_params() LT 4  then begin 
      print,'Syntax - ' + $
        'hires_upd_nam, fil, filnm, signm, contnm [v1.1]'
      return
  endif 
	img = readfits(fil,head)
    	fxaddpar, head, 'FILENAME', filnm
    	fxaddpar, head, 'SIGFILE ', signm
    	fxaddpar, head, 'CONTFILE', contnm
	mwrfits, img, fil, head, /create
return
end
