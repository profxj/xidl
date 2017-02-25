;+
; 
; NAME:
; xdimg_finddark
;    Version 1.1
;
; PURPOSE:
;    Finds the superdark file associated with the requested images.
;    Updates the structure with the flg_drk and img_drk keywords.
;
; CALLING SEQUENCE:
;   xdimg_finddark, struct, imgarr
;
; INPUTS:
;   struct -- dimg_strct defining the images of interest
;
; RETURNS:
;
; OUTPUTS:
;
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;     xdimg_finddark, struct, imgarr
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
;
; REVISION HISTORY:
;   21-June-2007 Written by LKP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


pro xdimg_finddark, struct, imgarr, SVSTRCT=SVSTRCT

  if  N_params() LT 1  then begin 
     print,'Syntax - ' + $
           'xdimg_finddark, struct, imgarr, /SVSTRCT (v1.1)'
     return
  endif 

for i=0,n_elements(imgarr)-1 do begin

   exp_str = strtrim(round(struct[imgarr[i]].exp*1000.)/1000., 2)
   exp_pos = strpos(exp_str, '.')
   exp_str = strmid(exp_str, 0, exp_pos+4)

   if struct[imgarr[i]].ccd eq 'NIRC2N' then camera = 'n' $
        else if struct[imgarr[i]].ccd eq 'NIRC2W' then camera = 'w' $
        else camera = ''

   coadds = strtrim(struct[imgarr[i]].coadds,2)
   multisam = strtrim(struct[imgarr[i]].multisam,2)
   sampmode = strtrim(struct[imgarr[i]].sampmode,2)

   lookfile =  'Darks/Dark'+'_'+camera+exp_str+'_'+coadds+'_'+sampmode+'_'+multisam+'.fits'

   a = findfile(lookfile, count=count)
   if count eq 1 then begin
      struct[imgarr[i]].flg_drk = 1
      struct[imgarr[i]].img_drk = lookfile
   endif else begin
      print, "We don't have darks for the following setup:"
      print, '     camera='+camera+',  itime='+exp_str+',  coadds='+coadds+',  sampmode=' $
             +sampmode+',  multisam='+multisam

   endelse

endfor

; Resave the updated structure so that you can pick up where you left off easily.
if keyword_set(SVSTRCT) then mwrfits, struct, 'struct.fits', /create

end

