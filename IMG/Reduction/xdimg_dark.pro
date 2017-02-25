;+
; 
; NAME:
; xdimg_dark
;    Version 1.1
;
; PURPOSE:
;    Creates super dark frame given direct image structure.
;
; CALLING SEQUENCE:
;   xdimg_dark, struct
;
; INPUTS:
;   struct -- dimg_strct defining the images of interest
;
; RETURNS:
;
; OUTPUTS:
;   OUTFIL - Name of Dark file [default is e.g. 'Darks/Dark_n10.000_1_3_8.fits']
;            The 10.000 is the exposure time, the "n" indicates the narrow camera,
;            the "1_3_8" means coadds=1, sampmode=3, multisam=8.
;
; OPTIONAL KEYWORDS:
;
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   xdimg_dark, nght1_strct
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
;
; REVISION HISTORY:
;   11-June-2007 Written by LKP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro xdimg_dark, struct

  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'xdimg_dark, struct (v1.1)'
      return
  endif 
  
;  Optional Keywords
  if not keyword_set( OUTFIL ) then begin
      ; Create directory if necessary
      a = findfile('Darks/..', count=count)
      if count EQ 0 then file_mkdir, 'Darks', 'Darks/Sub'
  endif
  
;  Find the Dark frames
  dark = where(struct.type EQ 'DRK' AND struct.flg_anly NE 0, ndark)

  if ndark eq 0 then begin
     print, 'There are no dark frames to process!'
     return
  endif
  
;  Loop over all dark frames to make a superdark for each type of dark available.
;  NOTE: This means, a separate superdark for each exposure time, ncoadd, nread, sampmode, camera!
  for i=0,ndark - 1 do begin
     exp = round(struct[dark[i]].exp * 1000.)/1000.  ; Round exposure time to nearest 1/1000 sec.
     exp_str = strtrim(exp,2)
     exp_pos = strpos(exp_str, '.')
     exp_str = strmid(exp_str, 0, exp_pos+4)

     coadds = struct[dark[i]].coadds
     sampmode = struct[dark[i]].sampmode
     multisam = struct[dark[i]].multisam

     if struct[dark[i]].ccd eq 'NIRC2N' then camera = 'n' $
        else if struct[dark[i]].ccd eq 'NIRC2W' then camera = 'w' $
        else camera = ''
     
     outfil = 'Darks/Dark'+'_'+camera+exp_str+'_'+strtrim(coadds,2)+'_'+strtrim(sampmode,2)+'_'+strtrim(multisam,2)+'.fits'
     
     ; Check whether we've already made a dark with this exposure time.
     a = findfile(outfil, count=count)
     if count GT 0 then continue else print, 'Working on  '+outfil
     
     ; Find all of the darks with the specific camera, exposure time, coadds, sampmode, multisam setup.
     darks_exp = where( round(struct[dark].exp * 1000.)/1000. eq exp AND $
                        struct[dark].ccd eq struct[dark[i]].ccd AND $
                        struct[dark].coadds eq coadds AND $
                        struct[dark].multisam eq multisam, ndark_exp)
     xdim = 1024 & ydim = 1024
     dark_stack = fltarr(xdim,ydim,ndark_exp)
     
     ; Get stack of all darks of a specific type.
     for j=0, ndark_exp-1 do begin

        darkfil = struct[dark[darks_exp[j]]].rootpth+struct[dark[darks_exp[j]]].img_root
        darkfil = mrdfits(darkfil, 0, header, /silent)

        ; Currently can only make stack for 1024x1024 images.
        darksz = size(darkfil, /dim)
        if darksz[0] ne xdim or darksz[1] ne ydim then begin
           print, 'Incompatible image sizes.  Aborting!'
           return
        endif

        dark_stack[*,*,j]=darkfil

     endfor

     ; For now we're doing nothing fancy.
     ; Just assume the median of the dark stack is the true dark, 
     ; without doing any sigma clipping.

     if ndark_exp gt 1 then begin
        dark_compress=djs_median(dark_stack, 3)  ;should I use median of some small region for combining?
     endif else begin
        print, 'Only one dark of this type.  Not median combining.'
        dark_compress=dark_stack[*,*,0]
     endelse

     ; Output
     if ndark_exp ge 1 then begin
        mwrfits, dark_compress, outfil, /create
        print, 'Used ' + strtrim(ndark_exp,2) + ' individual dark frames.'

     endif
        
  endfor

; Resave the updated structure so that you can pick up where you left off easily.
  mwrfits, struct, 'struct.fits', /create

  print, 'All done with Dark frames!'

end
