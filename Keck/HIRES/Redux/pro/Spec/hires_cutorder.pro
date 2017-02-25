;+ 
; NAME:
; hires_cutorder
;     Version 1.1
;
; PURPOSE:
;   Cut out an order (or two) and write to disk as a simple
;   multi-extension FITS file
;
; CALLING SEQUENCE:
;   hires_cutorder, file, ORDRS 
;
; INPUTS:
;    file - Filename to grab from
;    ORDRS - Orders to cut
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
;   17-Oct-2011 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro hires_cutorder, file, ORDRS, OUTROOT=outroot

  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'hires_cutorder, file, ordrs, OUTROOT=outroot [v1.0]'
      return
  endif 

  ;; Outroot
  if not keyword_set(OUTROOT) then begin
     ipos = strpos(file, '.fits')
     outroot = strmid(file,0,ipos)
  endif

  ;; Open the file
  spec2d = xmrdfits(file,1,/silent)
  nord = n_elements(ORDRS)
  gdo = lonarr(nord)
  for jj=0L,nord-1 do begin
     gdo[jj] = where(spec2d.phys_ordr EQ ordrs[jj])
  endfor
          
  sz = size(spec2d.fx,/dimensions)

  ;; Loop
  for jj=0L,nord-1 do begin
      qq = gdo[jj]
      ;; Grab good data
      gd = where(spec2d.wave[*, qq] GT 0. AND $
                 spec2d.var[*, qq] GT 0., ngd)
      ;; Start/end pix
      i0 = gd[0]
      i1 = gd[ngd-1]
      np = i1-i0+1
      gdp = i0 + lindgen(np)

      ;; Filename
      outfil = outroot + '_'+ x_padstr(strtrim(ordrs[jj],2),3,'0',/rever)+'.fits'

      ;; Create dummy file
      mwrfits, spec2d.fx[gdp,qq], outfil, /create
      mwrfits, sqrt((spec2d.var[gdp,qq]>0.)), outfil
      mwrfits, spec2d.wave[gdp,qq], outfil
      print, 'hires_cutorder: Writing: ', outfil
   endfor

  return 

end
      
