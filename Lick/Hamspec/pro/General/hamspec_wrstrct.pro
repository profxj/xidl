;+ 
; NAME:
; hamspec_wrstrct   
;    Version 1.1
;
; PURPOSE:
;    Write the hamspec structure to a FITS file and write an ASCII summary
;
; CALLING SEQUENCE:
;   
;  hamspec_wrstrct, hamspec, /ANONLY, OUTFIL=, FITS=
;
; INPUTS:
;   hamspec   - HIRES structure
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   ANONLY - Only print files with flg_anly NE 0   
;
; OPTIONAL OUTPUTS:
;   OUTFIL= - Output file (default: 'hamspec.list')
;   FITS=   - Name of fits output file (default: 'strct.fits')
;
; COMMENTS:
;
; EXAMPLES:
;   hamspec_wrstrct, hamspec, FITS='hamspec_13oct02.fits'
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   24-Jan-2005 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro hamspec_wrstrct, hamspec, ANONLY=anonly, OUTFIL=outfil, FITS=fits

;
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'hamspec_wrstrct, hamspec, /ANONLY, OUTFIL=, FITS= [v1.0]'
      return
  endif 
  
;  Optional Keywords
  
  if not keyword_set( OUTFIL ) then outfil = 'hamspec.list'
  if not keyword_set( FITS ) then fits = 'strct.fits'

  nimg = n_elements(hamspec)
  if nimg EQ 1 then begin
      print, 'hamspec_wrstrct: I dont think you passed a structure in'
      return
  endif
  
  close, /all
  openw, 1, outfil

  for q=0,nimg-1 do begin
      
      if hamspec[q].obj_id GE 0L then begin
          if hamspec[q].obj_id LT 10 then id_str = '0'+strtrim(hamspec[q].obj_id,2) $
          else id_str = strtrim(hamspec[q].obj_id,2)
      endif else id_str = '-1'
      ;; Print
      if (not keyword_set( ANONLY ) OR hamspec[q].flg_anly NE 0) then $
        printf, 1, $
        FORMAT='(i4,1x,i4,1x,a20,1x,i1,1x,a12,1x,i2,a4,a5,i5,i2,i2,1x,i5,1x,f5.2,1x,a7,' + $
        'i6,i6,1x,a7)', $
        q, $
        hamspec[q].frame, $
        hamspec[q].img_root, $
        hamspec[q].flg_anly, $
        strtrim(hamspec[q].Obj,2), $
        hamspec[q].setup, $
        id_str, $
        hamspec[q].type, $
        round(hamspec[q].exp), $
        hamspec[q].colbin, $
        hamspec[q].rowbin, $
        hamspec[q].width, $
        hamspec[q].length, $
        hamspec[q].block, $
        hamspec[q].echangl, $
        hamspec[q].xdangl, $
        hamspec[q].lamp
  endfor

  close, 1

  ;; Write fits file
  mwrfits, hamspec, fits, /create


end
