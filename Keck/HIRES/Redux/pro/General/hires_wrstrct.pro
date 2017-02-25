;+ 
; NAME:
; hires_wrstrct   
;    Version 1.1
;
; PURPOSE:
;    Write the hires structure to a FITS file and write an ASCII summary
;
; CALLING SEQUENCE:
;   
;  hires_wrstrct, hires, /ANONLY, OUTFIL=, FITS=
;
; INPUTS:
;   hires   - HIRES structure
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   ANONLY - Only print files with flg_anly NE 0   
;
; OPTIONAL OUTPUTS:
;   OUTFIL= - Output file (default: 'hires.list')
;   FITS=   - Name of fits output file (default: 'strct.fits')
;
; COMMENTS:
;
; EXAMPLES:
;   hires_wrstrct, hires, FITS='hires_13oct02.fits'
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   24-Jan-2005 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro hires_wrstrct, hires, ANONLY=anonly, OUTFIL=outfil, FITS=fits

;
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'hires_wrstrct, hires, /ANONLY, OUTFIL=, FITS= [v1.0]'
      return
  endif 
  
;  Optional Keywords
  
  if not keyword_set( OUTFIL ) then outfil = 'hires.list'
  if not keyword_set( FITS ) then fits = 'strct.fits'

  nimg = n_elements(hires)
  if nimg EQ 1 then begin
      print, 'hires_wrstrct: I dont think you passed a structure in'
      return
  endif
  
  close, /all
  openw, 1, outfil

  for q=0,nimg-1 do begin
      
      if hires[q].obj_id GE 0L then begin
          if hires[q].obj_id LT 10 then id_str = '0'+strtrim(hires[q].obj_id,2) $
          else id_str = strtrim(hires[q].obj_id,2)
      endif else id_str = '-1'
      ;; CCD
      case hires[q].chip of 
          -1: c_chip = 'S' ; Single
          1: c_chip = 'B'
          2: c_chip = 'G'
          3: c_chip = 'R'
          else: stop
      endcase
      ;; Print
      if (not keyword_set( ANONLY ) OR hires[q].flg_anly NE 0) then $
        printf, 1, $
        FORMAT='(i4,1x,i4,1x,a20,1x,i1,1x,a12,1x,i2,a4,a5,i5,a2,i2,i2,a4,a6,' + $
        'f8.4,f8.4,1x,a7,a6)', $
        q, $
        hires[q].frame, $
        hires[q].img_root, $
        hires[q].flg_anly, $
        strtrim(hires[q].Obj,2), $
        hires[q].setup, $
        id_str, $
        hires[q].type, $
        round(hires[q].exp), $
        c_chip, $
        hires[q].colbin, $
        hires[q].rowbin, $
        hires[q].decker, $
        hires[q].block, $
        hires[q].echangl, $
        hires[q].xdangl, $
        hires[q].lamp, $
        hires[q].lampfil
        
      
  endfor

  close, 1

  ;; Write fits file
  mwrfits, hires, fits, /create


end
