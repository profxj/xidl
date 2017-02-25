;+ 
; NAME:
; uves_wrstrct   
;    Version 1.1
;
; PURPOSE:
;    Write the uves structure to a FITS file and write an ASCII summary
;
; CALLING SEQUENCE:
;   
;  uves_wrstrct, uves, /ANONLY, OUTFIL=, FITS=
;
; INPUTS:
;   uves   - An ESI structure
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   ANONLY - Only print files with flg_anly NE 0   
;
; OPTIONAL OUTPUTS:
;   OUTFIL= - Output file (default: 'uves.list')
;   FITS=   - Name of fits output file (default: 'strct.fits')
;
; COMMENTS:
;
; EXAMPLES:
;   uves_wrstrct, uves, FITS='uves_13oct02.fits'
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   24-Jan-2005 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro uves_wrstrct, uves, ANONLY=anonly, OUTFIL=outfil, FITS=fits

;
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'uves_wrstrct, uves, /ANONLY, OUTFIL=, FITS= [v1.0]'
      return
  endif 
  
;  Optional Keywords
  
  if not keyword_set( OUTFIL ) then outfil = 'uves.list'
  if not keyword_set( FITS ) then fits = 'strct.fits'

  nimg = n_elements(uves)
  if nimg EQ 1 then begin
      print, 'uves_wrstrct: I dont think you passed a structure in'
      return
  endif
  
  close, /all
  openw, 1, outfil

  for q=0,nimg-1 do begin
      
      if uves[q].obj_id GE 0L then begin
          if uves[q].obj_id LT 10 then id_str = '0'+strtrim(uves[q].obj_id,2) $
          else id_str = strtrim(uves[q].obj_id,2)
      endif else id_str = '-1'
      ;; CCD
      case uves[q].side of 
          1: c_side = 'B'
          2: c_side = 'R'
          else: stop
      endcase
      ;; Print
      if (not keyword_set( ANONLY ) OR uves[q].flg_anly NE 0) then $
        printf, 1, $
        FORMAT='(i4,1x,i4,1x,a9,1x,i1,1x,a12,1x,i2,a4,a5,i5,a2,i2,i2,1x,f4.2,1x,a6,' + $
        'f8.4,f8.4,1x,a7,a6)', $
        q, $
        uves[q].frame, $
        uves[q].img_root, $
        uves[q].flg_anly, $
        strtrim(uves[q].Obj,2), $
        uves[q].setup, $
        id_str, $
        uves[q].type, $
        round(uves[q].exp), $
        c_side, $
        uves[q].colbin, $
        uves[q].rowbin, $
        uves[q].slitwid, $
        uves[q].block, $
        uves[q].echangl, $
        uves[q].xdangl, $
        uves[q].lamp, $
        uves[q].lampfil
        
      
  endfor

  close, 1

  ;; Write fits file
  mwrfits, uves, fits, /create


end
