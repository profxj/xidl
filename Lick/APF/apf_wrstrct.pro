;+ 
; NAME:
; apf_wrstrct   
;    Version 1.1
;
; PURPOSE:
;    Write the apf structure to a FITS file and write an ASCII summary
;
; CALLING SEQUENCE:
;   
;  apf_wrstrct, apf, /ANONLY, OUTFIL=, FITS=
;
; INPUTS:
;   apf   - HIRES structure
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   ANONLY - Only print files with flg_anly NE 0   
;
; OPTIONAL OUTPUTS:
;   OUTFIL= - Output file (default: 'apf.list')
;   FITS=   - Name of fits output file (default: 'strct.fits')
;
; COMMENTS:
;
; EXAMPLES:
;   apf_wrstrct, apf, FITS='apf_13oct02.fits'
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   24-Jan-2005 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro apf_wrstrct, apf, ANONLY=anonly, OUTFIL=outfil, FITS=fits

;
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'apf_wrstrct, apf, /ANONLY, OUTFIL=, FITS= [v1.0]'
      return
  endif 
  
;  Optional Keywords
  
  if not keyword_set( OUTFIL ) then outfil = 'apf.list'
  if not keyword_set( FITS ) then fits = 'strct.fits'

  nimg = n_elements(apf)
  if nimg EQ 1 then begin
      print, 'apf_wrstrct: I dont think you passed a structure in'
      return
  endif
  
  close, /all
  openw, 1, outfil

  for q=0,nimg-1 do begin
      
      if apf[q].obj_id GE 0L then begin
          if apf[q].obj_id LT 10 then id_str = '0'+strtrim(apf[q].obj_id,2) $
          else id_str = strtrim(apf[q].obj_id,2)
      endif else id_str = '-1'
      ;; CCD
      case apf[q].chip of 
          -1: c_chip = 'S' ; Single
          1: c_chip = 'B'
          2: c_chip = 'G'
          3: c_chip = 'R'
          else: stop
      endcase
      ;; Print
      if (not keyword_set( ANONLY ) OR apf[q].flg_anly NE 0) then $
        printf, 1, $
        FORMAT='(i4,1x,i5,1x,a20,1x,i1,1x,a12,1x,i2,a4,a5,i5,a2,i2,i2,a6,a6,' + $
        'f8.4,f8.4,1x,a7,a6)', $
        q, $
        apf[q].frame, $
        apf[q].img_root, $
        apf[q].flg_anly, $
        strtrim(apf[q].Obj,2), $
        apf[q].setup, $
        id_str, $
        apf[q].type, $
        round(apf[q].exp), $
        c_chip, $
        apf[q].colbin, $
        apf[q].rowbin, $
        apf[q].decker, $
        apf[q].block, $
        apf[q].echangl, $
        apf[q].xdangl, $
        apf[q].lamp, $
        apf[q].lampfil
        
      
  endfor

  close, 1

  ;; Write fits file
  mwrfits, apf, fits, /create


end
