;+ 
; NAME:
; imacsls_wrstrct   
;    Version 1.1
;
; PURPOSE:
;    Write the imacsls structure to a FITS file and write an ASCII summary
;
; CALLING SEQUENCE:
;   
;  imacsls_wrstrct, imacsls, /ANONLY, OUTFIL=, FITS=
;
; INPUTS:
;   imacsls   - IMACS long slit structure
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   /ANONLY - Only print files with flg_anly NE 0   
;
; OPTIONAL OUTPUTS:
;   OUTFIL= - Output file (default: imacsls.list)
;   FITS=   - Name of fits output file
;
; COMMENTS:
;
; EXAMPLES:
;   imacsls_wrstrct, imacsls, FITS='imacsls_13oct02.fits'
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   08-Dec-2003 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro imacsls_wrstrct, imacsls, ANONLY=anonly, OUTFIL=outfil, FITS=fits

;
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'imacsls_wrstrct, imacsls, /ANONLY, OUTFIL=, FITS= (v1.0)'
      return
  endif 
  
;  Optional Keywords
  
  if not keyword_set( OUTFIL ) then outfil = 'imacsls.list'

  nimg = n_elements(imacsls)
  if nimg EQ 1 then begin
      print, 'imacsls_wrstrct: I dont think you passed a structure in'
      return
  endif
  
  close, /all
  openw, 1, outfil

  for q=0,nimg-1 do begin
      
      if imacsls[q].obj_id GE 0L then begin
          if imacsls[q].obj_id LT 10 then id_str = '0'+strtrim(imacsls[q].obj_id,2) $
          else id_str = strtrim(imacsls[q].obj_id,2)
      endif else id_str = '-1'
      case imacsls[q].mode of
          0: mode = 'IMG'
          1: mode = 'SPEC'
          else:
      endcase
      if (not keyword_set( ANONLY ) OR imacsls[q].flg_anly NE 0) then $
        printf, 1, $
        FORMAT='(i4,1x,a7,1x,i1,1x,a12,1x,a4,1x,a4,1x,i2,a4,a10,f8.3,i5,f5.2,4i2,f5.2,2i2)',$
        q, $
        imacsls[q].img_root, $
        imacsls[q].flg_anly, $
        imacsls[q].Obj, $
        id_str, $
        mode, $
        imacsls[q].side, $
        imacsls[q].type, $
        imacsls[q].grising, $
        imacsls[q].grangle, $
        long(imacsls[q].exp), $
        imacsls[q].slit, $
        imacsls[q].rbin, $
        imacsls[q].cbin, $
        imacsls[q].arclamp, $
        imacsls[q].qtzlamp, $
        imacsls[q].AM, $
        imacsls[q].flg_ov, $
        imacsls[q].flg_final
        
      
  endfor

  close, 1

  if keyword_set( FITS ) then mwrfits, imacsls, fits, /create


end
