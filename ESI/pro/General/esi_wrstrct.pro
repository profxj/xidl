;+ 
; NAME:
; esi_wrstrct   
;    Version 1.1
;
; PURPOSE:
;    Write the esi structure to a FITS file and write an ASCII summary
;
; CALLING SEQUENCE:
;   
;  esi_wrstrct, esi, /ANONLY, OUTFIL=, FITS=
;
; INPUTS:
;   esi   - An ESI structure
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   ANONLY - Only print files with flg_anly NE 0   
;
; OPTIONAL OUTPUTS:
;   OUTFIL= - Output file (default: esi.list)
;   FITS=   - Name of fits output file
;
; COMMENTS:
;
; EXAMPLES:
;   esi_wrstrct, esi, FITS='esi_13oct02.fits'
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   13-Jul-2002 Written by JXP
;   29-Jan-2003  Polished (JXP)
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro esi_wrstrct, esi, ANONLY=anonly, OUTFIL=outfil, FITS=fits

;
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'esi_wrstrct, esi, /ANONLY, OUTFIL=, FITS= (v1.1)'
      return
  endif 
  
;  Optional Keywords
  
  if not keyword_set( OUTFIL ) then outfil = 'esi.list'

  nimg = n_elements(esi)
  if nimg EQ 1 then begin
      print, 'esi_wrstrct: I dont think you passed a structure in'
      return
  endif
  
  close, /all
  openw, 1, outfil

  for q=0,nimg-1 do begin
      
      if esi[q].obj_id GE 0L then begin
          if esi[q].obj_id LT 10 then id_str = '0'+strtrim(esi[q].obj_id,2) $
          else id_str = strtrim(esi[q].obj_id,2)
      endif else id_str = '-1'
      case esi[q].mode of
          0: mode = 'IMG'
          1: mode = 'LWD'
          2: mode = 'ECH'
          else:
      endcase
      if (not keyword_set( ANONLY ) OR esi[q].flg_anly NE 0) then $
        printf, 1, $
        FORMAT='(i4,1x,a7,1x,i1,1x,a12,1x,a3,1x,a4,a5,i5,f5.2,i2,i2,a2,f5.2,2i2)',$
        q, $
        esi[q].img_root, $
        esi[q].flg_anly, $
        esi[q].Obj, $
        id_str, $
        mode, $
        esi[q].type, $
        long(esi[q].exp), $
        esi[q].slit, $
        esi[q].arclamp, $
        esi[q].qtzlamp, $
        esi[q].imfilt, $
        esi[q].AM, $
        esi[q].flg_ov, $
        esi[q].flg_final
        
      
  endfor

  close, 1

  if keyword_set( FITS ) then mwrfits, esi, fits, /create


end
