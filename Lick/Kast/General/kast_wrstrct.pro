;+ 
; NAME:
; kast_wrstrct   
;    Version 1.0
;
; PURPOSE:
;    Write the kast structure to a FITS file and write an ASCII summary
;
; CALLING SEQUENCE:
;   
;  kast_wrstrct, kast, /ANONLY, OUTFIL=, FITS=
;
; INPUTS:
;   kast   - An ESI structure
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   ANONLY - Only print files with flg_anly NE 0   
;
; OPTIONAL OUTPUTS:
;   OUTFIL= - Output file (default: kast.list)
;   FITS=   - Name of fits output file
;
; COMMENTS:
;
; EXAMPLES:
;   kast_wrstrct, kast, FITS='kast_13oct02.fits'
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   01-Mar-2003 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro kast_wrstrct, kast, ANONLY=anonly, OUTFIL=outfil, FITS=fits

;
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'kast_wrstrct, kast, /ANONLY, OUTFIL=, FITS= (v1.0)'
      return
  endif 
  
;  Optional Keywords
  
  if not keyword_set( OUTFIL ) then outfil = 'kast.list'

  nimg = n_elements(kast)
  if nimg EQ 1 then begin
      print, 'kast_wrstrct: I dont think you passed a structure in'
      return
  endif
  
  close, /all
  openw, 1, outfil

  for q=0,nimg-1 do begin
      
      if kast[q].obj_id GE 0L then begin
          if kast[q].obj_id LT 10 then id_str = '0'+strtrim(kast[q].obj_id,2) $
          else id_str = strtrim(kast[q].obj_id,2)
      endif else id_str = '-1'
      case kast[q].mode of
          0: mode = 'IMG'
          1: mode = 'SPEC'
          else:
      endcase
      if (not keyword_set( ANONLY ) OR kast[q].flg_anly NE 0) then $
        printf, 1, $
        FORMAT='(i4,1x,a7,1x,i1,1x,a12,1x,a4,1x,a4,1x,a4,a10,i5,f5.2,i2,i2,f5.2,2i2)',$
        q, $
        kast[q].img_root, $
        kast[q].flg_anly, $
        kast[q].Obj, $
        id_str, $
        mode, $
        kast[q].type, $
        kast[q].grising, $
        long(kast[q].exp), $
        kast[q].slit, $
        kast[q].arclamp, $
        kast[q].qtzlamp, $
        kast[q].AM, $
        kast[q].flg_ov, $
        kast[q].flg_final
        
      
  endfor

  close, 1

  if keyword_set( FITS ) then mwrfits, kast, fits, /create


end
