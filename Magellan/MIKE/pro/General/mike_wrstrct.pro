;+ 
; NAME:
; mike_wrstrct   
;    Version 2.0
;
; PURPOSE:
;    Write the mike structure to an ASCII summary file and produce a
;    FITS file (recommended).
;
; CALLING SEQUENCE:
;   
;  mike_wrstrct, mike, /ANONLY, OUTFIL=, FITS=
;
; INPUTS:
;   mike   - An ESI structure
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   /ANONLY - Only print files with flg_anly NE 0   
;
; OPTIONAL OUTPUTS:
;   OUTFIL= - Output file (default: 'mike.list')
;   FITS=   - Name of fits output file (default: 'strct.fits')
;
; COMMENTS:
;
; EXAMPLES:
;   mike_wrstrct, mike, FITS='mike_13oct02.fits'
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

pro mike_wrstrct, mike, ANONLY=anonly, OUTFIL=outfil, FITS=fits

;
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'mike_wrstrct, mike, /ANONLY, OUTFIL=, FITS= (v1.1)'
      return
  endif 
  
;  Optional Keywords
  
  if not keyword_set( OUTFIL ) then outfil = 'mike.list'
  if not keyword_set( FITS ) then fits = 'strct.fits'

  nimg = n_elements(mike)
  if nimg EQ 1 then begin
      print, 'mike_wrstrct: I dont think you passed a structure in'
      return
  endif
  
  close, /all
  openw, 1, outfil

  for q=0,nimg-1 do begin
      
      if mike[q].obj_id GE 0L then begin
          if mike[q].obj_id LT 10 then id_str = '0'+strtrim(mike[q].obj_id,2) $
          else id_str = strtrim(mike[q].obj_id,2)
      endif else id_str = '-1'
;      case mike[q].mode of
;          0: mode = 'IMG'
;          1: mode = 'LWD'
;          2: mode = 'ECH'
;          else:
;      endcase
      if (not keyword_set( ANONLY ) OR mike[q].flg_anly NE 0) then $
        printf, 1, $
        FORMAT='(i4,1x,a7,1x,i1,1x,a12,1x,a4,a5,i5,f5.2,i2,i2,f5.2,2i2)',$
        q, $
        mike[q].img_root, $
        mike[q].flg_anly, $
        mike[q].Obj, $
        id_str, $
;        mode, $
        mike[q].type, $
        long(mike[q].exp), $
        mike[q].slit, $
        mike[q].arclamp, $
        mike[q].qtzlamp, $
        mike[q].AM, $
        mike[q].colbin, $
        mike[q].rowbin
      
  endfor

  close, 1

  ;; Write fits file
  mwrfits, mike, fits, /create


end
