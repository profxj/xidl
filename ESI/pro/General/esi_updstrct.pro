;+ 
; NAME:
; esi_updstrct   
;     Version 1.0
;
; PURPOSE:
;    Creates and outputs a structure for a series of ESI
;    spectroscopic images
;
; CALLING SEQUENCE:
;   
;  esi_updstrct, struct, LIST=list, MKDIR=mkdir, 
;               NOFILE=nofile, NOLIST=nolist
;
; INPUTS:
;
; RETURNS:
;
; OUTPUTS:
;   struct     -  Creates an IDL structure for direct images 
;         -  ASCII file summarizing the structure
;
; OPTIONAL KEYWORDS:
;   LIST       - Image list
;   MKDIR      - Make directories
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   esi_updstrct, nght1_strct, /MKDIR
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   01-Aug-2002 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro esi_updstrct, esi, LIST=list, FITS=fits, EDIT=edit

;
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'esi_updstrct, esi, LIST=, FITS=, /EDIT (v1.0)'
      return
  endif 
  
;  Optional Keywords
  if not keyword_set( FRMCRD ) then frmcrd = 'FRAMENO'
  
;  Grab the files
  img = findfile('Raw/esi*.fits*') 
  nimg = n_elements(img)

;  Loop on Indv images
  flg_fst = 0
  for q=0,nimg-1 do begin
      head = xheadfits(img[q], /silent)
      frameno = sxpar(head, frmcrd)
      a = where(esi.frame EQ frameno, na)
      if na EQ 0 then begin
          if flg_fst EQ 0 then begin
              newlst = [img[q]]
              flg_fst = 1
          endif else newlst = [newlst, img[q]]
      endif
  endfor
          
; Run esi_strct

  if flg_fst EQ 1 then esi_strct, newesi, IMG=newlst, /NOFILE $
  else begin
      print, 'esi_updstrct: No new files!'
      return
  endelse

; Update

  esi = [esi, newesi]
  srt = sort(esi.frame)
  esi = temporary(esi[srt])

; Output

  esi_wrstrct, esi, FITS=fits

  return
end
